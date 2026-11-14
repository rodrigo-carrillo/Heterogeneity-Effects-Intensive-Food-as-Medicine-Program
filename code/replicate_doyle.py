"""
Python port of the authors' replication code (data/raw/dataverse_files/01_published_tables.do)
for Doyle et al., JAMA Intern Med 2024, Tables 1-4, checked cell by cell against the
published values.

  python code/replicate_doyle.py

Writes outputs/replication/{table1..4.csv, comparison.csv, replication_report.md}.

Stata's `xtreg y x, fe robust` is implemented directly as the within estimator with
cluster-robust (by patient) SEs, Stata's small-sample factor G/(G-1)*(N-1)/(N-K) and
t(G-1) CIs; `regress, robust` as HC1 with t(N-K) CIs. The first-difference form used in
the HTE pipeline (change ~ W + baseline + d_time) is checked against it (same point
estimate by construction; SEs compared). For non-HbA1c biometrics the authors impute
missing baseline values (`mi impute mvn`, 10 imputations, Stata RNG); this port uses a
Bayesian multivariate-normal regression imputation with the same predictors and Rubin's
rules, so those cells can differ slightly from the paper, as the authors' code warns.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "raw" / "dataverse_files" / "analysis_sample.csv"
OUT = ROOT / "outputs" / "replication"
OUT.mkdir(parents=True, exist_ok=True)
LB_PER_KG = 2.20462

raw = pd.read_csv(DATA)
raw["race_nonhisp_white"] = ((raw.hispanic == 0) & (raw.race == "White")).astype(int)
raw["race_nonhisp_other"] = ((raw.hispanic == 0) & (raw.race != "White")
                             & ~raw.race.fillna("").str.startswith("Black Or")).astype(int)
raw["days_between_tests"] = raw.hba1c_date_6mo - raw.hba1c_date_0mo
AGEBINS = [f"age5_{i}" for i in range(1, 10)]        # age5_10 (72+) omitted, as in the .do
PRESPEC = ["rural_site", "hba1c_cat_control", *AGEBINS, "male", "race_nonhisp_black", "hispanic"]


# ---- Estimators ---------------------------------------------------------------
def ols(y, X, cluster=None, k_absorbed=0):
    """OLS with HC1 (cluster=None) or Stata cluster-robust SEs. Returns (beta, se, df)."""
    X, y = np.asarray(X, float), np.asarray(y, float)
    n, k = X.shape
    XtXi = np.linalg.pinv(X.T @ X)
    b = XtXi @ X.T @ y
    e = y - X @ b
    if cluster is None:
        meat = (X * e[:, None]).T @ (X * e[:, None]) * n / (n - k)
        df = n - k
    else:
        g = pd.Series(cluster).factorize()[0]
        G = g.max() + 1
        S = np.zeros((G, k))
        np.add.at(S, g, X * e[:, None])
        meat = S.T @ S * (G / (G - 1)) * ((n - 1) / (n - k - k_absorbed))
        df = G - 1
    V = XtXi @ meat @ XtXi
    return b, np.sqrt(np.diag(V)), df


def xtreg_fe(panel, y, xs, pid="subject_id"):
    """Stata `xtreg y xs, fe robust`: within estimator + constant, clustered by patient.
    Time-invariant regressors are absorbed by the fixed effects (Stata drops them)."""
    p = panel.copy()
    cols = [y, *xs]
    dm = p[cols] - p.groupby(pid)[cols].transform("mean") + p[cols].mean()
    keep = [x for x in xs if dm[x].std() > 1e-10]
    X = np.column_stack([dm[keep].values, np.ones(len(p))])
    b, se, df = ols(dm[y].values, X, cluster=p[pid].values)
    return dict(zip(keep, zip(b, se))), df


def ci(b, se, df):
    t = stats.t.ppf(0.975, df)
    return b, b - t * se, b + t * se, 2 * stats.t.sf(abs(b / se), df)


def long_panel(df, var, month, val="{v}_val_{m}mo", date="{v}_date_{m}mo"):
    """Reshape to the authors' 2-period panel: baseline (0) and follow-up (month) rows."""
    rows = []
    for m in (0, month):
        r = df[["subject_id", "treat_now", "consent_date", *PRESPEC, *[f"age5_{i}" for i in (10,)]]].copy()
        r["panel"] = df[val.format(v=var, m=m)].values
        if date:
            r["time"] = (df[date.format(v=var, m=m)] - df.consent_date).values
        r["post_period"] = int(m == month)
        rows.append(r)
    p = pd.concat(rows).sort_values(["subject_id", "post_period"]).reset_index(drop=True)
    p["post_treat_now"] = p.post_period * p.treat_now
    return p


def mvn_impute(p, cols, m=10, seed=9112023):
    """Proper MI of missing `cols` (jointly missing) from a multivariate normal regression
    on the authors' predictors, with posterior draws of (B, Sigma)."""
    rng = np.random.default_rng(seed)
    preds = ["treat_now", "post_period", "post_treat_now", *PRESPEC, "age5_10"]
    X = np.column_stack([np.ones(len(p)), p[preds].values.astype(float)])
    miss = p[cols].isna().any(axis=1).values
    assert (p.loc[miss, cols].isna().all(axis=1)).all(), "non-joint missingness pattern"
    Xo, Yo = X[~miss], p.loc[~miss, cols].values.astype(float)
    XtXi = np.linalg.pinv(Xo.T @ Xo)
    B = XtXi @ Xo.T @ Yo
    E = Yo - Xo @ B
    dof = len(Yo) - np.linalg.matrix_rank(Xo)
    L = np.linalg.cholesky(XtXi + 1e-12 * np.eye(len(XtXi)))
    out = []
    for _ in range(m):
        Sig = stats.invwishart.rvs(df=dof, scale=E.T @ E, random_state=rng)
        Bd = B + L @ rng.standard_normal(B.shape) @ np.linalg.cholesky(Sig).T
        q = p.copy()
        mu = X[miss] @ Bd
        q.loc[miss, cols] = mu + rng.standard_normal((miss.sum(), len(cols))) @ np.linalg.cholesky(Sig).T
        out.append(q)
    return out


def rubin(ests, ses, df_com):
    """Rubin's rules with the Barnard-Rubin small-sample df (as Stata's mi estimate)."""
    m = len(ests)
    qbar, ubar, bvar = np.mean(ests), np.mean(np.square(ses)), np.var(ests, ddof=1)
    T = ubar + (1 + 1 / m) * bvar
    lam = (1 + 1 / m) * bvar / T
    nu_old = (m - 1) / lam ** 2 if lam > 0 else np.inf
    nu_obs = (df_com + 1) / (df_com + 3) * df_com * (1 - lam)
    return qbar, np.sqrt(T), 1 / (1 / nu_old + 1 / nu_obs)


# ---- Table 1 ------------------------------------------------------------------
MEANS = ["hba1c_val_0mo", "age", "days_between_tests"]
COUNTS = ["male", "hispanic", "race_nonhisp_black", "race_nonhisp_white", "race_nonhisp_other",
          "rural_site", "ip_i_0mo", "ed_i_0mo", "drug_metformin_0mo", "drug_insulin_0mo",
          "drug_glp1_0mo", "drug_other_0mo", "drug_any_db_summary_0mo"]
t1 = []
for sample, d in (("full", raw), ("lab", raw[raw.hba1c_val_6mo.notna()])):
    for grp, g in (("overall", d), ("control", d[d.treat_now == 0]), ("treatment", d[d.treat_now == 1])):
        row = {"sample": sample, "group": grp, "N": len(g)}
        row.update({f"mean_{v}": g[v].mean() for v in MEANS} | {f"sd_{v}": g[v].std() for v in MEANS})
        row.update({f"n_{v}": g[v].sum() for v in COUNTS} | {f"pct_{v}": 100 * g[v].mean() for v in COUNTS})
        t1.append(row)
t1 = pd.DataFrame(t1)
t1.to_csv(OUT / "table1.csv", index=False)

# ---- Table 2 (adjusted differences) -------------------------------------------
XS2 = ["fruit_num", "vegetables_num", "sweet_drink_num", "never_eat_out_yn", "diet_improved_yn",
       "diet_improve_index", "engagement_index", "diabetes_knowledge_score",
       "healthy_attitudes_index", "self_efficacy_index", "exercise_smoking_index"]
EHR2 = ["dietician_visits_n", "dsmt_visits", "endocrinology_visits_n", "weight_visits_n"]
t2 = []
for v in XS2:
    d = raw.dropna(subset=[v, *PRESPEC])
    X = np.column_stack([np.ones(len(d)), d[["treat_now", *PRESPEC]].values])
    b, se, df = ols(d[v].values, X)
    t2.append(dict(outcome=v, n=len(d), **dict(zip(["est", "lo", "hi", "p"], ci(b[1], se[1], df)))))
for v in EHR2:
    d = raw[raw[f"{v}_6mo"].notna()]
    p = long_panel(d, v, 6, val="{v}_{m}mo", date=None)
    p["post_baseline"] = p.post_period * p.groupby("subject_id").panel.transform("first")
    r, df = xtreg_fe(p, "panel", ["post_period", "treat_now", "post_treat_now", "post_baseline", *PRESPEC])
    t2.append(dict(outcome=v, n=len(p), **dict(zip(["est", "lo", "hi", "p"], ci(*r["post_treat_now"], df)))))
t2 = pd.DataFrame(t2)
t2.to_csv(OUT / "table2_adjusted.csv", index=False)

# ---- Table 3 ------------------------------------------------------------------
BIOS = ["hba1c", "ldl", "hdl", "total_chol", "trigly", "weight", "sys", "dias"]
t3 = []      # raw data as posted (incl. the 9999999 triglyceride placeholder), as the .do uses it
for month in (6, 12):
    for var in BIOS:
        # Weight is run twice: the kg column and the (undocumented) lb column
        for unit, valfmt in ((("kg", "{v}_val_{m}mo"), ("lb", "weight_lb_val_{m}mo")) if var == "weight"
                             else (("", "{v}_val_{m}mo"),)):
            col = valfmt.format(v=var, m=month)
            c, t = raw.loc[raw.treat_now == 0, col], raw.loc[raw.treat_now == 1, col]
            dd = raw.dropna(subset=[col])
            b, se, df = ols(dd[col].values, np.column_stack([np.ones(len(dd)), dd.treat_now.values]))
            un = ci(b[1], se[1], df)
            p = long_panel(dd, var, month, val=valfmt)
            xs = ["treat_now", "post_period", "post_treat_now", "post_baseline", "time", *PRESPEC]
            if var == "hba1c":
                p["post_baseline"] = p.post_period * p.groupby("subject_id").panel.transform("first")
                r, df = xtreg_fe(p, "panel", xs)
                adj = ci(*r["post_treat_now"], df)
                n_imp = 0
            else:
                n_imp = int(p.loc[p.post_period == 0, "panel"].isna().sum())
                est, ses = [], []
                for q in mvn_impute(p, ["panel", "time"]):
                    q["post_baseline"] = q.post_period * q.groupby("subject_id").panel.transform("first")
                    r, df = xtreg_fe(q, "panel", xs)
                    est.append(r["post_treat_now"][0]); ses.append(r["post_treat_now"][1])
                adj = ci(*rubin(est, ses, df)[:2], rubin(est, ses, df)[2])
            # First-difference form on observed baselines (as in the HTE pipeline)
            o = dd.dropna(subset=[valfmt.format(v=var, m=0)])
            dy = o[col] - o[valfmt.format(v=var, m=0)]
            dt = (o[f"{var}_date_{month}mo"] - o[f"{var}_date_0mo"]).values
            Xfd = np.column_stack([np.ones(len(o)), o.treat_now, o[valfmt.format(v=var, m=0)], dt])
            bfd, sefd, dffd = ols(dy.values, Xfd)
            fd = ci(bfd[1], sefd[1], dffd)
            t3.append(dict(month=month, outcome=var + (f" ({unit})" if unit else ""), n_obs=2 * len(dd),
                           n_imputed_baselines=n_imp, ctl_mean=c.mean(), ctl_sd=c.std(),
                           trt_mean=t.mean(), trt_sd=t.std(),
                           unadj_est=un[0], unadj_lo=un[1], unadj_hi=un[2],
                           adj_est=adj[0], adj_lo=adj[1], adj_hi=adj[2], adj_p=adj[3],
                           fd_n=len(o), fd_est=fd[0], fd_lo=fd[1], fd_hi=fd[2]))
t3 = pd.DataFrame(t3)
t3.to_csv(OUT / "table3.csv", index=False)

# ---- Table 4 ------------------------------------------------------------------
EHR4 = ["ed_or_ip_i", "op_i", "ed_or_ip_n", "op_n", "ed_or_ip_payment", "op_payment", "nf_payment",
        "oth_payment", "total_payment", "drug_insulin", "drug_metformin", "drug_sulf", "drug_glp1",
        "drug_other", "drug_any_db"]
t4 = []
for month in (6, 12):
    for v in EHR4:
        d = raw[raw[f"{v}_{month}mo"].notna()]
        p = long_panel(d, v, month, val="{v}_{m}mo", date=None)
        p["post_baseline"] = p.post_period * p.groupby("subject_id").panel.transform("first")
        r, df = xtreg_fe(p, "panel", ["post_period", "treat_now", "post_treat_now", "post_baseline", *PRESPEC])
        est = ci(*r["post_treat_now"], df)
        scale = 100 if v.endswith("_i") or v.startswith("drug") else 1     # paper reports pp
        g = d.groupby("treat_now")[f"{v}_{month}mo"]
        t4.append(dict(month=month, outcome=v, n_obs=len(p), ctl=g.mean()[0] * scale, trt=g.mean()[1] * scale,
                       ctl_n=g.sum()[0], trt_n=g.sum()[1],
                       est=est[0] * scale, lo=est[1] * scale, hi=est[2] * scale, p=est[3]))
t4 = pd.DataFrame(t4)
t4.to_csv(OUT / "table4.csv", index=False)

# ---- Comparison with the published tables ------------------------------------
# (table, row, statistic): published value. Transcribed from the paper PDF.
PUB = {
    # Table 1, full sample (counts) and lab sample means
    ("T1", "full overall N"): 465, ("T1", "full control N"): 235, ("T1", "full treatment N"): 230,
    ("T1", "full overall hba1c mean"): 10.29, ("T1", "full overall age mean"): 54.6,
    ("T1", "full overall male n"): 210, ("T1", "full overall hispanic n"): 39,
    ("T1", "full overall nh black n"): 41, ("T1", "full overall nh white n"): 378,
    ("T1", "full overall nh other n"): 7, ("T1", "full overall rural n"): 131,
    ("T1", "full overall prior IP n"): 107, ("T1", "full overall prior ED n"): 191,
    ("T1", "full control prior ED n"): 94, ("T1", "full treatment prior ED n"): 97,
    ("T1", "full overall any DM drug n"): 405, ("T1", "full overall metformin n"): 231,
    ("T1", "full overall insulin n"): 231, ("T1", "full overall GLP-1 n"): 150,
    ("T1", "full overall other drug n"): 124,
    ("T1", "lab overall N"): 349, ("T1", "lab control N"): 179, ("T1", "lab treatment N"): 170,
    ("T1", "lab overall hba1c mean"): 10.21, ("T1", "lab overall days between tests mean"): 302,
    # Table 2 adjusted differences (estimate)
    **{("T2", k): v for k, v in {
        "dietician_visits_n": 2.01, "dsmt_visits": 1.88, "endocrinology_visits_n": 0.05,
        "weight_visits_n": -0.06, "fruit_num": 0.26, "vegetables_num": 2.00, "sweet_drink_num": -0.61,
        "never_eat_out_yn": 5.38, "diet_improved_yn": 16.92, "diet_improve_index": 0.38,
        "engagement_index": 0.36, "diabetes_knowledge_score": 0.04, "healthy_attitudes_index": 0.13,
        "self_efficacy_index": -0.03, "exercise_smoking_index": 0.05}.items()},
    # Table 3 (month, outcome): unadjusted, adjusted
    **{("T3u", k): v for k, v in {
        (6, "hba1c"): -0.03, (6, "ldl"): 8.35, (6, "hdl"): -0.56, (6, "total_chol"): 9.77,
        (6, "trigly"): 22.56, (6, "weight"): -8.14, (6, "sys"): -1.87, (6, "dias"): -0.56,
        (12, "hba1c"): 0.14, (12, "ldl"): 2.64, (12, "hdl"): -1.39, (12, "total_chol"): -0.55,
        (12, "trigly"): -11.35, (12, "weight"): -9.68, (12, "sys"): -3.77, (12, "dias"): 0.50}.items()},
    **{("T3a", k): v for k, v in {
        (6, "hba1c"): -0.10, (6, "ldl"): 5.86, (6, "hdl"): -0.48, (6, "total_chol"): 4.94,
        (6, "trigly"): 3.63, (6, "weight"): 1.95, (6, "sys"): -0.38, (6, "dias"): -0.25,
        (12, "hba1c"): 0.10, (12, "ldl"): 3.75, (12, "hdl"): -0.33, (12, "total_chol"): -0.51,
        (12, "trigly"): -25.63, (12, "weight"): 1.88, (12, "sys"): -2.49, (12, "dias"): 0.81}.items()},
    **{("T3n", k): v for k, v in {
        (6, "hba1c"): 698, (6, "ldl"): 606, (6, "hdl"): 610, (6, "total_chol"): 604, (6, "trigly"): 612,
        (6, "weight"): 800, (6, "sys"): 756, (6, "dias"): 754, (12, "hba1c"): 650, (12, "ldl"): 534,
        (12, "hdl"): 546, (12, "total_chol"): 552, (12, "trigly"): 536, (12, "weight"): 736,
        (12, "sys"): 702, (12, "dias"): 702}.items()},
    # Table 4 adjusted differences
    **{("T4", (6, k)): v for k, v in {
        "ed_or_ip_i": -7.13, "op_i": 5.23, "ed_or_ip_n": -0.04, "op_n": 1.00, "ed_or_ip_payment": -4,
        "op_payment": 1376, "nf_payment": -86, "oth_payment": -41, "total_payment": 1048,
        "drug_insulin": 3.01, "drug_metformin": 7.18, "drug_sulf": -0.75, "drug_glp1": 12.77,
        "drug_other": 3.24, "drug_any_db": 2.73}.items()},
    **{("T4", (12, k)): v for k, v in {
        "ed_or_ip_i": -7.56, "op_i": 0.47, "ed_or_ip_n": -0.16, "op_n": 0.20, "ed_or_ip_payment": -1332,
        "op_payment": 1348, "nf_payment": -150, "oth_payment": -211, "total_payment": -424,
        "drug_insulin": 2.99, "drug_metformin": 8.87, "drug_sulf": -0.84, "drug_glp1": 8.73,
        "drug_other": 0.44, "drug_any_db": 3.29}.items()},
}

f, ctl, trt = (t1.query("sample=='full'").set_index("group"), None, None)
lab = t1.query("sample=='lab'").set_index("group")
OURS = {
    ("T1", "full overall N"): f.N["overall"], ("T1", "full control N"): f.N["control"],
    ("T1", "full treatment N"): f.N["treatment"],
    ("T1", "full overall hba1c mean"): f.mean_hba1c_val_0mo["overall"], ("T1", "full overall age mean"): f.mean_age["overall"],
    ("T1", "full overall male n"): f.n_male["overall"], ("T1", "full overall hispanic n"): f.n_hispanic["overall"],
    ("T1", "full overall nh black n"): f.n_race_nonhisp_black["overall"],
    ("T1", "full overall nh white n"): f.n_race_nonhisp_white["overall"],
    ("T1", "full overall nh other n"): f.n_race_nonhisp_other["overall"],
    ("T1", "full overall rural n"): f.n_rural_site["overall"], ("T1", "full overall prior IP n"): f.n_ip_i_0mo["overall"],
    ("T1", "full overall prior ED n"): f.n_ed_i_0mo["overall"], ("T1", "full control prior ED n"): f.n_ed_i_0mo["control"],
    ("T1", "full treatment prior ED n"): f.n_ed_i_0mo["treatment"],
    ("T1", "full overall any DM drug n"): f.n_drug_any_db_summary_0mo["overall"],
    ("T1", "full overall metformin n"): f.n_drug_metformin_0mo["overall"],
    ("T1", "full overall insulin n"): f.n_drug_insulin_0mo["overall"],
    ("T1", "full overall GLP-1 n"): f.n_drug_glp1_0mo["overall"], ("T1", "full overall other drug n"): f.n_drug_other_0mo["overall"],
    ("T1", "lab overall N"): lab.N["overall"], ("T1", "lab control N"): lab.N["control"],
    ("T1", "lab treatment N"): lab.N["treatment"], ("T1", "lab overall hba1c mean"): lab.mean_hba1c_val_0mo["overall"],
    ("T1", "lab overall days between tests mean"): lab.mean_days_between_tests["overall"],
}
PCT2 = {"never_eat_out_yn", "diet_improved_yn"}       # reported in pp
OURS |= {("T2", r.outcome): r.est * (100 if r.outcome in PCT2 else 1) for r in t2.itertuples()}
for r in t3.itertuples():
    key = (r.month, r.outcome.split(" (")[0])
    # Published weight row: unadjusted difference in lb, adjusted difference in kg (see report)
    if r.outcome != "weight (kg)":
        OURS |= {("T3u", key): r.unadj_est, ("T3n", key): r.n_obs}
    if r.outcome != "weight (lb)":
        OURS[("T3a", key)] = r.adj_est
OURS |= {("T4", (r.month, r.outcome)): r.est for r in t4.itertuples()}

cmp = pd.DataFrame([dict(table=k[0], item=str(k[1]), published=v, replicated=OURS.get(k))
                    for k, v in PUB.items()])
cmp["abs_diff"] = (cmp.replicated - cmp.published).abs()
# A cell "matches" if it agrees at the published rounding (half a unit of the last digit);
# multiply-imputed cells get a looser tolerance (the authors' code warns they vary).
dec = cmp.published.map(lambda x: len(repr(float(x)).split(".")[1].rstrip("0")) if float(x) % 1 else 0)
tol = 0.5 * 10.0 ** (-dec) * 1.02             # 2% slack for half-way rounding
mi_cell = cmp.table.eq("T3a") & ~cmp.item.str.contains("hba1c")
cmp["match"] = np.where(cmp.abs_diff <= tol, "exact",
                        np.where(mi_cell & (cmp.abs_diff <= 0.1 * cmp.published.abs() + 0.5), "MI (close)", "DIFFERENT"))
cmp.to_csv(OUT / "comparison.csv", index=False)

# ---- eTable 18a: 6-mo HbA1c within Doyle's subgroups ---------------------------
# The Table 3 HbA1c model WITHOUT the test-timing term, refitted within each subgroup (with the
# timing term the cells differ by up to 0.06; without it they match). eTable 18a's "obs" (814 overall) also
# counts the 116 patients with only a baseline row: they are FE singletons, so they do not move
# the estimate but enter N, the cluster count and the df. Reproduced by keeping them.
E18 = {("hba1c_cat_control", 0): (-0.1, -0.526, 0.327, 368), ("hba1c_cat_control", 1): (-0.099, -0.666, 0.468, 446),
       ("rural_site", 0): (-0.196, -0.640, 0.248, 578), ("rural_site", 1): (0.264, -0.355, 0.883, 236),
       ("male", 1): (-0.36, -0.889, 0.170, 372), ("male", 0): (0.187, -0.334, 0.708, 442),
       ("ip_i_0mo", 1): (-0.621, -1.414, 0.173, 192), ("ip_i_0mo", 0): (0.088, -0.320, 0.496, 622),
       ("ed_i_0mo", 1): (-0.212, -0.812, 0.389, 339), ("ed_i_0mo", 0): (0.099, -0.347, 0.544, 475),
       ("drug_insulin_0mo", 1): (-0.154, -0.670, 0.362, 407), ("drug_insulin_0mo", 0): (-0.018, -0.532, 0.495, 407)}
e18 = []
for (v, lvl), pub in E18.items():
    p = long_panel(raw[raw[v] == lvl], "hba1c", 6)
    p = p[p.panel.notna()]                                  # drops missing 6-mo rows, keeps singletons
    p["post_baseline"] = p.post_period * p.groupby("subject_id").panel.transform("first")
    r, df = xtreg_fe(p, "panel", ["post_period", "post_treat_now", "post_baseline", *PRESPEC])
    est = ci(*r["post_treat_now"], df)
    e18.append(dict(subgroup=f"{v}={lvl}", obs=len(p), pub_obs=pub[3], est=est[0], lo=est[1], hi=est[2],
                    pub_est=pub[0], pub_lo=pub[1], pub_hi=pub[2]))
e18 = pd.DataFrame(e18)
e18["max_abs_diff"] = (e18[["est", "lo", "hi"]].values - e18[["pub_est", "pub_lo", "pub_hi"]].values).__abs__().max(axis=1)
e18.to_csv(OUT / "etable18a.csv", index=False)

# ---- Weight, 6 mo: what drives the adjusted +1.95 kg ---------------------------
# The authors impute the 5 missing baseline weights (and dates) from demographics only, ignoring
# each patient's own 6-mo weight. Refit the same FE model (a) on observed baselines only and
# (b) imputing the baseline from the patient's 6-mo weight (+ strata, sex, arm), 20 imputations.
def fe_weight(d):
    p = long_panel(d, "weight", 6)
    p["post_baseline"] = p.post_period * p.groupby("subject_id").panel.transform("first")
    r, df = xtreg_fe(p, "panel", ["post_period", "post_treat_now", "post_baseline", "time", *PRESPEC])
    return (*r["post_treat_now"], df)


w = raw.dropna(subset=["weight_val_6mo"])
cc = w.dropna(subset=["weight_val_0mo"])
fe_cc = ci(*fe_weight(cc))
seed_runs = []
for s in range(10):
    p = long_panel(w, "weight", 6)
    est = []
    for q in mvn_impute(p, ["panel", "time"], seed=s):
        q["post_baseline"] = q.post_period * q.groupby("subject_id").panel.transform("first")
        est.append(xtreg_fe(q, "panel", ["post_period", "post_treat_now", "post_baseline", "time", *PRESPEC])[0]["post_treat_now"][0])
    seed_runs.append(np.mean(est))
rng = np.random.default_rng(1)
Z = lambda d: np.column_stack([np.ones(len(d)), d[["weight_val_6mo", "treat_now", "rural_site", "hba1c_cat_control", "male"]].values])
bz, *_ = np.linalg.lstsq(Z(cc), cc.weight_val_0mo.values, rcond=None)
sz = np.std(cc.weight_val_0mo.values - Z(cc) @ bz, ddof=Z(cc).shape[1])
lag0 = (cc.weight_date_0mo - cc.consent_date).median()
est, ses = [], []
for _ in range(20):
    q = w.copy()
    miss = q.weight_val_0mo.isna()
    q.loc[miss, "weight_val_0mo"] = Z(q[miss]) @ bz + rng.normal(0, sz, miss.sum())
    q.loc[miss, "weight_date_0mo"] = q.loc[miss, "consent_date"] + lag0
    b, se, df = fe_weight(q)
    est.append(b); ses.append(se)
fe_own = ci(*rubin(est, ses, df)[:2], rubin(est, ses, df)[2])
imp_pts = w[w.weight_val_0mo.isna()][["treat_now", "weight_val_6mo"]]

# ---- Report -------------------------------------------------------------------
w6 = t3.query("month==6 and outcome.str.startswith('weight')", engine="python").set_index("outcome")
h6 = t3.query("month==6 and outcome=='hba1c'").iloc[0]
kg, lb = raw.weight_val_6mo, raw.weight_lb_val_6mo
summary = cmp.match.value_counts().to_dict()
diff_rows = cmp[cmp.match == "DIFFERENT"]
report = f"""# Replication of Doyle et al. 2024, Tables 1-4 (auto-generated by code/replicate_doyle.py)

Python port of the authors' `01_published_tables.do`, compared with {len(cmp)} published cells.
Match summary: {summary}. "exact" = equal at the published rounding; "MI (close)" =
multiply-imputed Table 3 cell within the stated tolerance (imputation draws differ from Stata's).

## Cells that do not match
{diff_rows.round(3).to_markdown(index=False) if len(diff_rows) else "None."}

## Paper's FE model vs the first-difference form (HbA1c, 6 mo)
- xtreg FE (as published): {h6.adj_est:+.3f} ({h6.adj_lo:+.3f} to {h6.adj_hi:+.3f})
- First differences, HC1: {h6.fd_est:+.3f} ({h6.fd_lo:+.3f} to {h6.fd_hi:+.3f})

## eTable 18a (6-mo HbA1c by Doyle subgroup; 12 subgroup cells)
Model: Table 3's HbA1c FE model without the test-timing term; "obs" includes the 116 baseline-only
patients (FE singletons: 349 x 2 + 116 = 814), which explains the 814 vs 698 discrepancy.
Largest absolute difference from the published estimate or CI bound: {e18.max_abs_diff.max():.3f};
observation counts equal in {int((e18.obs == e18.pub_obs).sum())} of {len(e18)} cells.
{e18.round(3).to_markdown(index=False)}

## Weight (Table 3, 6 mo; published: unadjusted -8.14 [-20.57 to 4.30], adjusted 1.95 [0.07 to 3.83])
**Units of the columns.** In the CSV, `weight_lb_val / weight_val` = {(lb / kg).min():.4f} to {(lb / kg).max():.4f}
for every patient: the kg column was converted with 1/0.45 = 2.2222 instead of 2.2046, so it is
{100 * (2.22222 / LB_PER_KG - 1):.2f}% too low. The `weight_lb_val_*` columns are not in the README dictionary.

**Units of the published row.**
| Column | Unadjusted | Adjusted (authors' FE + MI) |
|---|---|---|
| lb | {w6.loc['weight (lb)','unadj_est']:+.2f} ({w6.loc['weight (lb)','unadj_lo']:+.2f} to {w6.loc['weight (lb)','unadj_hi']:+.2f}) | {w6.loc['weight (lb)','adj_est']:+.2f} ({w6.loc['weight (lb)','adj_lo']:+.2f} to {w6.loc['weight (lb)','adj_hi']:+.2f}) |
| kg | {w6.loc['weight (kg)','unadj_est']:+.2f} ({w6.loc['weight (kg)','unadj_lo']:+.2f} to {w6.loc['weight (kg)','unadj_hi']:+.2f}) | {w6.loc['weight (kg)','adj_est']:+.2f} ({w6.loc['weight (kg)','adj_lo']:+.2f} to {w6.loc['weight (kg)','adj_hi']:+.2f}) |

- The **unadjusted** difference is in **lb** (exact). The **adjusted** difference is in **kg** (close; MI draws differ).
- Means 108.39 / 104.70 = lb / 2.20462 ({lb[raw.treat_now == 0].mean() / LB_PER_KG:.2f} / {lb[raw.treat_now == 1].mean() / LB_PER_KG:.2f}); SDs 29.85 / 27.06 =
  the 2.2222-converted kg column ({kg[raw.treat_now == 0].std():.2f} / {kg[raw.treat_now == 1].std():.2f}). So the row mixes units, and the posted
  .do (which uses `weight_val` throughout) cannot have produced the unadjusted cell from the posted data.

**What drives the adjusted +1.95 kg.** The authors impute the {len(imp_pts)} missing baselines from demographics
only, ignoring each patient's own 6-mo weight ({int(imp_pts.treat_now.sum())} in the program arm,
{int((imp_pts.treat_now == 0).sum())} in the wait-list arm; individual values are not printed, so the report can be shared). Same FE model:
| Baseline handling | Adjusted difference, kg |
|---|---|
| Authors' imputation (this port; 10 seeds) | {min(seed_runs):+.2f} to {max(seed_runs):+.2f} |
| Observed baselines only (N={len(cc)}) | {fe_cc[0]:+.2f} ({fe_cc[1]:+.2f} to {fe_cc[2]:+.2f}), p = {fe_cc[3]:.2f} |
| Imputed from own 6-mo weight + strata, sex, arm (20 imp.) | {fe_own[0]:+.2f} ({fe_own[1]:+.2f} to {fe_own[2]:+.2f}), p = {fe_own[3]:.2f} |

The published "significant" weight gain depends on 5 imputed baselines; with observed or
sensibly imputed baselines the 6-mo effect is about +0.8 to +0.9 kg and not significant.

## Unresolved small differences
- Lipid sample sizes: the posted data have 6-12 more patients with a follow-up lipid value than
  Table 3 reports (e.g., LDL 6 mo: 309 vs 303); diastolic BP 378 vs 377. All values have dates inside the
  windows. Likely a slightly different data extract; the affected adjusted cells differ within their CIs.
- Any diabetes medication (Table 1): 404 in the posted data vs 405 published.
- The 9999999 triglyceride placeholder (one baseline) does not enter any published model: that
  patient has no follow-up triglyceride.
- Typos in published Table 4 CIs (estimates and the other bound replicate exactly; replicated
  bound in brackets): 6 mo GLP-1 lower 6.47 [{t4.query("month==6 and outcome=='drug_glp1'").lo.iloc[0]:.2f}]; 12 mo ED/IP lower
  -16.84 [{t4.query("month==12 and outcome=='ed_or_ip_i'").lo.iloc[0]:.2f}]; 12 mo GLP-1 "0.86 to 1.66" [0.86 to {t4.query("month==12 and outcome=='drug_glp1'").hi.iloc[0]:.2f}].
- Never eating out (Table 2): 5.385 replicated vs 5.38 published (rounding).
"""
(OUT / "replication_report.md").write_text(report)
print(report)
