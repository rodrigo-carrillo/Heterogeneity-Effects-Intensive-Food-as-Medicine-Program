"""Claim 3 (weight x prior ED) and the subgroup part of claim 4 (HbA1c), independent implementation.
Methods: (a) OLS change ~ T*S (+ strata) HC1; (b) stratified difference in means within subgroup
(stratum-propensity IPW / Hajek); (c) cross-fitted AIPW with known stratum propensity and
linear outcome models fit separately by arm; (d) within-stratum permutation test of interaction."""
import numpy as np, pandas as pd
from scipy import stats

rng = np.random.default_rng(20261021)
D = pd.read_csv("data/raw/dataverse_files/analysis_sample.csv")
D["stratum"] = D.rural_site * 2 + D.hba1c_cat_control
D["days_a1c_consent"] = D.consent_date - D.hba1c_date_0mo
SUBS = ["hba1c_cat_control", "rural_site", "male", "ip_i_0mo", "ed_i_0mo",
        "drug_insulin_0mo", "drug_metformin_0mo"]
COV = ["age", "hba1c_val_0mo", "days_a1c_consent", "bmi_val_0mo", "sys_val_0mo", "dias_val_0mo",
       "op_n_0mo", "ed_or_ip_n_0mo", "male", "hispanic", "race_nonhisp_black", "rural_site",
       "hba1c_cat_control", "drug_insulin_0mo", "drug_metformin_0mo", "drug_glp1_0mo",
       "drug_sulf_0mo", "drug_other_0mo", "drug_any_db_0mo", "ip_i_0mo", "ed_i_0mo"]


def hc1(X, y):
    A = np.linalg.inv(X.T @ X); b = A @ X.T @ y; e = y - X @ b; n, k = X.shape
    return b, A @ (X.T * e**2) @ X @ A * n / (n - k)


def ols_int(d, s, y, strata=True):
    cols = [np.ones(len(d)), d.treat_now, d[s], d.treat_now * d[s]]
    if strata:
        for k in [1, 2, 3]:
            cols.append((d.stratum == k).astype(float))
    X = np.column_stack(cols); X = X[:, _indep(X)]
    b, V = hc1(X, d[y].values)
    return b[3], np.sqrt(V[3, 3]), b[1], b[1] + b[3]


def _indep(X):
    keep = []
    for j in range(X.shape[1]):
        if np.linalg.matrix_rank(X[:, keep + [j]]) == len(keep) + 1: keep.append(j)
    return keep


def stratum_prop(d):
    return d.groupby("stratum").treat_now.transform("mean").values


def aipw_scores(d, y, covs, K=5, reps=1):
    """Cross-fitted AIPW scores; propensity = known stratum share (sample share in stratum)."""
    e = stratum_prop(d); T = d.treat_now.values; Y = d[y].values
    X = d[covs].values.astype(float)
    X = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1)
    out = np.zeros(len(d))
    for r in range(reps):
        folds = rng.permutation(np.arange(len(d)) % K)
        mu1 = np.zeros(len(d)); mu0 = np.zeros(len(d))
        for k in range(K):
            tr, te = folds != k, folds == k
            for t, mu in [(1, mu1), (0, mu0)]:
                m = tr & (T == t)
                Z = np.column_stack([np.ones(m.sum()), X[m]])
                # light ridge for stability
                lam = 1.0; P = np.eye(Z.shape[1]); P[0, 0] = 0
                b = np.linalg.solve(Z.T @ Z + lam * P, Z.T @ Y[m])
                mu[te] = np.column_stack([np.ones(te.sum()), X[te]]) @ b
        out += mu1 - mu0 + T * (Y - mu1) / e - (1 - T) * (Y - mu0) / (1 - e)
    return out / reps


def sub_from_scores(g, s):
    res = {}
    for v in [0, 1]:
        x = g[s == v]; res[v] = (x.mean(), x.std(ddof=1) / np.sqrt(len(x)), len(x))
    diff = res[1][0] - res[0][0]; se = np.hypot(res[1][1], res[0][1])
    return res, diff, se, 2 * stats.norm.sf(abs(diff / se))


def perm_test(d, s, y, B=10000):
    obs = ols_int(d, s, y)[0]
    T = d.treat_now.values.copy(); strata = d.stratum.values
    idx = [np.where(strata == k)[0] for k in np.unique(strata)]
    Y = d[y].values; S = d[s].values
    st = np.column_stack([(strata == k).astype(float) for k in [1, 2, 3]])
    cnt = 0
    for _ in range(B):
        Tp = T.copy()
        for ii in idx: Tp[ii] = rng.permutation(T[ii])
        X = np.column_stack([np.ones(len(d)), Tp, S, Tp * S, st]); X = X[:, _indep(X)]
        b = np.linalg.lstsq(X, Y, rcond=None)[0]
        cnt += abs(b[3]) >= abs(obs)
    return (cnt + 1) / (B + 1)


def holm(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0
    for r, i in enumerate(o):
        run = max(run, (m - r) * p[i]); adj[i] = min(1, run)
    return adj


def analyze(d, y, covs, label, B):
    print(f"\n######## {label}: N={len(d)}, treated={int(d.treat_now.sum())}")
    print("overall stratified-ATE (AIPW):", end=" ")
    g = aipw_scores(d, y, covs, reps=20)
    print(f"{g.mean():.3f} (se {g.std(ddof=1)/np.sqrt(len(g)):.3f})")
    rows = []
    for s in SUBS:
        bI, seI, b0, b1 = ols_int(d, s, y)
        pO = 2 * stats.norm.sf(abs(bI / seI))
        res, dA, seA, pA = sub_from_scores(g, d[s].values)
        # unadjusted stratified Hajek within subgroup
        e = stratum_prop(d); T = d.treat_now.values; Y = d[y].values
        haj = {}
        for v in [0, 1]:
            m = d[s].values == v
            w1 = T[m] / e[m]; w0 = (1 - T[m]) / (1 - e[m])
            haj[v] = (w1 * Y[m]).sum() / w1.sum() - (w0 * Y[m]).sum() / w0.sum()
        pP = perm_test(d, s, y, B=B)
        rows.append(dict(var=s, n0=res[0][2], n1=res[1][2],
                         aipw_S0=res[0][0], aipw_S0_lo=res[0][0] - 1.96 * res[0][1], aipw_S0_hi=res[0][0] + 1.96 * res[0][1],
                         aipw_S1=res[1][0], aipw_S1_lo=res[1][0] - 1.96 * res[1][1], aipw_S1_hi=res[1][0] + 1.96 * res[1][1],
                         hajek_S0=haj[0], hajek_S1=haj[1],
                         ols_int=bI, p_ols=pO, aipw_int=dA, p_aipw=pA, p_perm=pP))
    R = pd.DataFrame(rows)
    for c in ["p_ols", "p_aipw", "p_perm"]: R["holm_" + c] = holm(R[c])
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print(R.round(4).to_string(index=False))
    return R


# ---- Claim 3: weight change (kg as stored in weight_val, i.e. lb*0.45)
W = D[D.weight_val_0mo.notna() & D.weight_val_6mo.notna()].copy()
W["dw"] = W.weight_val_6mo - W.weight_val_0mo
W["dw_true_kg"] = (W.weight_lb_val_6mo - W.weight_lb_val_0mo) * 0.45359237
for c in COV: W[c] = W[c].fillna(W[c].median())
covW = COV + ["weight_val_0mo"]
print("ED split in weight sample:", W.ed_i_0mo.value_counts().to_dict(),
      "treated by ED:", W.groupby("ed_i_0mo").treat_now.sum().to_dict())
RW = analyze(W, "dw", covW, "Weight change (kg, data factor 0.45)", B=10000)
RW.to_csv("verification/independent_agent/subgroups_weight.csv", index=False)
# sensitivity: OLS interaction with baseline covariates added
X = np.column_stack([np.ones(len(W)), W.treat_now, W.ed_i_0mo, W.treat_now * W.ed_i_0mo,
                     W[[c for c in covW if c != "ed_i_0mo"]].values])
b, V = hc1(X, W.dw.values)
print("OLS T x ED with all baseline covariates: int=%.3f se=%.3f p=%.4f" % (b[3], np.sqrt(V[3, 3]), 2 * stats.norm.sf(abs(b[3] / np.sqrt(V[3, 3])))))
# robustness to outliers: drop |dw| > 20 kg
Wt = W[W.dw.abs() <= 20]
bI, seI, _, _ = ols_int(Wt, "ed_i_0mo", "dw")
print(f"Drop |dw|>20kg (n={len(Wt)}): OLS int {bI:.3f} p={2*stats.norm.sf(abs(bI/seI)):.4f}")
print("dw distribution:", W.dw.describe().round(2).to_dict())

# ---- HbA1c subgroups (claim 4)
H = D[D.hba1c_val_6mo.notna()].copy()
H["da1c"] = H.hba1c_val_6mo - H.hba1c_val_0mo
for c in COV: H[c] = H[c].fillna(H[c].median())
RH = analyze(H, "hba1c_val_6mo", COV, "HbA1c at 6 mo (level; baseline in covariates)", B=5000)
RH.to_csv("verification/independent_agent/subgroups_hba1c.csv", index=False)
RH2 = analyze(H, "da1c", COV, "HbA1c change", B=5000)
