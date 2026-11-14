"""Agent B independent rebuild of the HTE analysis (spec: docs/02_analysis_plan.md, sap-v3).

Usage: python b_analysis.py <outcome> [--pc]      (outcome in hba1c, weight, sbp, ed_ip)
Writes out/<outcome>.json (and out/<outcome>_pc.json with --pc).
"""
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import StratifiedKFold

from b_data import load_raw, analytic, SUBGROUPS

warnings.filterwarnings("ignore")
OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)

SEED = 20261029
N_REPEATS = 5
N_FOLDS = 5
N_JOBS = 3
N_TREES_CF = 4000
N_BOOT_RATE = 1000
N_RERAND = 2000
N_PERM_SUB = 2000


# ----------------------------------------------------------------------------- helpers
def folds_for(strat_label, rep):
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED + 1000 * rep)
    fold = np.empty(len(strat_label), dtype=int)
    for k, (_, te) in enumerate(skf.split(np.zeros(len(strat_label)), strat_label)):
        fold[te] = k
    return fold


def rf(seed, leaf=5):
    return RandomForestRegressor(n_estimators=500, min_samples_leaf=leaf, max_features=1 / 3,
                                 random_state=seed, n_jobs=N_JOBS)


def ols_hc1(y, Z):
    """OLS with HC1 covariance. Returns beta, cov."""
    n, k = Z.shape
    ZtZi = np.linalg.pinv(Z.T @ Z)
    b = ZtZi @ Z.T @ y
    r = y - Z @ b
    meat = (Z * r[:, None] ** 2).T @ Z
    cov = ZtZi @ meat @ ZtZi * n / (n - k)
    return b, cov


def wald(b, cov, idx):
    bb = b[idx]
    V = cov[np.ix_(idx, idx)]
    stat = float(bb @ np.linalg.pinv(V) @ bb)
    return stat, float(stats.chi2.sf(stat, len(idx)))


def holm(p):
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m)
    run = 0.0
    for i, j in enumerate(order):
        run = max(run, (m - i) * p[j])
        adj[j] = min(1.0, run)
    return adj


def within_stratum_permute(W, strata, rng):
    Wp = W.copy()
    for s in np.unique(strata):
        idx = np.where(strata == s)[0]
        Wp[idx] = W[rng.permutation(idx)]
    return Wp


class KnownPropensity(BaseEstimator, ClassifierMixin):
    """Returns the known stratum propensity; strata identified from two binary columns of X."""

    def __init__(self, i_rural=0, i_cat=1, table=None):
        self.i_rural = i_rural
        self.i_cat = i_cat
        self.table = table

    def fit(self, X, y, **kw):
        self.classes_ = np.array([0, 1])
        return self

    def predict_proba(self, X):
        X = np.asarray(X)
        s = (2 * np.round(X[:, self.i_rural]) + np.round(X[:, self.i_cat])).astype(int)
        e = np.array([self.table[int(k)] for k in s])
        return np.column_stack([1 - e, e])

    def predict(self, X):
        return (self.predict_proba(X)[:, 1] > 0.5).astype(int)


# ----------------------------------------------------------------------------- AIPW
def aipw_scores(X, Y, W, e, fold, seed):
    mu0 = np.empty(len(Y))
    mu1 = np.empty(len(Y))
    for k in np.unique(fold):
        tr, te = fold != k, fold == k
        m1 = rf(seed + k).fit(X[tr & (W == 1)], Y[tr & (W == 1)])
        m0 = rf(seed + 50 + k).fit(X[tr & (W == 0)], Y[tr & (W == 0)])
        mu1[te] = m1.predict(X[te])
        mu0[te] = m0.predict(X[te])
    G = mu1 - mu0 + W * (Y - mu1) / e - (1 - W) * (Y - mu0) / (1 - e)
    return G


def mean_ci(g):
    n = len(g)
    m = float(np.mean(g))
    se = float(np.std(g, ddof=1) / np.sqrt(n))
    return m, se, m - 1.96 * se, m + 1.96 * se, float(2 * stats.norm.sf(abs(m / se)))


# ----------------------------------------------------------------------------- causal forest OOF
def cf_oof(X, Y, W, fold, rep, i_rural, i_cat, e_table):
    from econml.dml import CausalForestDML
    tau = np.empty(len(Y))
    for k in np.unique(fold):
        tr, te = fold != k, fold == k
        cf = CausalForestDML(
            model_y=rf(SEED + 7 * rep + k),
            model_t=KnownPropensity(i_rural, i_cat, e_table),
            discrete_treatment=True, n_estimators=N_TREES_CF, min_samples_leaf=10,
            max_samples=0.45, honest=True, cv=5, inference=False,
            random_state=SEED + 100 * rep + k, n_jobs=N_JOBS)
        cf.fit(Y[tr], W[tr], X=X[tr])
        tau[te] = cf.effect(X[te]).ravel()
    return tau


def blp_test(G, tau):
    tc = tau - tau.mean()
    Z = np.column_stack([np.ones(len(G)), tc])
    b, cov = ols_hc1(G, Z)
    t = b[1] / np.sqrt(cov[1, 1])
    return dict(beta1=float(b[0]), beta2=float(b[1]), se2=float(np.sqrt(cov[1, 1])),
                p_one_sided=float(stats.norm.sf(t)))


def autoc_stat(G_sorted):
    n = len(G_sorted)
    toc = np.cumsum(G_sorted) / np.arange(1, n + 1) - G_sorted.mean()
    return toc.mean()


def rate_autoc(G, tau, rng):
    order = np.argsort(-tau, kind="stable")
    Gs = G[order]
    est = autoc_stat(Gs)
    n = len(G)
    boots = np.empty(N_BOOT_RATE)
    for b in range(N_BOOT_RATE):
        idx = np.sort(rng.integers(0, n, n))  # sorted positions keep the tau ranking
        boots[b] = autoc_stat(Gs[idx])
    se = boots.std(ddof=1)
    return dict(autoc=float(est), se=float(se), p_one_sided=float(stats.norm.sf(est / se)))


def gates(G, tau):
    q = pd.qcut(pd.Series(tau).rank(method="first"), 3, labels=False).values
    return [float(G[q == g].mean()) for g in range(3)]


# ----------------------------------------------------------------------------- parametric interaction test
def interaction_design(Xs, W):
    Xc = Xs - Xs.mean(0)
    return np.column_stack([np.ones(len(W)), W, Xc, W[:, None] * Xc])


def interaction_wald(Y, Xs, W):
    Z = interaction_design(Xs, W)
    b, cov = ols_hc1(Y, Z)
    p = Xs.shape[1]
    idx = list(range(2 + p, 2 + 2 * p))
    return wald(b, cov, idx)


def interaction_test(Y, Xs, W, strata, rng):
    stat, p_asym = interaction_wald(Y, Xs, W)
    # constant-effect null (Ding, Feller & Miratrix 2016): plug-in tau from main-effects ANCOVA
    Z0 = np.column_stack([np.ones(len(W)), W, Xs - Xs.mean(0)])
    b0, _ = ols_hc1(Y, Z0)
    tau0 = b0[1]
    Y0 = Y - tau0 * W
    cnt = 0
    for _ in range(N_RERAND):
        Wp = within_stratum_permute(W, strata, rng)
        s, _ = interaction_wald(Y0 + tau0 * Wp, Xs, Wp)
        cnt += s >= stat
    return dict(wald=stat, df=Xs.shape[1], p_asymptotic=p_asym, tau0=float(tau0),
                p_randomization=float((1 + cnt) / (1 + N_RERAND)), n_rerand=N_RERAND)


# ----------------------------------------------------------------------------- main
def run(outcome, pc=False):
    t0 = time.time()
    d, quirks, e_full = load_raw()
    df, X_cols = analytic(d, outcome)
    X = df[X_cols].to_numpy(float)
    Y = df["Y"].to_numpy(float)
    W = df["treat_now"].to_numpy(int)
    strata = df["stratum"].to_numpy(int)
    e_an = df.groupby("stratum")["treat_now"].mean().to_dict()
    e = df["stratum"].map(e_an).to_numpy(float)          # primary: analytic-sample stratum share
    e_alt = df["stratum"].map(e_full).to_numpy(float)    # sensitivity: full-465 stratum share
    strat_label = strata * 2 + W
    i_rural, i_cat = X_cols.index("rural_site"), X_cols.index("hba1c_cat_control")

    res = dict(outcome=outcome, n=len(Y), n_treat=int(W.sum()), modifiers=X_cols,
               e_analytic=e_an, e_full465=e_full, y_sd=float(Y.std(ddof=1)), quirks=quirks)

    pcs = [None]
    if pc:
        older = (df["age"] >= df["age"].median()).to_numpy(int)
        res["pc_older_def"] = f"age >= median ({df['age'].median()})"
        res["pc_older_n"] = int(older.sum())
        pcs = [0.25, 0.5, 1.0]
        res["pc"] = {}

    for delta in pcs:
        Yuse = Y if delta is None else Y + delta * res["y_sd"] * W * older
        G_reps, tau_reps, blps, rates, gts, ate_reps = [], [], [], [], [], []
        for rep in range(N_REPEATS):
            fold = folds_for(strat_label, rep)
            G = aipw_scores(X, Yuse, W, e, fold, SEED + 10 * rep)
            tau = cf_oof(X, Yuse, W, fold, rep, i_rural, i_cat, e_an)
            G_reps.append(G)
            tau_reps.append(tau)
            ate_reps.append(mean_ci(G))
            blps.append(blp_test(G, tau))
            rates.append(rate_autoc(G, tau, np.random.default_rng(SEED + rep)))
            gts.append(gates(G, tau))
            print(outcome, delta, "rep", rep, f"{time.time() - t0:.0f}s", blps[-1]["p_one_sided"],
                  rates[-1]["p_one_sided"], flush=True)
        G_bar = np.mean(G_reps, axis=0)
        p_blp = float(np.median([b["p_one_sided"] for b in blps]))
        p_rate = float(np.median([r["p_one_sided"] for r in rates]))
        omni = dict(
            blp_beta2_median=float(np.median([b["beta2"] for b in blps])),
            blp_beta1_median=float(np.median([b["beta1"] for b in blps])),
            blp_p_median=p_blp, blp_reps=blps,
            autoc_median=float(np.median([r["autoc"] for r in rates])),
            rate_p_median=p_rate, rate_reps=rates,
            gates_median=np.median(np.array(gts), axis=0).tolist(),
            decision_detected=bool(p_blp < 0.05 and p_rate < 0.05),
            tau_oof_sd_mean=float(np.mean([t.std() for t in tau_reps])),
        )
        if delta is not None:
            res["pc"][str(delta)] = omni
            continue
        res["omnibus"] = omni
        # ATE: repeat-averaged scores (primary) and median-aggregated (secondary)
        m, se, lo, hi, p = mean_ci(G_bar)
        res["ate_aipw"] = dict(est=m, se=se, lo=lo, hi=hi, p=p)
        ests = np.array([a[0] for a in ate_reps])
        ses = np.array([a[1] for a in ate_reps])
        med = np.median(ests)
        se_med = float(np.median(np.sqrt(ses ** 2 + (ests - med) ** 2)))
        res["ate_aipw_median_agg"] = dict(est=float(med), se=se_med, lo=float(med - 1.96 * se_med),
                                          hi=float(med + 1.96 * se_med), reps=ests.tolist())
        # ATE sensitivity: full-465 stratum propensity
        G_alt = np.mean([aipw_scores(X, Y, W, e_alt, folds_for(strat_label, r), SEED + 10 * r)
                         for r in range(N_REPEATS)], axis=0)
        m2, se2, lo2, hi2, _ = mean_ci(G_alt)
        res["ate_aipw_e_full465"] = dict(est=m2, se=se2, lo=lo2, hi=hi2)
        res["ate_unadjusted_diff_means"] = float(Y[W == 1].mean() - Y[W == 0].mean())
        # BLP of AIPW scores on covariates (joint Wald)
        Xs = (X - X.mean(0)) / X.std(0)
        b, cov = ols_hc1(G_bar, np.column_stack([np.ones(len(Y)), Xs]))
        res["blp_covariates_joint"] = dict(zip(["wald", "p"], wald(b, cov, list(range(1, Xs.shape[1] + 1)))))
        # parametric interaction test
        res["interaction_test"] = interaction_test(Y, X, W, strata, np.random.default_rng(SEED + 77))
        # subgroups
        rows = []
        for name, col in SUBGROUPS.items():
            g = df[col].to_numpy(int)
            lv = {}
            for level in (0, 1):
                mk = g == level
                mm, ss, ll, hh, _ = mean_ci(G_bar[mk])
                lv[level] = (mm, ss)
                rows.append(dict(outcome=outcome, subgroup=name, variable=col, level=level,
                                 n=int(mk.sum()), n_treat=int(W[mk].sum()),
                                 est=mm, se=ss, lo=ll, hi=hh))
            z = (lv[1][0] - lv[0][0]) / np.sqrt(lv[1][1] ** 2 + lv[0][1] ** 2)
            pint = float(2 * stats.norm.sf(abs(z)))
            for r in rows[-2:]:
                r["diff_1_minus_0"] = lv[1][0] - lv[0][0]
                r["p_interaction"] = pint
        ps = [rows[2 * i]["p_interaction"] for i in range(len(SUBGROUPS))]
        ph = holm(ps)
        for i in range(len(SUBGROUPS)):
            rows[2 * i]["p_holm"] = rows[2 * i + 1]["p_holm"] = float(ph[i])
        res["subgroups"] = rows
        if outcome == "weight":
            res["weight_x_prior_ed"] = weight_ed_checks(df, Y, W, strata)
    res["runtime_s"] = time.time() - t0
    fn = OUT / (f"{outcome}_pc.json" if pc else f"{outcome}.json")
    fn.write_text(json.dumps(res, indent=1, default=float))
    print("done", fn, f"{res['runtime_s']:.0f}s")


def weight_ed_checks(df, Y, W, strata):
    g = df["ed_i_0mo"].to_numpy(float)
    S = pd.get_dummies(df["stratum"], prefix="s", drop_first=True).to_numpy(float)

    def fit(Yv, Wv):
        Z = np.column_stack([np.ones(len(Yv)), Wv, g, Wv * g, S])
        b, cov = ols_hc1(Yv, Z)
        return b, cov

    b, cov = fit(Y, W)
    t_obs = b[3] / np.sqrt(cov[3, 3])
    out = dict(coef_interaction=float(b[3]), se_hc1=float(np.sqrt(cov[3, 3])),
               lo=float(b[3] - 1.96 * np.sqrt(cov[3, 3])), hi=float(b[3] + 1.96 * np.sqrt(cov[3, 3])),
               p_hc1=float(2 * stats.t.sf(abs(t_obs), len(Y) - 4 - S.shape[1])),
               effect_no_ed=float(b[1]), effect_ed=float(b[1] + b[3]))
    rng = np.random.default_rng(SEED + 99)
    # (a) sharp-null permutation of W within strata (no effect for anyone), stat = |t|
    # (b) constant-effect null: Y0 = Y - tau0*W with tau0 from the no-interaction model
    Z0 = np.column_stack([np.ones(len(Y)), W, g, S])
    b0, _ = ols_hc1(Y, Z0)
    Y0 = Y - b0[1] * W
    ca = cb = 0
    for _ in range(N_PERM_SUB):
        Wp = within_stratum_permute(W, strata, rng)
        bp, cp = fit(Y, Wp)
        ca += abs(bp[3] / np.sqrt(cp[3, 3])) >= abs(t_obs)
        bq, cq = fit(Y0 + b0[1] * Wp, Wp)
        cb += abs(bq[3] / np.sqrt(cq[3, 3])) >= abs(t_obs)
    out["p_perm_sharp_null"] = float((1 + ca) / (1 + N_PERM_SUB))
    out["p_perm_constant_effect"] = float((1 + cb) / (1 + N_PERM_SUB))
    out["n_perm"] = N_PERM_SUB
    return out


if __name__ == "__main__":
    run(sys.argv[1], pc="--pc" in sys.argv)
