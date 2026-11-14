"""E05: is the HbA1c (and weight) omnibus null an artefact of the causal forest's settings?
Same folds/AIPW scores as the pipeline (code/hte_lib.aipw_scores, 5-fold stratified on stratum x arm),
alternative out-of-fold CATE learners, BLP + RATE as in §7. 3 repeats.
  lin_T   : OLS T-learner (separate linear models per arm, all X)
  ridge_T : ridge T-learner (alpha=10 on standardized X)
  rf_T    : random-forest T-learner (min leaf 5)
  cf_loose: CausalForestDML, honest=False, min leaf 5, max_samples 0.5 (less regularized)
  cf_pipe : pipeline forest (fast: 500 trees)"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
from common import load_raw, analytic, HERE, ROOT  # noqa: E402
sys.path.insert(0, str(ROOT / "code"))
import hte_lib  # noqa: E402
from econml.dml import CausalForestDML  # noqa: E402

outcome = sys.argv[1] if len(sys.argv) > 1 else "hba1c"
raw = load_raw()
df, F = analytic(raw, outcome)
X = df[F].to_numpy(float); Y = df.Y.to_numpy(float); W = df.treat_now.to_numpy(int)
e = df.E_HAT.to_numpy(); strata = df.stratum.to_numpy(int)
cfg = {"n_trees": 500, "n_trees_nuis": 200, "k_folds": 5, "min_leaf_cf": 10,
       "stratum_cols": [F.index("rural_site"), F.index("hba1c_cat_control")]}


def toc(b, p):
    o = np.argsort(-p, kind="stable"); bb = b[o]
    return np.cumsum(bb) / np.arange(1, len(bb) + 1) - bb.mean()


def tests(gam, tau, seed):
    m = sm.OLS(gam, sm.add_constant(tau - tau.mean())).fit(cov_type="HC1")
    a = toc(-gam, -tau).mean()
    rb = np.random.default_rng(seed)
    bs = [toc(-gam[i], -tau[i]).mean() for i in (rb.integers(0, len(gam), len(gam)) for _ in range(500))]
    return m.params[1], 1 - norm.cdf(m.params[1] / m.bse[1]), a, 1 - norm.cdf(a / np.std(bs, ddof=1))


rows = []
for rep in range(3):
    skf = StratifiedKFold(5, shuffle=True, random_state=20260621 + rep)
    taus = {k: np.zeros(len(Y)) for k in ("lin_T", "ridge_T", "rf_T", "cf_loose", "cf_pipe")}
    gam = np.zeros(len(Y))
    for k, (tr, te) in enumerate(skf.split(X, strata * 2 + W)):
        s = 20260621 + 100 * rep + k
        gam[te] = hte_lib.aipw_scores(Y, W, X, e, tr, te, s, 200)
        sc = StandardScaler().fit(X[tr])
        Xtr, Xte = sc.transform(X[tr]), sc.transform(X[te])
        for name, mk in (("lin_T", lambda: LinearRegression()), ("ridge_T", lambda: Ridge(alpha=10.0)),
                         ("rf_T", lambda: RandomForestRegressor(300, min_samples_leaf=5, random_state=s, n_jobs=1))):
            mu = {a: mk().fit(Xtr[W[tr] == a], Y[tr][W[tr] == a]).predict(Xte) for a in (0, 1)}
            taus[name][te] = mu[1] - mu[0]
        cfl = CausalForestDML(model_y=RandomForestRegressor(200, min_samples_leaf=5, random_state=s, n_jobs=1),
                              model_t=hte_lib.StratumPropensity(cols=tuple(cfg["stratum_cols"])),
                              discrete_treatment=True, cv=5, n_estimators=500, min_samples_leaf=5,
                              max_samples=0.5, honest=False, random_state=s, n_jobs=1).fit(Y[tr], W[tr], X=X[tr])
        taus["cf_loose"][te] = cfl.effect(X[te])
        taus["cf_pipe"][te] = hte_lib.make_forest(s, cfg).fit(Y[tr], W[tr], X=X[tr]).effect(X[te])
    for name, t in taus.items():
        b2, p1, a, pa = tests(gam, t, rep)
        rows.append(dict(outcome=outcome, repeat=rep, learner=name, tau_sd=t.std(), blp_b2=b2, blp_p1=p1,
                         autoc=a, autoc_p1=pa))
        print(rows[-1], flush=True)
R = pd.DataFrame(rows)
R.to_csv(HERE / f"e05_alt_learners_{outcome}.csv", index=False)
print(R.groupby("learner")[["tau_sd", "blp_b2", "blp_p1", "autoc", "autoc_p1"]].median().round(3).to_string())
