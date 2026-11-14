"""
Model-fitting helpers for code/hte_food_as_medicine.py, plus a worker entry point.

Why a separate worker process: EconML's compiled GRF tree builder segfaulted
intermittently (2026-10-21; EconML 0.16 and 0.17, Python 3.12 and 3.13, macOS),
always after ~100 forest fits in one long-lived process and never when the same
fits were repeated in a fresh process. Each cross-fitting repeat therefore runs
in its own short-lived process (``python code/hte_lib.py in.npz out.npz cfg``),
and a crashed repeat is simply rerun. Seeds are fixed, so reruns are identical.
"""
import json
import sys

import numpy as np
from econml.dml import CausalForestDML
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import StratifiedKFold


class StratumPropensity(BaseEstimator, ClassifierMixin):
    """P(W=1 | stratum) estimated as the treated fraction in each stratum.

    `cols` are the column indices (in the feature matrix) of the variables that
    define the randomization strata (here rural_site and hba1c_cat_control)."""
    def __init__(self, cols=(0, 1)):
        self.cols = cols

    def _key(self, X):
        X = np.asarray(X)
        return [tuple(r) for r in np.round(X[:, list(self.cols)]).astype(int)]

    def fit(self, X, y, sample_weight=None):
        y = np.asarray(y).ravel()
        self.classes_ = np.unique(y)
        keys = self._key(X)
        self.p_ = {k: np.mean(y[[kk == k for kk in keys]] == self.classes_[1]) for k in set(keys)}
        self.p_all_ = np.mean(y == self.classes_[1])
        return self

    def predict_proba(self, X):
        p = np.array([self.p_.get(k, self.p_all_) for k in self._key(X)])
        return np.column_stack([1 - p, p])

    def predict(self, X):
        return self.classes_[(self.predict_proba(X)[:, 1] > 0.5).astype(int)]


def make_forest(seed, cfg):
    """Causal forest (generalized random forest; Athey, Tibshirani & Wager 2019)
    with honest splitting and cross-fitted nuisance models. n_jobs=1: the
    multithreaded builder crashed more often; parallelism is across processes."""
    return CausalForestDML(
        model_y=RandomForestRegressor(n_estimators=cfg["n_trees_nuis"], min_samples_leaf=5,
                                      random_state=seed, n_jobs=1),
        model_t=StratumPropensity(cols=tuple(cfg["stratum_cols"])),
        discrete_treatment=True, cv=cfg["k_folds"], n_estimators=cfg["n_trees"],
        min_samples_leaf=cfg["min_leaf_cf"], max_samples=0.45, honest=True,
        inference=True, random_state=seed, n_jobs=1)


def aipw_scores(Y, W, X, e, train, test, seed, n_trees_nuis):
    """Out-of-fold AIPW scores for the test indices:
    Gamma = mu1 - mu0 + W(Y - mu1)/e - (1 - W)(Y - mu0)/(1 - e),
    with arm-specific random-forest outcome models fit on the train indices."""
    mu = {}
    for arm in (0, 1):
        idx = train[W[train] == arm]
        m = RandomForestRegressor(n_estimators=n_trees_nuis, min_samples_leaf=5,
                                  random_state=seed, n_jobs=1).fit(X[idx], Y[idx])
        mu[arm] = m.predict(X[test])
    Wt, Yt, et = W[test], Y[test], e[test]
    return mu[1] - mu[0] + Wt * (Yt - mu[1]) / et - (1 - Wt) * (Yt - mu[0]) / (1 - et)


def crossfit_repeat(Y, W, X, e, strata, repeat, seed0, cfg):
    """One K-fold partition (stratified on stratum x arm): out-of-fold CATE and AIPW scores."""
    n = len(Y)
    tau, gamma = np.zeros(n), np.zeros(n)
    skf = StratifiedKFold(n_splits=cfg["k_folds"], shuffle=True, random_state=seed0 + repeat)
    for k, (tr, te) in enumerate(skf.split(X, strata * 2 + W)):
        s = seed0 + 100 * repeat + k
        forest = make_forest(s, cfg).fit(Y[tr], W[tr], X=X[tr])
        tau[te] = forest.effect(X[te])
        gamma[te] = aipw_scores(Y, W, X, e, tr, te, s, cfg["n_trees_nuis"])
    return tau, gamma


if __name__ == "__main__":
    # Worker: python hte_lib.py <in.npz> <out.npz> '<json cfg incl. repeat, seed0>'
    inp, out, cfg = sys.argv[1], sys.argv[2], json.loads(sys.argv[3])
    d = np.load(inp)
    tau, gamma = crossfit_repeat(d["Y"], d["W"].astype(int), d["X"], d["e"],
                                 d["strata"].astype(int), cfg["repeat"], cfg["seed0"], cfg)
    np.savez(out, tau=tau, gamma=gamma)
