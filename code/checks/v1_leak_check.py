"""
Reproduce the v1 leakage diagnosis (docs/03_v1_review.md, issue #1).

v1 (code/archive/v1/hte_food_as_medicine.ipynb, the spec that produced
outputs/v1_2026-06/) scored the BLP calibration test on in-sample forest
predictions, tau_hat = cf.effect(X), for the same patients whose outcomes feed
the DR scores. This script refits v1's exact specification (covariates, median
imputation + missingness flags, forest settings, RF propensity, DR scores, seed)
and runs the BLP twice:
  1. in-sample tau_hat, as v1 did;
  2. out-of-fold tau_hat: 5-fold, each patient's CATE from a forest that never
     saw that patient.
Documented result (2026-10-21, macOS): in-sample slope +11.2 (p<0.001) vs
out-of-fold slope -2.8. Forest numbers can move slightly across platforms.

  python code/checks/v1_leak_check.py
"""
# Rewritten 2026-10-22 from the v1 notebook and docs/03 after the original was
# lost when the project folder was uploaded to OneDrive.
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from econml.dml import CausalForestDML
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.model_selection import KFold

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw" / "dataverse_files" / "analysis_sample.csv"
RNG = 20260621  # v1 seed

# ---- v1 data preparation (notebook cells 1-4; lipids commented out) ----------
df = pd.read_csv(DATA)
df = df[df["hba1c_val_6mo"].notna()].copy()
Y = (df["hba1c_val_6mo"] - df["hba1c_val_0mo"]).values
W = df["treat_now"].astype(int).values
continuous = ["age", "hba1c_val_0mo", "bmi_val_0mo", "sys_val_0mo", "dias_val_0mo",
              "op_n_0mo", "ed_or_ip_n_0mo"]
binary = ["male", "hispanic", "race_nonhisp_black", "rural_site", "drug_insulin_0mo",
          "drug_metformin_0mo", "drug_glp1_0mo", "drug_sulf_0mo", "drug_other_0mo",
          "drug_any_db_0mo", "ip_i_0mo", "ed_i_0mo"]
X = df[continuous + binary]
miss = X.isna().astype(int).add_suffix("_miss")
miss = miss[[c for c in miss.columns if miss[c].mean() >= 0.05]]
X_imp = pd.DataFrame(SimpleImputer(strategy="median").fit_transform(X),
                     columns=X.columns, index=X.index)
Xv = pd.concat([X_imp, miss], axis=1).values
print(f"N = {len(Y)} (treated {W.sum()}), {Xv.shape[1]} features")


def v1_forest():
    """v1 notebook cell 5."""
    return CausalForestDML(
        n_estimators=4000, min_samples_leaf=10, max_depth=None, max_samples=0.45,
        discrete_treatment=True, cv=5,
        model_y=RandomForestRegressor(n_estimators=500, min_samples_leaf=5, random_state=RNG),
        model_t=RandomForestClassifier(n_estimators=500, min_samples_leaf=5, random_state=RNG),
        random_state=RNG)


def dr_scores(Y, W, X, n_splits=5, seed=RNG):
    """v1 notebook cell 6: cross-fitted AIPW scores, RF propensity clipped to [0.05, 0.95]."""
    dr = np.zeros_like(Y, dtype=float)
    for tr, te in KFold(n_splits=n_splits, shuffle=True, random_state=seed).split(X):
        rf = dict(n_estimators=500, min_samples_leaf=5, random_state=seed)
        mu1 = RandomForestRegressor(**rf).fit(X[tr][W[tr] == 1], Y[tr][W[tr] == 1]).predict(X[te])
        mu0 = RandomForestRegressor(**rf).fit(X[tr][W[tr] == 0], Y[tr][W[tr] == 0]).predict(X[te])
        e = np.clip(RandomForestClassifier(**rf).fit(X[tr], W[tr]).predict_proba(X[te])[:, 1],
                    0.05, 0.95)
        dr[te] = mu1 - mu0 + W[te] * (Y[te] - mu1) / e - (1 - W[te]) * (Y[te] - mu0) / (1 - e)
    return dr


def blp(dr, tau, label):
    fit = sm.OLS(dr, sm.add_constant(tau - tau.mean())).fit(cov_type="HC1")
    print(f"{label:34s} beta1 {fit.params[0]:+.3f}   beta2 {fit.params[1]:+.2f} "
          f"(SE {fit.bse[1]:.2f}, two-sided p {fit.pvalues[1]:.4f})")
    return fit.params[1]


dr = dr_scores(Y, W, Xv)
print(f"Mean DR score (ATE): {dr.mean():+.3f}")

# 1. v1: in-sample predictions
tau_in = v1_forest().fit(Y=Y, T=W, X=Xv).effect(Xv)

# 2. out-of-fold predictions, same forest settings
tau_oof = np.zeros_like(Y, dtype=float)
for k, (tr, te) in enumerate(KFold(n_splits=5, shuffle=True, random_state=RNG).split(Xv)):
    tau_oof[te] = v1_forest().fit(Y=Y[tr], T=W[tr], X=Xv[tr]).effect(Xv[te])
    print(f"  out-of-fold fold {k + 1}/5 done")

print("\nBLP calibration test (DR ~ 1 + centred tau_hat, HC1)")
b_in = blp(dr, tau_in, "in-sample tau_hat (v1):")
b_oof = blp(dr, tau_oof, "out-of-fold tau_hat (v2+):")
print(f"\nDocumented 2026-10-21: in-sample +11.2 (p<0.001), out-of-fold -2.8.")
print(f"This run:              in-sample {b_in:+.1f}, out-of-fold {b_oof:+.1f}.")
