"""
Agent A: independent recomputation of the two manuscript P values for weight x prior-year ED
that no script in code/ writes ("linear-model and permutation P = .002 and .003";
docs/04 table: OLS HC1 0.0025, within-stratum permutation 0.0028, N=392).

    python check_weight_ed_ols_perm.py <run_root>

Uses <run_root>/outputs/v3/weight/patient_level_cate.csv (Y, treat_now, subject_id = the
N=392 analytic sample) merged with the raw data for ed_i_0mo and the randomization strata.
OLS: Y ~ T + ED + T:ED + stratum dummies, HC1. Permutation: T permuted within stratum
(10 000 times, seed fixed), statistic = |OLS interaction t|.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

root = Path(sys.argv[1])
pl = pd.read_csv(root / "outputs/v3/weight/patient_level_cate.csv")
raw = pd.read_csv(root / "data/raw/dataverse_files/analysis_sample.csv")
d = pl.merge(raw[["subject_id", "ed_i_0mo", "rural_site", "hba1c_cat_control"]], on="subject_id",
             how="left", validate="one_to_one")
d["stratum"] = d.rural_site * 2 + d.hba1c_cat_control
d["T"] = d.treat_now
assert d.ed_i_0mo.notna().all()
f = "Y ~ T * ed_i_0mo + C(stratum)"
m = smf.ols(f, data=d).fit(cov_type="HC1")
print(f"N={len(d)}, ED=0 n={int((d.ed_i_0mo == 0).sum())}, ED=1 n={int((d.ed_i_0mo == 1).sum())}")
print(f"OLS HC1 interaction {m.params['T:ed_i_0mo']:.4f} (SE {m.bse['T:ed_i_0mo']:.4f}), "
      f"p = {m.pvalues['T:ed_i_0mo']:.4f}")

t_obs = abs(m.tvalues["T:ed_i_0mo"])
rng = np.random.default_rng(20261029)
groups = [np.where(d.stratum.values == s)[0] for s in np.unique(d.stratum)]
B, hits = 10000, 0
X = d.copy()
for b in range(B):
    T = d["T"].values.copy()
    for g in groups:
        T[g] = rng.permutation(T[g])
    X["T"] = T
    mb = smf.ols(f, data=X).fit(cov_type="HC1")
    hits += abs(mb.tvalues["T:ed_i_0mo"]) >= t_obs
p_perm = (hits + 1) / (B + 1)
print(f"Within-stratum permutation (B={B}, |t| statistic): p = {p_perm:.4f}")
