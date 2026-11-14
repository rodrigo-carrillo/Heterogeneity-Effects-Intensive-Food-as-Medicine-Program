"""Alternative estimators and fragility checks for the weight x prior-year ED interaction
(post hoc: added after the signal was seen; not in docs/02).

    python code/weight_ed_checks.py   ->  outputs/v3/weight/weight_ed_checks.{json,txt}

Uses the N = 392 analytic sample of the main weight run (outputs/v3/weight/patient_level_cate.csv,
so run code/run_all.py first) merged with prior-year ED and the randomization strata.
  OLS          : Y ~ T * ED + stratum dummies, HC1 SE.
  Permutation  : T permuted within stratum, statistic |t| of the OLS interaction, 10 000 draws,
                 under the sharp null (no effect for anyone) and under a constant-effect null
                 (Y - tau*T, tau from the OLS without the interaction), which is the null that
                 matters for an interaction.
  Fragility    : OLS after excluding |weight change| > 10 kg.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "outputs" / "v3" / "weight"
SEED, B = 20261102, 10_000
F = "Y ~ T * ed_i_0mo + C(stratum)"

pl = pd.read_csv(RUN / "patient_level_cate.csv")
raw = pd.read_csv(ROOT / "data" / "raw" / "dataverse_files" / "analysis_sample.csv")
d = pl.merge(raw[["subject_id", "ed_i_0mo", "rural_site", "hba1c_cat_control"]],
             on="subject_id", how="left", validate="one_to_one")
d["stratum"] = 2 * d.rural_site + d.hba1c_cat_control
d["T"] = d.treat_now.astype(int)
assert len(d) == 392 and d.ed_i_0mo.notna().all()


def ols(data):
    m = smf.ols(F, data=data).fit(cov_type="HC1")
    k = "T:ed_i_0mo"
    return m.params[k], m.bse[k], m.tvalues[k], m.pvalues[k]


b, se, t_obs, p = ols(d)
tau = smf.ols("Y ~ T + ed_i_0mo + C(stratum)", data=d).fit().params["T"]

rng = np.random.default_rng(SEED)
groups = [np.flatnonzero(d.stratum.to_numpy() == s) for s in np.unique(d.stratum)]
T0, Y0 = d["T"].to_numpy(), d["Y"].to_numpy()
hits = {"sharp": 0, "constant": 0}
X = d.copy()
for _ in range(B):
    T = T0.copy()
    for g in groups:
        T[g] = rng.permutation(T[g])
    X["T"] = T
    X["Y"] = Y0
    hits["sharp"] += abs(ols(X)[2]) >= abs(t_obs)
    X["Y"] = Y0 - tau * T0 + tau * T
    hits["constant"] += abs(ols(X)[2]) >= abs(t_obs)

trim = d[(d.Y.abs() <= 10)]
bt, set_, _, pt = ols(trim)

res = {
    "n": len(d), "n_no_ed": int((d.ed_i_0mo == 0).sum()), "n_ed": int((d.ed_i_0mo == 1).sum()),
    "ols_interaction": b, "ols_se": se, "ols_p_hc1": p,
    "perm_B": B, "perm_seed": SEED, "constant_effect_tau": tau,
    "perm_p_sharp_null": (hits["sharp"] + 1) / (B + 1),
    "perm_p_constant_effect_null": (hits["constant"] + 1) / (B + 1),
    "trim10_n": len(trim), "trim10_excluded": len(d) - len(trim),
    "trim10_interaction": bt, "trim10_se": set_, "trim10_p_hc1": pt,
}
(RUN / "weight_ed_checks.json").write_text(json.dumps(res, indent=2))
lines = [f"{k}: {v:.6g}" if isinstance(v, float) else f"{k}: {v}" for k, v in res.items()]
(RUN / "weight_ed_checks.txt").write_text("\n".join(lines) + "\n")
sys.stdout.write("\n".join(lines) + "\n")
