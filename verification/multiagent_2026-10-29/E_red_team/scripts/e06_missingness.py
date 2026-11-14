"""E06: differential inclusion (complete-case) by arm, overall and by arm x subgroup, per outcome.
Also: does arm-differential inclusion vary with baseline covariates (joint arm x X test)?"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import chi2

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
from common import load_raw, analytic, PRESPEC, HERE  # noqa: E402

raw = load_raw()
always = ["age", "hba1c_val_0mo", "days_a1c_to_consent", "op_n_0mo", "ed_or_ip_n_0mo", "male",
          "hispanic", "race_nonhisp_black", "rural_site", "hba1c_cat_control", "drug_insulin_0mo",
          "drug_metformin_0mo", "drug_glp1_0mo", "drug_sulf_0mo", "drug_other_0mo", "ip_i_0mo", "ed_i_0mo"]
rows = []
for oc in ("hba1c", "weight", "sbp", "ed_ip"):
    for ldl in (False, True):
        df, F = analytic(raw, oc, ldl)
        r = raw.copy()
        r["inc"] = r.subject_id.isin(df.subject_id).astype(int)
        inc_by_arm = r.groupby("treat_now").inc.mean()
        m0 = smf.logit("inc ~ treat_now + " + " + ".join(always), data=r).fit(disp=0, method="bfgs", maxiter=2000)
        try:
            m1 = smf.logit("inc ~ treat_now * (" + " + ".join(always) + ")", data=r).fit(disp=0, method="bfgs", maxiter=2000)
            lr = 2 * (m1.llf - m0.llf); dfree = m1.df_model - m0.df_model
        except Exception:
            lr, dfree = np.nan, 1
        row = dict(outcome=oc, ldl=ldl, n_inc=int(r.inc.sum()), inc_ctrl=inc_by_arm[0], inc_trt=inc_by_arm[1],
                   p_arm=smf.logit("inc ~ treat_now", data=r).fit(disp=0).pvalues["treat_now"],
                   p_arm_x_X_joint=chi2.sf(lr, dfree))
        for v in PRESPEC:
            try:
                mm = smf.logit(f"inc ~ treat_now * {v}", data=r).fit(disp=0, method="bfgs", maxiter=2000)
                row[f"p_arm_x_{v}"] = mm.pvalues[f"treat_now:{v}"]
            except Exception:
                row[f"p_arm_x_{v}"] = np.nan
        rows.append(row)
R = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(R.round(3).T.to_string())
R.to_csv(HERE / "e06_missingness.csv", index=False)
# Detail for weight and HbA1c x prior ED
for oc in ("weight", "hba1c"):
    df, F = analytic(raw, oc)
    r = raw.copy(); r["inc"] = r.subject_id.isin(df.subject_id).astype(int)
    print(f"\n{oc}: inclusion rate by prior ED x arm")
    print(r.groupby(["ed_i_0mo", "treat_now"]).inc.agg(["size", "mean"]).round(3).to_string())
