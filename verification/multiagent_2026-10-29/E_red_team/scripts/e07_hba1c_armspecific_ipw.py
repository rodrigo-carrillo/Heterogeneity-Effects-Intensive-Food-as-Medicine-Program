"""E07: HbA1c complete-case inclusion differs by arm x prior ED (E06). The pipeline's S3 inclusion
model is additive in arm and X, so it cannot correct arm-specific selection. Refit the 7 Doyle
subgroup interactions and the overall effect with ARM-SPECIFIC inclusion weights (logistic
model of inclusion on X fitted separately in each arm), and compare with unweighted."""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
from common import load_raw, analytic, PRESPEC, HERE  # noqa: E402

always = ["age", "hba1c_val_0mo", "days_a1c_to_consent", "op_n_0mo", "ed_or_ip_n_0mo", "male",
          "hispanic", "race_nonhisp_black", "rural_site", "hba1c_cat_control", "drug_insulin_0mo",
          "drug_metformin_0mo", "drug_glp1_0mo", "drug_sulf_0mo", "drug_other_0mo", "ip_i_0mo", "ed_i_0mo"]
raw = load_raw()
out = []
for oc in ("hba1c", "weight"):
    df, F = analytic(raw, oc)
    r = raw.copy()
    r["inc"] = r.subject_id.isin(df.subject_id).astype(int)
    p = np.zeros(len(r))
    for a in (0, 1):
        m = r.treat_now == a
        Xs = StandardScaler().fit_transform(r.loc[m, always].astype(float))
        p[m.to_numpy()] = LogisticRegression(C=1.0, max_iter=5000).fit(Xs, r.loc[m, "inc"]).predict_proba(Xs)[:, 1]
    r["ipw"] = 1 / p
    d = r[r.inc == 1].merge(df[["subject_id", "Y"]], on="subject_id")
    for k in (1, 2, 3):
        d[f"st{k}"] = (d.stratum == k).astype(int)
    d["ipw"] = d.ipw / d.ipw.mean()
    print(f"\n{oc}: N={len(d)}, arm-specific IPW range {d.ipw.min():.2f}-{d.ipw.max():.2f}")
    for wlab, wts in (("unweighted", None), ("arm-specific IPW", d.ipw)):
        f0 = smf.wls("Y ~ treat_now + st1 + st2 + st3", d, weights=wts if wts is not None else 1.0).fit(cov_type="HC1")
        print(f"  {wlab:18s} overall effect {f0.params['treat_now']:+.3f} (p {f0.pvalues['treat_now']:.3f})")
        for v in PRESPEC:
            fm = smf.wls(f"Y ~ treat_now*{v} + st1 + st2 + st3", d, weights=wts if wts is not None else 1.0).fit(cov_type="HC1")
            out.append(dict(outcome=oc, weights=wlab, subgroup=v, inter=fm.params[f"treat_now:{v}"],
                            p=fm.pvalues[f"treat_now:{v}"]))
R = pd.DataFrame(out)
pd.set_option("display.width", 200)
print(R.pivot_table(index=["outcome", "subgroup"], columns="weights", values=["inter", "p"]).round(3).to_string())
R.to_csv(HERE / "e07_armspecific_ipw.csv", index=False)
