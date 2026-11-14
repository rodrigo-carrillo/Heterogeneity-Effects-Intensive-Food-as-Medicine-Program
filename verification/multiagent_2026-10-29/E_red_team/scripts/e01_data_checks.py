"""E01: data-level checks (dates, windows, missingness by arm, balance in ED subgroups)."""
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[4]
d = pd.read_csv(ROOT/"data/raw/dataverse_files/analysis_sample.csv")
pd.set_option("display.width", 200)
print("treat_now values:", d.treat_now.value_counts().to_dict())
# baseline dates relative to consent (positive = after consent -> post-randomization)
for v in ["hba1c","weight","bmi","sys","dias","ldl"]:
    x = d[f"{v}_date_0mo"] - d["consent_date"]
    print(f"{v}_date_0mo - consent: n={x.notna().sum()} min={x.min()} median={x.median()} max={x.max()} n_after_consent={(x>0).sum()}")
for v in ["hba1c","weight","sys"]:
    x = d[f"{v}_date_6mo"] - d["consent_date"]
    print(f"{v}_date_6mo - consent by arm:")
    print(x.groupby(d.treat_now).describe().round(0).to_string())
# weight conversions
r = d.weight_val_0mo / (d.weight_lb_val_0mo/2.20462)
print("posted kg / (lb/2.20462):", r.describe().round(4).to_dict())
# BMI vs weight consistency
d["wkg0"] = d.weight_lb_val_0mo/2.20462; d["wkg6"] = d.weight_lb_val_6mo/2.20462
d["dw"] = d.wkg6 - d.wkg0
print("bmi_date_0mo == weight_date_0mo:", (d.bmi_date_0mo==d.weight_date_0mo).mean())
print("weight change kg:", d.dw.describe().round(2).to_dict())
print("largest |dw|:"); print(d.loc[d.dw.abs().sort_values(ascending=False).index[:12], ["subject_id","treat_now","ed_i_0mo","wkg0","wkg6","dw","weight_date_0mo","weight_date_6mo","consent_date"]].round(1).to_string())
# missingness of 6-mo weight by arm x prior ED
d["has_w"] = d.dw.notna().astype(int)
print("\n6-mo weight change observed, by arm x prior ED:")
print(d.groupby(["ed_i_0mo","treat_now"]).has_w.agg(["size","mean"]).round(3).to_string())
import statsmodels.formula.api as smf
m = smf.logit("has_w ~ treat_now*ed_i_0mo", data=d).fit(disp=0)
print("logit observed ~ arm*ED: interaction p =", round(m.pvalues["treat_now:ed_i_0mo"],3))
for oc in ["hba1c","sys"]:
    d["has"] = (d[f"{oc}_val_6mo"]-d[f"{oc}_val_0mo"]).notna().astype(int)
    print(oc, d.groupby(["ed_i_0mo","treat_now"]).has.mean().round(3).to_dict())
# ED variables relationships
print("\ned_i_0mo vs ed_or_ip_n_0mo crosstab (count capped 5):")
print(pd.crosstab(d.ed_i_0mo, d.ed_or_ip_n_0mo.clip(upper=5)))
print("ed_i_0mo vs ip_i_0mo:"); print(pd.crosstab(d.ed_i_0mo, d.ip_i_0mo))
