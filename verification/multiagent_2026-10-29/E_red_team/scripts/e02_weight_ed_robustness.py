"""E02: robustness of weight x prior-year ED (cheap estimators; no forests)."""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
from common import load_raw, analytic, HERE  # noqa: E402

pd.set_option("display.width", 220)
raw = load_raw()
df, F = analytic(raw, "weight")
print("N", len(df), "treated", int(df.treat_now.sum()))
ST = "st1 + st2 + st3"
TERM = "treat_now:ed_i_0mo"
rows = []


def inter(formula, data, label, term=TERM):
    m = smf.ols(formula, data=data).fit(cov_type="HC1")
    b, se, p = m.params[term], m.bse[term], m.pvalues[term]
    print(f"{label:62s} n={int(m.nobs):4d} inter={b:+.2f} (SE {se:.2f}) p={p:.4f}")
    rows.append(dict(spec=label, n=int(m.nobs), est=b, se=se, p=p))


base = f"Y ~ treat_now*ed_i_0mo + {ST}"
inter(base, df, "A. OLS + strata (docs/04 spec)")
inter(base + " + weight_val_0mo", df, "B. + baseline weight")
inter(base + " + bmi_val_0mo", df, "C. + baseline BMI")
allc = " + ".join([c for c in F if c not in ("ed_i_0mo", "rural_site", "hba1c_cat_control")])
inter(base + " + " + allc, df, "D. + all X main effects")
inter(base + " + treat_now*weight_val_0mo", df, "E. + treat x baseline weight")
inter(base + " + treat_now*bmi_val_0mo", df, "F. + treat x BMI")
inter(base + " + treat_now*drug_insulin_0mo", df, "G. + treat x insulin")
inter(base + " + treat_now*ip_i_0mo", df, "H. + treat x prior IP")
inter(base + " + treat_now*hba1c_val_0mo", df, "I. + treat x baseline HbA1c")
inter(base + " + treat_now*op_n_0mo", df, "J. + treat x outpatient visits")
inter(base + " + treat_now*(drug_insulin_0mo + ip_i_0mo + bmi_val_0mo + age + male"
      " + hba1c_val_0mo + op_n_0mo)", df, "K. + treat x 7 other modifiers jointly")
inter(base + " + treat_now*(" + allc + ")", df, "L. + treat x ALL other X (fully interacted)")
d2 = df.copy()
d2["Ylev"] = d2.weight_val_6mo
inter(f"Ylev ~ treat_now*ed_i_0mo + weight_val_0mo + {ST}", d2, "M. ANCOVA on 6-mo level")
d2["Ypct"] = 100 * d2.Y / d2.weight_val_0mo
inter(f"Ypct ~ treat_now*ed_i_0mo + {ST}", d2, "N. % weight change (pp)")
for thr in (15, 10):
    inter(base, df[df.Y.abs() <= thr], f"O. exclude |change| > {thr} kg")
for q in (0.025, 0.05):
    lo, hi = df.Y.quantile([q, 1 - q])
    d3 = df.copy()
    d3["Y"] = d3.Y.clip(lo, hi)
    inter(base, d3, f"P. winsorize change at {q:.1%}/{1-q:.1%}")
X = df[["treat_now", "ed_i_0mo", "st1", "st2", "st3"]].astype(float).copy()
X["int"] = df.treat_now * df.ed_i_0mo
X = sm.add_constant(X)
hub = sm.RLM(df.Y, X, M=sm.robust.norms.HuberT()).fit()
print(f"{'Q. Huber robust regression':62s} inter={hub.params['int']:+.2f} (SE {hub.bse['int']:.2f}) p={hub.pvalues['int']:.4f}")
rows.append(dict(spec="Q. Huber robust regression", n=len(df), est=hub.params["int"],
                 se=hub.bse["int"], p=hub.pvalues["int"]))
qr = smf.quantreg(base, df).fit(q=0.5)
print(f"{'R. median regression':62s} inter={qr.params[TERM]:+.2f} p={qr.pvalues[TERM]:.4f}")
rows.append(dict(spec="R. median regression", n=len(df), est=qr.params[TERM], se=qr.bse[TERM], p=qr.pvalues[TERM]))
df["d6"] = df.weight_date_6mo - df.consent_date
df["d0"] = df.consent_date - df.weight_date_0mo
inter(base, df[df.d6 <= 210], "S. 6-mo weight <=210 d after consent")
inter(base, df[df.d6 <= 190], "S2. 6-mo weight <=190 d after consent")
inter(base, df[df.d0 <= 90], "T. baseline weight <=90 d before consent")
inter(base, df[df.d0 <= 30], "T2. baseline weight <=30 d before consent")
inter(base + " + d6 + d0", df, "U. + measurement timing")
for v, lab in [("ed_or_ip_i_0mo", "any ED or IP"), ("ip_i_0mo", "any IP")]:
    d4 = df.copy(); d4["ed_i_0mo"] = d4[v]
    inter(base, d4, f"V. modifier = {lab}")
d4 = df.copy(); d4["ed_i_0mo"] = (d4.ed_or_ip_n_0mo >= 2).astype(int)
inter(base, d4, "V. modifier = >=2 ED/IP visits")
d4 = df.copy(); d4["ed_i_0mo"] = np.log1p(d4.ed_or_ip_n_0mo)
inter(base, d4, "V. modifier = log(1+ED/IP count), continuous")
d4 = df[(df.ed_i_0mo == 0) | (df.ed_or_ip_n_0mo >= 2)]
inter(base, d4, "V. ED=0 vs ED with >=2 ED/IP visits (drop single-visit)")
r = raw.copy()
r["Y"] = r.weight_val_6mo - r.weight_val_0mo
r = r[r.Y.notna()].copy()
for k in (1, 2, 3):
    r[f"st{k}"] = (r.stratum == k).astype(int)
inter(base, r, "W. all with weight change (N=395, no CC on X)")
for s in (0, 1):
    inter(f"Y ~ treat_now*ed_i_0mo + hba1c_cat_control", df[df.rural_site == s], f"X. site rural={s}")
for s in (0, 1):
    inter(f"Y ~ treat_now*ed_i_0mo + rural_site", df[df.hba1c_cat_control == s], f"X. HbA1c stratum={s}")
for s in (0, 1):
    inter(base, df[df.male == s], f"X. male={s}")
res = pd.DataFrame(rows)
res.to_csv(HERE / "e02_weight_ed_specs.csv", index=False)

print("\nMean weight change by prior ED x arm (kg):")
print(df.groupby(["ed_i_0mo", "treat_now"]).Y.agg(["size", "mean", "std", "median"]).round(2).to_string())

# Influence
fit = smf.ols(base, df).fit()
k = list(fit.params.index).index(TERM)
dfb = fit.get_influence().dfbetas[:, k]
order = np.argsort(-np.abs(dfb))
print("\nTop-8 |DFBETAS| for the interaction:")
print(df.iloc[order[:8]][["subject_id", "treat_now", "ed_i_0mo", "weight_val_0mo", "Y"]]
      .assign(dfbeta=dfb[order[:8]]).round(2).to_string())
loo = np.array([smf.ols(base, df.drop(index=i)).fit(cov_type="HC1").pvalues[TERM] for i in range(len(df))])
print(f"Leave-one-out interaction p: max {loo.max():.4f}, median {np.median(loo):.4f}")
d5 = df.copy(); removed = 0; hit_holm = hit_raw = None
sign = np.sign(fit.params[TERM])
while removed < 80:
    mm = smf.ols(base, d5).fit()
    p = smf.ols(base, d5).fit(cov_type="HC1").pvalues[TERM]
    if hit_holm is None and p > 0.05 / 7:
        hit_holm = removed
    if p > 0.05:
        hit_raw = removed
        break
    db = mm.get_influence().dfbetas[:, k]
    d5 = d5.drop(index=d5.index[int(np.argmax(db * sign))]); removed += 1
print(f"Greedy removal of most supportive patients: p > 0.05/7 after {hit_holm}; p > 0.05 after {hit_raw} of {len(df)}")
# compare: same exercise for a "typical" null subgroup to calibrate fragility is not meaningful; report as is.

print("\nBaseline balance (SMD treat vs control) within prior-ED subgroups (|SMD|>0.2 shown):")
bal = []
for v in F + ["ldl_val_0mo"]:
    row = {"var": v}
    for e in (0, 1):
        s = df[df.ed_i_0mo == e]
        t, c = s.loc[s.treat_now == 1, v].dropna(), s.loc[s.treat_now == 0, v].dropna()
        sd = np.sqrt((t.var() + c.var()) / 2)
        row[f"smd_ed{e}"] = (t.mean() - c.mean()) / sd if sd > 0 else 0
    bal.append(row)
bal = pd.DataFrame(bal)
print(bal[(bal.smd_ed0.abs() > 0.2) | (bal.smd_ed1.abs() > 0.2)].round(3).to_string(index=False))
print("max |SMD| ED=0:", bal.smd_ed0.abs().max().round(3), " ED=1:", bal.smd_ed1.abs().max().round(3))
bal.to_csv(HERE / "e02_balance_within_ed.csv", index=False)
print("\nPrior ED=1 vs 0, baseline means:")
print(df.groupby("ed_i_0mo")[F].mean().T.round(2).to_string())
