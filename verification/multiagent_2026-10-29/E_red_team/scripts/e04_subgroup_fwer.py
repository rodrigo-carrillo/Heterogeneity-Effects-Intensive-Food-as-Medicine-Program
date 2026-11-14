"""E04: how often does a 'weight x prior ED'-like hit arise by chance in the subgroup search
actually performed (7 Doyle + 3 added subgroups, 4 outcomes, 2 covariate sets)?

Treatment is permuted within the 4 randomization strata on all 465 patients (same permutation
for every outcome, which preserves the real correlation between outcomes), so the null is
'no effect of any kind'. Interaction p: OLS of Y on W x S + stratum dummies, HC1 SE (docs/04
shows this gives nearly the same p as the pipeline's AIPW: 0.0025 vs 0.0019).
Also: randomization p of the weight x ED interaction under a CONSTANT-effect null (Y - tau*W),
which is the null that matters for an interaction (the docs' permutation test used the sharp
zero-effect null)."""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
from common import load_raw, analytic, PRESPEC, HERE  # noqa: E402

B = 2000
raw = load_raw()
raw["age_ge65"] = (raw.age >= 65).astype(int)
raw["nonhisp_white"] = ((raw.hispanic == 0) & (raw.race == "White")).astype(int)
raw["a1c_test_gt90d"] = (raw.days_a1c_to_consent > 90).astype(int)
ADDED = ["age_ge65", "nonhisp_white", "a1c_test_gt90d"]
RUNS = [(oc, ldl) for oc in ("hba1c", "weight", "sbp", "ed_ip") for ldl in (False, True)]
samples = {}
for oc, ldl in RUNS:
    df, F = analytic(raw, oc, ldl)
    keep = raw.index[raw.subject_id.isin(df.subject_id)]
    samples[(oc, ldl)] = keep.to_numpy()


def inter_p(y, w, s, st):
    """Vectorised OLS interaction p (HC1) for y ~ 1 + w + s + w*s + stratum dummies."""
    X = np.column_stack([np.ones(len(y)), w, s, w * s, st])
    XtX_inv = np.linalg.pinv(X.T @ X)
    b = XtX_inv @ X.T @ y
    r = y - X @ b
    n, k = X.shape
    meat = (X * r[:, None] ** 2).T @ X
    V = XtX_inv @ meat @ XtX_inv * n / (n - k)
    z = b[3] / np.sqrt(V[3, 3])
    from scipy.stats import norm
    return 2 * norm.sf(abs(z)), b[3]


Yall = {}
for oc in ("hba1c", "weight", "sbp", "ed_ip"):
    d, _ = analytic(raw.assign(), oc)  # just to get Y definition
r = raw.copy()
r["Y_hba1c"] = r.hba1c_val_6mo - r.hba1c_val_0mo
r["Y_weight"] = r.weight_val_6mo - r.weight_val_0mo
r["Y_sbp"] = r.sys_val_6mo - r.sys_val_0mo
r["Y_ed_ip"] = 100.0 * r.ed_or_ip_i_6mo
strata = r.stratum.to_numpy()
STD = np.column_stack([(strata == k).astype(float) for k in (1, 2, 3)])
W0 = r.treat_now.to_numpy().astype(float)


def all_ps(W):
    out = {}
    for (oc, ldl), idx in samples.items():
        y = r[f"Y_{oc}"].to_numpy()[idx]
        for fam, vars_ in (("doyle", PRESPEC), ("added", ADDED)):
            ps = [inter_p(y, W[idx], r[v].to_numpy()[idx].astype(float), STD[idx])[0] for v in vars_]
            out[(oc, ldl, fam)] = np.array(ps)
    return out


obs = all_ps(W0)
print("Observed (OLS) Doyle-family p by run; Holm min:")
for k, v in obs.items():
    if k[2] == "doyle":
        print(k, np.round(v, 4), "Holm min", multipletests(v, method="holm")[1].min().round(4))

rng = np.random.default_rng(2026)
rows = []
for b in range(B):
    Wp = W0.copy()
    for s in np.unique(strata):
        idx = np.where(strata == s)[0]
        Wp[idx] = rng.permutation(W0[idx])
    ps = all_ps(Wp)
    row = {}
    for k, v in ps.items():
        row["holm_" + "_".join(map(str, k))] = multipletests(v, method="holm")[1].min()
        row["raw_" + "_".join(map(str, k))] = v.min()
    rows.append(row)
R = pd.DataFrame(rows)
R.to_csv(HERE / "e04_perm_minp.csv", index=False)

print(f"\n{B} within-stratum permutations (sharp null), OLS interaction p")
for oc in ("hba1c", "weight", "sbp", "ed_ip"):
    print(f"  {oc:7s} main: P(Doyle-family min Holm < .05) = {(R[f'holm_{oc}_False_doyle'] < .05).mean():.3f}")
main_doyle = R[[f"holm_{oc}_False_doyle" for oc in ("hba1c", "weight", "sbp", "ed_ip")]].min(axis=1)
print(f"P(any of 4 outcomes, main runs, Doyle family Holm < .05) = {(main_doyle < .05).mean():.3f}")
print(f"P(any of 4 outcomes, main runs, Doyle family Holm <= .0134 [observed]) = {(main_doyle <= .0134).mean():.3f}")
allcols = [c for c in R.columns if c.startswith("holm_")]
print(f"P(any Holm < .05 in any family x 8 runs) = {(R[allcols].min(axis=1) < .05).mean():.3f}")
rawcols = [c for c in R.columns if c.startswith("raw_")]
print(f"P(any raw subgroup p < .05 among 10 subgroups x 8 runs) = {(R[rawcols].min(axis=1) < .05).mean():.3f}")
print(f"P(any raw subgroup p <= .0025 among 10 subgroups x 8 runs) = {(R[rawcols].min(axis=1) <= .0025).mean():.3f}")

# Constant-effect randomization test for weight x ED (Ding-Feller-Miratrix style plug-in)
idx = samples[("weight", False)]
y = r["Y_weight"].to_numpy()[idx]; w = W0[idx]; s = r["ed_i_0mo"].to_numpy()[idx].astype(float); st = STD[idx]
p_obs, b_obs = inter_p(y, w, s, st)
X = np.column_stack([np.ones(len(y)), w, st]); tau0 = np.linalg.lstsq(X, y, rcond=None)[0][1]
res = {}
for tau in (tau0, 0.0, -0.5, 2.0):  # sensitivity of the plug-in over a range covering the ATE CI
    y0 = y - tau * w
    cnt = 0
    rng2 = np.random.default_rng(7)
    for _ in range(5000):
        wp = w.copy()
        for k in np.unique(strata[idx]):
            ii = np.where(strata[idx] == k)[0]
            wp[ii] = rng2.permutation(w[ii])
        cnt += abs(inter_p(y0 + tau * wp, wp, s, st)[1]) >= abs(b_obs)
    res[round(tau, 3)] = (cnt + 1) / 5001
print(f"\nWeight x ED interaction b = {b_obs:+.2f}, HC1 p = {p_obs:.4f}")
print("Randomization p under constant effect tau (5000 within-stratum perms, within analytic sample):", res)
