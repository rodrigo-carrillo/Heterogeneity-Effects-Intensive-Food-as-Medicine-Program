"""Independent verifier: claims 1, 2, 5 (replication of published FE models, weight units, ED var).
Written from scratch, blind to the project's code/ and outputs/."""
import numpy as np, pandas as pd
from scipy import stats

D = pd.read_csv("data/raw/dataverse_files/analysis_sample.csv")


def ols(y, X):
    XtX_inv = np.linalg.inv(X.T @ X)
    b = XtX_inv @ X.T @ y
    e = y - X @ b
    return b, e, XtX_inv


def fe_xtreg(df, var, month=6, timevar=True, yname=None, basename=None):
    """Stata-style xtreg y post treat post_treat post_baseline [time] ..., fe vce(robust).
    Within-transformation on the stacked 2-period panel, cluster-robust SE with Stata's
    small-sample factor (G/(G-1))*((N-1)/(N-K)), t with G-1 df."""
    y0n = yname.format(0) if yname else f"{var}_val_0mo"
    y1n = yname.format(month) if yname else f"{var}_val_{month}mo"
    d = df[df[y1n].notna()].copy()
    rows = []
    for per, yn in [(0, y0n), (1, y1n)]:
        r = pd.DataFrame({"id": d.subject_id, "y": d[yn], "post": per,
                          "treat": d.treat_now})
        r["post_treat"] = per * d.treat_now.values
        r["post_base"] = per * d[y0n].values
        if timevar:
            dn = f"{var}_date_{0 if per == 0 else month}mo"
            r["time"] = (d[dn] - d.consent_date).values
        rows.append(r)
    p = pd.concat(rows).dropna()  # complete case (no MI)
    # keep only ids with both periods
    p = p[p.groupby("id").id.transform("size") == 2]
    xs = ["post", "post_treat", "post_base"] + (["time"] if timevar else [])
    W = p[["y"] + xs] - p.groupby("id")[["y"] + xs].transform("mean")
    X = W[xs].values; y = W["y"].values
    b, e, Ainv = ols(y, X)
    G = p.id.nunique(); N = len(p); K = len(xs) + 1  # + constant
    meat = np.zeros((len(xs), len(xs)))
    for g, idx in p.groupby("id").indices.items():
        s = X[idx].T @ e[idx]; meat += np.outer(s, s)
    V = Ainv @ meat @ Ainv * (G / (G - 1)) * ((N - 1) / (N - K))
    j = xs.index("post_treat"); se = np.sqrt(V[j, j]); tq = stats.t.ppf(.975, G - 1)
    return dict(b=b[j], lo=b[j] - tq * se, hi=b[j] + tq * se, G=G,
                p=2 * stats.t.sf(abs(b[j] / se), G - 1))


def fd(df, var, month=6, timevar=True, yname=None, scale=1.0):
    """First difference: dY ~ 1 + treat + Y0 + dtime, HC1."""
    y0n = yname.format(0) if yname else f"{var}_val_0mo"
    y1n = yname.format(month) if yname else f"{var}_val_{month}mo"
    d = df[df[y1n].notna() & df[y0n].notna()].copy()
    cols = [np.ones(len(d)), d.treat_now.values, d[y0n].values * scale]
    if timevar:
        cols.append(((d[f"{var}_date_{month}mo"] - d.consent_date) -
                     (d[f"{var}_date_0mo"] - d.consent_date)).values)
    X = np.column_stack(cols); y = ((d[y1n] - d[y0n]) * scale).values
    m = ~np.isnan(X).any(1); X, y = X[m], y[m]
    b, e, Ainv = ols(y, X); n, k = X.shape
    V = Ainv @ (X.T * e**2) @ X @ Ainv * n / (n - k)
    se = np.sqrt(V[1, 1]); tq = stats.t.ppf(.975, n - k)
    return dict(b=b[1], lo=b[1] - tq * se, hi=b[1] + tq * se, n=n)


def fmt(r, s=1.0):
    return f"{r['b']*s:.3f} ({r['lo']*s:.3f} to {r['hi']*s:.3f})"


print("== Claim 1 ==")
r = fe_xtreg(D, "hba1c"); print("HbA1c FE xtreg :", fmt(r), "G=", r["G"], "p=%.3f" % r["p"])
print("HbA1c FD       :", fmt(fd(D, "hba1c")), "n=", fd(D, "hba1c")["n"])
r = fe_xtreg(D, "ed_or_ip_i", timevar=False, yname="ed_or_ip_i_{}mo")
print("ED/IP FE xtreg (x100):", fmt(r, 100), "G=", r["G"])
r = fd(D, "ed_or_ip_i", timevar=False, yname="ed_or_ip_i_{}mo")
print("ED/IP FD (x100)      :", fmt(r, 100), "n=", r["n"])
# Also check outpatient visits count (published 1.00 (0.28 to 1.72)) as extra replication
print("op_n FE:", fmt(fe_xtreg(D, "op_n", timevar=False, yname="op_n_{}mo")))

print("\n== Claim 2: weight units ==")
d = D
print("weight_val_0mo mean", d.weight_val_0mo.mean(), " weight_lb_val_0mo mean", d.weight_lb_val_0mo.mean())
print("ratio lb/val:", (d.weight_lb_val_0mo / d.weight_val_0mo).describe())
h = np.sqrt(d.weight_val_0mo / d.bmi_val_0mo)
print("implied height (m) if kg:", h.describe())
for g in [0, 1]:
    s = d.loc[d.treat_now == g, "weight_val_6mo"]
    print(f"6mo weight_val treat={g}: mean {s.mean():.2f} sd {s.std():.2f} n={s.notna().sum()}")
# unadjusted (regress y treat, robust)
for nm in ["weight_val_6mo", "weight_lb_val_6mo"]:
    dd = d[d[nm].notna()]
    X = np.column_stack([np.ones(len(dd)), dd.treat_now]); y = dd[nm].values
    b, e, A = ols(y, X); n = len(y); V = A @ (X.T * e**2) @ X @ A * n / (n - 2)
    se = np.sqrt(V[1, 1]); tq = stats.t.ppf(.975, n - 2)
    print(f"unadjusted {nm}: {b[1]:.2f} ({b[1]-tq*se:.2f} to {b[1]+tq*se:.2f})")
rk = fe_xtreg(D, "weight"); print("adj FE complete-case kg :", fmt(rk), "G=", rk["G"], "p=%.3f" % rk["p"])
rl = fe_xtreg(D, "weight", yname="weight_lb_val_{}mo"); print("adj FE complete-case lb :", fmt(rl), "G=", rl["G"], "p=%.3f" % rl["p"])
print("lb result converted to kg:", fmt(rl, 1 / 2.20462))
print("published 1.95 lb -> kg:", 1.95 / 2.20462, 0.07 / 2.20462, 3.83 / 2.20462)

print("\n== Claim 5: ed_i_0mo counts ==")
print("overall", int(D.ed_i_0mo.sum()), D.groupby("treat_now").ed_i_0mo.sum().to_dict())
lab = D[D.hba1c_val_6mo.notna()]
print("lab sample", int(lab.ed_i_0mo.sum()), lab.groupby("treat_now").ed_i_0mo.sum().to_dict())
print("ip_i_0mo", int(D.ip_i_0mo.sum()), D.groupby("treat_now").ip_i_0mo.sum().to_dict())
print("consistency: ed_i==0 & ed_or_ip_i==0 implies ...",
      pd.crosstab(D.ed_i_0mo, D.ed_or_ip_i_0mo), pd.crosstab(D.ed_or_ip_i_0mo, D.ip_i_0mo))
