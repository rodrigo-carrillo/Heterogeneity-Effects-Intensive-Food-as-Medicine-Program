"""
Summarizes the Stata runs made by run_stata_checks.sh.

  python summarize_stata.py <work dir> <results dir> <project root>

Writes:
  stata_vs_python.csv   largest difference between Stata and code/replicate_doyle.py, by quantity
  imputed_cells.csv     Table 3 adjusted cells with imputation: published value vs the authors'
                        code under 40 sort seeds (kg and lb runs) and the proper per-imputation
                        baseline (10 sort seeds)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

work, out, root = map(Path, sys.argv[1:4])
LB_PER_KG = 2.20462

# ---- Stata vs the Python port (non-imputed numbers) ---------------------------
st = pd.read_csv(work / "default_a/output/T3_biometrics.csv")
py = pd.read_csv(root / "outputs/replication/table3.csv")
py = py[py.outcome != "weight (lb)"].assign(biometric=lambda d: d.outcome.str.replace(" (kg)", "", regex=False))
m = st.merge(py, on=["biometric", "month"], suffixes=("_st", "_py"))
pairs = [("ctl_mean_st", "ctl_mean_py"), ("ctl_sd_st", "ctl_sd_py"), ("trt_mean_st", "trt_mean_py"),
         ("trt_sd_st", "trt_sd_py"), ("unadj_diff", "unadj_est"), ("unadj_ci_low", "unadj_lo"),
         ("unadj_ci_high", "unadj_hi")]
rows = [("T3 " + a.replace("_st", ""), len(m), np.abs(m[a] - m[b]).max()) for a, b in pairs]
h = m[m.biometric == "hba1c"]
rows += [("T3 HbA1c " + a, len(h), np.abs(h[a] - h[b]).max())
         for a, b in (("adj_diff", "adj_est"), ("adj_ci_low", "adj_lo"), ("adj_ci_high", "adj_hi"))]
st4 = pd.read_csv(work / "default_a/output/T4_regressions.csv").rename(columns={"months": "month", "variable": "outcome"})
m4 = st4.merge(pd.read_csv(root / "outputs/replication/table4.csv"), on=["outcome", "month"])
sc = np.where(m4.outcome.str.endswith("_i") | m4.outcome.str.startswith("drug"), 100, 1)   # paper reports pp
rows += [(f"T4 {a} (relative)", len(m4), np.max(np.abs(m4[a] * sc - m4[b]) / np.maximum(1, np.abs(m4[b]))))
         for a, b in (("adj_diff", "est"), ("adj_ci_low", "lo"), ("adj_ci_high", "hi"))]
st2 = pd.read_csv(work / "default_a/output/T2_6mo_adjusted_diffs.csv")
m2 = st2.merge(pd.read_csv(root / "outputs/replication/table2_adjusted.csv").rename(columns={"outcome": "variable"}), on="variable")
rows += [("T2 " + a, len(m2), np.abs(m2[a] - m2[b]).max()) for a, b in (("adj_diff", "est"), ("ci_low", "lo"), ("ci_high", "hi"))]
pd.DataFrame(rows, columns=["quantity", "n_cells", "max_abs_diff"]).to_csv(out / "stata_vs_python.csv", index=False)

# ---- Imputed Table 3 cells -------------------------------------------------------
PUB = {(6, "ldl"): (5.86, -3.71, 15.42), (6, "hdl"): (-0.48, -2.85, 1.90), (6, "total_chol"): (4.94, -7.56, 17.44),
       (6, "trigly"): (3.63, -47.30, 54.57), (6, "weight"): (1.95, 0.07, 3.83), (6, "sys"): (-0.38, -4.16, 3.40),
       (6, "dias"): (-0.25, -2.28, 1.79), (12, "ldl"): (3.75, -6.20, 13.71), (12, "hdl"): (-0.33, -4.06, 3.40),
       (12, "total_chol"): (-0.51, -16.44, 15.43), (12, "trigly"): (-25.63, -83.77, 32.50),
       (12, "weight"): (1.88, -0.14, 3.89), (12, "sys"): (-2.49, -6.04, 1.05), (12, "dias"): (0.81, -1.20, 2.83)}


def load(prefix, n):
    return pd.concat([pd.read_csv(work / f"{prefix}{k}/output/T3_biometrics.csv").assign(k=k) for k in range(1, n + 1)])


ss, lb, pr = load("ss", 40), load("lb", 40), load("proper", 10)
dflt = st.set_index(["month", "biometric"])
rows = []
for (mo, v), (pe, pl, ph) in PUB.items():
    x = ss[(ss.month == mo) & (ss.biometric == v)]
    p = pr[(pr.month == mo) & (pr.biometric == v)]
    r = dict(month=mo, outcome=v, published=pe, published_lo=pl, published_hi=ph,
             default_run=dflt.loc[(mo, v)].adj_diff,
             sortseed_min=x.adj_diff.min(), sortseed_median=x.adj_diff.median(), sortseed_max=x.adj_diff.max(),
             sortseed_distinct=x.adj_diff.round(6).nunique(),
             published_within_range=bool(x.adj_diff.min() <= pe <= x.adj_diff.max()),
             runs_matching_published=int((x.adj_diff.round(2) == pe).sum()),
             runs_ci_excludes_0=int(((x.adj_ci_low > 0) | (x.adj_ci_high < 0)).sum()),
             proper=p.adj_diff.iloc[0], proper_lo=p.adj_ci_low.iloc[0], proper_hi=p.adj_ci_high.iloc[0],
             proper_distinct=p.adj_diff.round(6).nunique())
    if v == "weight":
        y = lb[(lb.month == mo) & (lb.biometric == v)]
        r |= dict(lb_run_kg_min=y.adj_diff.min() / LB_PER_KG, lb_run_kg_max=y.adj_diff.max() / LB_PER_KG,
                  lb_runs_matching_published=int((y.adj_diff / LB_PER_KG).round(2).eq(pe).sum()),
                  lb_unadj=y.unadj_diff.iloc[0], lb_unadj_lo=y.unadj_ci_low.iloc[0], lb_unadj_hi=y.unadj_ci_high.iloc[0])
    rows.append(r)
pd.DataFrame(rows).to_csv(out / "imputed_cells.csv", index=False)
