"""
================================================================================
Heterogeneous Treatment Effects (HTE) in the Food-as-Medicine RCT
(Doyle et al., JAMA Intern Med 2024;184(2):154-163; NCT03718832)
================================================================================
Version 3 (2026-10-21). Adds secondary outcomes (weight, SBP, ED/inpatient use)
and runs with baseline LDL. v2 = HbA1c only; v1 superseded
(code/archive/v1/, docs/03_v1_review.md). Complete-case throughout (PI decision
2026-10-21): every patient analysed has the outcome and every candidate modifier
observed; nothing is imputed.

Question
--------
The trial found no average effect of an intensive food-as-medicine program on
6-month HbA1c (adjusted difference -0.10 pp, 95% CI -0.46 to 0.25). Do some
patients, defined by baseline characteristics, benefit (or are harmed) even
though the average effect is null?  Methods follow the ACCORD causal-forest /
causal-tree papers (Jiao 2024 Diabetes Obes Metab, doi:10.1111/dom.15453;
Kianmehr 2022 J Diabetes Complications, doi:10.1016/j.jdiacomp.2022.108287),
adapted to a small (N=349) wait-list trial. The same pipeline is applied to
the trial's 6-month secondary outcomes (--outcome).

Design constraints that drive every choice below
------------------------------------------------
1. Wait-list design: controls start the program at ~6 months, so the only
   randomized treatment-vs-usual-care contrast is the 6-MONTH outcome.
2. N = 349 with a 6-month HbA1c (~30x smaller than ACCORD). Everything is
   exploratory / hypothesis-generating.
3. Randomization was 1:1, stratified by site (rural/urban) x HbA1c (>=9.5% vs
   <9.5%). The propensity score is KNOWN by design, so it is not estimated
   with machine learning; we use the within-stratum treated fraction.
4. Every heterogeneity test is computed on OUT-OF-FOLD (cross-fitted) CATE
   predictions. v1 scored in-sample predictions against the same outcomes,
   which manufactured a spurious "significant" BLP slope (+13, p<0.001).

Pipeline
--------
  §1  Load data, define the analytic sample, outcome Y, treatment W, strata
  §2  Reproduce the published primary result (sanity check)
  §3  Baseline balance in the analytic sample
  §4  Candidate effect modifiers X (baseline only)
  §5  Average treatment effect (AIPW + causal forest)
  §6  Cross-fitted CATEs and doubly-robust (AIPW) scores, repeated K-fold
  §7  Omnibus heterogeneity tests: BLP calibration, GATES, RATE (AUTOC)
  §8  Which variables? Forest importance + BLP of DR scores on covariates
  §9  Pre-specified subgroups (trial strata and Doyle's subgroups)
  §10 ACCORD-style subgroups: tertiles + honest (sample-split) causal tree
  §11 Sensitivity analyses (outcome level, LDL, attrition IPW, stale A1c)
  §12 Figures
  §13 Write results (CSV, JSON, Markdown summary)

Outcomes (all at 6 months, the only randomized contrast; see OUTCOMES in §1)
  hba1c   HbA1c(6 mo) - HbA1c(baseline), pp          [primary]
  weight  weight(6 mo) - weight(baseline), kg
  sbp     SBP(6 mo) - SBP(baseline), mm Hg
  ed_ip   any ED visit or inpatient admission in the 6 months after consent,
          scaled 0/100 so effects are risk differences in percentage points
Sign convention: NEGATIVE effect = the outcome is lower with the program than
with usual care (less HbA1c, weight, SBP or ED/IP use).

How to run
----------
  python code/hte_food_as_medicine.py                    # HbA1c, full run (~10 min)
  python code/hte_food_as_medicine.py --fast             # quick smoke test (~1 min)
  python code/hte_food_as_medicine.py --outcome weight   # secondary outcome
  python code/hte_food_as_medicine.py --ldl
        # baseline LDL added to the modifiers (patients with LDL measured)
  python code/run_all.py                                 # every outcome x covariate set
Outputs go to outputs/v3/<run id>/ (outputs/v3_fast/ with --fast), where the
run id is e.g. "hba1c" or "weight_ldl_cc".
The notebook notebooks/hte_food_as_medicine.ipynb is generated from this file
by code/build_notebook.py; edit this file, not the notebook.
================================================================================
"""

# %% [markdown]
# ## §0 Setup and configuration

# %%
import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                    # file output only; works headless
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeRegressor
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# --- Command-line options (ignored gracefully inside Jupyter) ----------------
parser = argparse.ArgumentParser(description="HTE analysis, food-as-medicine RCT")
parser.add_argument("--fast", action="store_true",
                    help="fewer trees/repeats for a quick smoke test")
parser.add_argument("--outcome", default="hba1c", choices=["hba1c", "weight", "sbp", "ed_ip"],
                    help="6-month outcome to analyse (default: hba1c)")
parser.add_argument("--ldl", action="store_true",
                    help="include baseline LDL (the single lipid marker) in the MAIN modifier set")
ARGS, _ = parser.parse_known_args()

SEED = 20260621
FAST = ARGS.fast
N_TREES = 500 if FAST else 4000          # trees in each causal forest
N_TREES_NUIS = 200 if FAST else 500      # trees in each nuisance forest
K_FOLDS = 5                              # cross-fitting folds
N_REPEATS = 2 if FAST else 5             # repeated cross-fitting (median-aggregated)
N_REPEATS_SENS = 1 if FAST else 3        # repeats for sensitivity analyses
N_BOOT = 200 if FAST else 2000           # bootstrap draws (RATE)
N_STABILITY = 20 if FAST else 200        # honest-tree split-stability repeats
MIN_LEAF_CF = 10                         # min leaf size, causal forest
MIN_LEAF_TREE = 30                       # min leaf size, interpretable tree
N_TOP_SPLITS = 5                         # honest single-variable splits (Jiao 2024: top 5)
PC_SD = (0.25, 0.5, 1.0)                 # positive-control contrasts, in SDs of Y


def find_project_root() -> Path:
    """Walk up from this file (or the notebook's cwd) to the project root."""
    start = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
    for p in [start, *start.parents]:
        if (p / "data" / "raw" / "dataverse_files" / "analysis_sample.csv").exists():
            return p
    raise FileNotFoundError("Could not find data/raw/dataverse_files/analysis_sample.csv")


ROOT = find_project_root()
sys.path.insert(0, str(ROOT / "code"))
import hte_lib                           # model fitting + crash-isolated worker (see its docstring)
from hte_lib import aipw_scores, make_forest

DATA_CSV = ROOT / "data" / "raw" / "dataverse_files" / "analysis_sample.csv"
RUN_ID = ARGS.outcome + ("_ldl_cc" if ARGS.ldl else "")     # "_cc": complete case, kept for continuity
OUT = ROOT / "outputs" / ("v3_fast" if FAST else "v3") / RUN_ID
OUT.mkdir(parents=True, exist_ok=True)


class Tee:
    """Mirror console output to outputs/.../console_log.txt."""
    def __init__(self, path):
        self.file = open(path, "w", encoding="utf-8")
        self.stdout = sys.stdout
    def write(self, s):
        self.stdout.write(s); self.file.write(s)
        self.flush()                     # keep the log complete even if a run dies
    def flush(self):
        self.stdout.flush(); self.file.flush()


if not isinstance(sys.stdout, Tee):
    sys.stdout = Tee(OUT / "console_log.txt")

RESULTS = {"config": {"run_id": RUN_ID, "outcome": ARGS.outcome, "fast": FAST,
                      "ldl_in_main": ARGS.ldl, "complete_case": True,
                      "seed": SEED, "n_trees": N_TREES, "k_folds": K_FOLDS,
                      "n_repeats": N_REPEATS}}


def header(title):
    print("\n" + "=" * 78 + f"\n{title}\n" + "=" * 78)


def fmt_ci(est, lo, hi, d=2):
    return f"{est:+.{d}f} ({lo:+.{d}f} to {hi:+.{d}f})"


# Plot styling: one accent hue for single-series charts, recessive axes/grid.
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
ACCENT, ACCENT2 = "#2a78d6", "#eb6834"
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
})

# %% [markdown]
# ## §1 Data, analytic sample, outcome, treatment, strata
# OUTCOMES holds everything outcome-specific. `val`/`date` are column stems
# (e.g. weight_val_0mo, weight_val_6mo). kind "change" = 6-mo minus baseline;
# kind "binary" = 0/100 indicator. `paper_*` = Doyle 2024 Table 3/4 (6 mo):
# adjusted = longitudinal FE model; unadjusted = difference in 6-mo level.
# Candidate effect modifiers (baseline only) are defined here because the
# complete-case sample is selected on them.

# %%
header("§1 Data and analytic sample")
OUTCOMES = {
    "hba1c": dict(label="HbA1c", desc="HbA1c change", unit="pp", kind="change",
                  val="hba1c_val", date="hba1c_date", extra_x=[],
                  better="negative = program lowered HbA1c more",
                  paper_adj=(-0.10, -0.46, 0.25), paper_unadj=(-0.03, -0.42, 0.35), paper_n=349),
    "weight": dict(label="Weight", desc="weight change", unit="kg", kind="change",
                   val="weight_val", date="weight_date", extra_x=["weight_val_0mo"],
                   better="negative = less weight gain with program",
                   # Table 3: unadjusted difference is in lb (labelled kg); the adjusted
                   # +1.95 kg rests on 5 imputed baselines (see §2, outputs/replication/)
                   paper_adj=(1.95, 0.07, 3.83), paper_unadj=(-8.14, -20.57, 4.30), paper_n=400,
                   paper_note="unadjusted published in lb; adjusted kg with 5 imputed baselines"),
    "sbp": dict(label="Systolic BP", desc="SBP change", unit="mm Hg", kind="change",
                val="sys_val", date="sys_date", extra_x=[],
                better="negative = program lowered SBP more",
                paper_adj=(-0.38, -4.16, 3.40), paper_unadj=(-1.87, -5.85, 2.11), paper_n=378),
    "ed_ip": dict(label="Any ED/inpatient use", desc="any ED/IP use", unit="pp", kind="binary",
                  val="ed_or_ip_i", date=None, extra_x=[],
                  better="negative = fewer patients with ED/IP use",
                  # paper reports no CI for the unadjusted difference (30.43% vs 35.74%)
                  paper_adj=(-7.13, -15.30, 1.04), paper_unadj=(-5.31, np.nan, np.nan), paper_n=465),
}
OC = OUTCOMES[ARGS.outcome]
UNIT = OC["unit"]
EFFECT_LAB = f"Effect on {OC['desc']}, {UNIT}"

CONTINUOUS = ["age", "hba1c_val_0mo", "days_a1c_to_consent", "bmi_val_0mo",
              "sys_val_0mo", "dias_val_0mo", "op_n_0mo", "ed_or_ip_n_0mo"] + OC["extra_x"]
BINARY = ["male", "hispanic", "race_nonhisp_black", "rural_site", "hba1c_cat_control",
          "drug_insulin_0mo", "drug_metformin_0mo", "drug_glp1_0mo", "drug_sulf_0mo",
          "drug_other_0mo", "drug_any_db_0mo", "ip_i_0mo", "ed_i_0mo"]
# One lipid marker only (PI decision 2026-10-21): the one with the fewest missing
# baseline values. Missing of 465: LDL 92, TG 101, total cholesterol 102, HDL 104
# (checked in §1 below); LDL also has the fewest missing in every outcome sample.
LIPID_CANDIDATES = ["ldl_val_0mo", "trigly_val_0mo", "total_chol_val_0mo", "hdl_val_0mo"]
LIPIDS = ["ldl_val_0mo"]
MAIN_COLS = CONTINUOUS + BINARY + (LIPIDS if ARGS.ldl else [])

LABELS = {
    "age": "Age (y)", "hba1c_val_0mo": "Baseline HbA1c (%)",
    "days_a1c_to_consent": "Days, baseline A1c test to consent", "bmi_val_0mo": "BMI",
    "sys_val_0mo": "Systolic BP", "dias_val_0mo": "Diastolic BP",
    "op_n_0mo": "Prior-year outpatient visits", "ed_or_ip_n_0mo": "Prior-year ED/IP visits",
    "male": "Male", "hispanic": "Hispanic", "race_nonhisp_black": "Non-Hispanic Black",
    "rural_site": "Rural site", "hba1c_cat_control": "HbA1c stratum >=9.5%",
    "drug_insulin_0mo": "Insulin", "drug_metformin_0mo": "Metformin",
    "drug_glp1_0mo": "GLP-1 RA", "drug_sulf_0mo": "Sulfonylurea",
    "drug_other_0mo": "Other glucose-lowering", "drug_any_db_0mo": "Any diabetes drug",
    "ip_i_0mo": "Any prior-year inpatient", "ed_i_0mo": "Any prior-year ED",
    "hdl_val_0mo": "HDL", "ldl_val_0mo": "LDL", "total_chol_val_0mo": "Total cholesterol",
    "trigly_val_0mo": "Triglycerides", "weight_val_0mo": "Baseline weight (kg)",
}

raw = pd.read_csv(DATA_CSV)
print(f"Loaded {DATA_CSV.relative_to(ROOT)}: {raw.shape[0]} patients x {raw.shape[1]} columns")
print(f"Run: {RUN_ID}  (outcome {OC['desc']} at 6 mo; LDL in X: {ARGS.ldl}; complete case, no imputation)")

# Data cleaning: one baseline triglyceride is the placeholder 9999999 mg/dL.
bad_tg = raw["trigly_val_0mo"] >= 9_999_999
raw.loc[bad_tg, "trigly_val_0mo"] = np.nan
print(f"Set {int(bad_tg.sum())} placeholder triglyceride value(s) (9999999) to missing")
# The posted kg weights are lb / 2.2222 (factor 1/0.45), 0.8% too low; rebuild them from
# the (undocumented) lb columns with the exact factor. Effects scale; p-values do not change.
for m in (0, 6, 12):
    raw[f"weight_val_{m}mo"] = raw[f"weight_lb_val_{m}mo"] / 2.20462
print("Weight in kg recomputed as weight_lb_val / 2.20462 (posted kg column used 2.2222)")
lipid_miss = raw[LIPID_CANDIDATES].isna().sum().sort_values()
print(f"Baseline lipids missing (of {len(raw)}): {lipid_miss.to_dict()}; marker used: {LIPIDS[0]}")
assert lipid_miss.index[0] == LIPIDS[0], "the least-missing lipid is no longer LDL"

# Derived baseline variables (dates are stored as Stata day numbers).
# Days from the qualifying (baseline) HbA1c test to consent/randomization.
# A stale baseline value carries more regression-to-the-mean, so we use it as
# a candidate modifier and adjust for it (the paper adjusts for test timing).
raw["days_a1c_to_consent"] = raw["consent_date"] - raw["hba1c_date_0mo"]
# Randomization stratum: site x HbA1c category (0..3)
raw["stratum"] = 2 * raw["rural_site"] + raw["hba1c_cat_control"]

# Outcome. Y is the analysis outcome; Y_level the 6-mo value (sensitivity S1).
if OC["kind"] == "change":
    raw["Y"] = raw[f"{OC['val']}_6mo"] - raw[f"{OC['val']}_0mo"]
    raw["Y_level"] = raw[f"{OC['val']}_6mo"]
    raw["days_base_to_consent"] = raw["consent_date"] - raw[f"{OC['date']}_0mo"]
else:
    raw["Y"] = 100.0 * raw[f"{OC['val']}_6mo"]
    raw["Y_level"] = raw["ed_or_ip_n_6mo"]          # S1: count of ED/IP visits instead
# Analytic sample (complete case): outcome (and its baseline) observed, and every
# candidate modifier in X observed (with --ldl, baseline LDL too). No imputation.
raw["in_sample"] = (raw["Y"].notna() & raw[MAIN_COLS].notna().all(axis=1)).astype(int)

df = raw[raw["in_sample"] == 1].copy().reset_index(drop=True)
W = df["treat_now"].astype(int).to_numpy()
Y = df["Y"].to_numpy()
print(f"Analytic sample: N={len(df)} (treatment {W.sum()}, control {(1 - W).sum()}); "
      f"paper's 6-mo sample for this outcome: {OC['paper_n']}")
print(f"Excluded from the {len(raw)} randomized-and-retained: {1 - raw['in_sample'].mean():.1%} overall; "
      f"by arm {(1 - raw.groupby('treat_now')['in_sample'].mean()).round(3).to_dict()}")

# Known design propensity: treated fraction within each randomization stratum.
strat_tab = df.groupby("stratum")["treat_now"].agg(["size", "mean"])
strat_tab.index = ["urban, A1c<9.5", "urban, A1c>=9.5", "rural, A1c<9.5", "rural, A1c>=9.5"]
print("\nRandomization strata in the analytic sample (n, treated fraction):")
print(strat_tab.round(3).to_string())
E_HAT = df.groupby("stratum")["treat_now"].transform("mean").to_numpy()

for arm, lab in [(0, "Control"), (1, "Treatment")]:
    s = df.loc[W == arm, "Y"]
    print(f"{lab:9s}: mean {OC['desc']} {s.mean():+.3f} (SD {s.std():.2f}) {UNIT}")
print(f"Unadjusted difference in Y: "
      f"{df.loc[W == 1, 'Y'].mean() - df.loc[W == 0, 'Y'].mean():+.3f} {UNIT}")

# %% [markdown]
# ## §2 Reproduce the published estimate for this outcome
# The paper fits a 2-period patient fixed-effects model (Stata xtreg, fe) with
# post x treatment, post x baseline value and measurement timing. With two
# periods this is algebraically a FIRST-DIFFERENCE regression,
#     change ~ treatment + baseline value + (time from consent, 6 mo - 0 mo),
# and the time-invariant covariates it lists (strata, age bins, sex,
# race/ethnicity) are absorbed by the fixed effects. We fit that regression
# ("paper model"; it reproduces the HbA1c estimate exactly) plus an ANCOVA
# that really adjusts for the prespecified covariates. The paper multiply-
# imputes missing baselines; we use observed baselines.
# Weight: Table 3's unadjusted difference (-8.14) is in POUNDS although labelled
# kg. Its adjusted +1.95 kg is in kg, but it depends on the paper's imputation of
# 5 missing baselines from demographics alone; with observed baselines it is
# about +0.8 kg (code/replicate_doyle.py ports the authors' code and shows this).

# %%
header("§2 Sanity check: reproduce the published estimate")
age_bins = " + ".join(f"age5_{i}" for i in range(1, 10))      # age5_10 = reference
adj = f"rural_site + hba1c_cat_control + {age_bins} + male + race_nonhisp_black + hispanic"
if OC["kind"] == "change":
    df["days_consent_to_6mo"] = df[f"{OC['date']}_6mo"] - df["consent_date"]
    df["d_time"] = df[f"{OC['date']}_6mo"] - df[f"{OC['date']}_0mo"]
    f_paper = f"Y ~ treat_now + {OC['val']}_0mo + d_time"
    f_ancova = (f"Y_level ~ treat_now + {OC['val']}_0mo + days_base_to_consent"
                f" + days_consent_to_6mo + {adj}")
else:
    df["Y_level_bin"] = df["Y"]                     # 0/100 indicator: same scale as the paper
    f_paper = f"Y_level_bin ~ treat_now + {OC['val']}_0mo"
    f_ancova = f"Y_level_bin ~ treat_now + {OC['val']}_0mo + {adj}"
paper_m = smf.ols(f_paper, data=df).fit(cov_type="HC1")
ancova = smf.ols(f_ancova, data=df).fit(cov_type="HC1")
b, (lo, hi) = ancova.params["treat_now"], ancova.conf_int().loc["treat_now"]
unadj = smf.ols(f"{f_ancova.split(' ~')[0]} ~ treat_now", data=df).fit(cov_type="HC1")
pa, pu = OC["paper_adj"], OC["paper_unadj"]
pm = (paper_m.params["treat_now"], *paper_m.conf_int().loc["treat_now"])
print(f"Unadjusted difference in 6-mo level (analytic sample): {fmt_ci(unadj.params['treat_now'], *unadj.conf_int().loc['treat_now'])}"
      f"   paper: {fmt_ci(*pu)}")
print(f"Paper model (first differences):     {fmt_ci(*pm)}   paper: {fmt_ci(*pa)}")
print(f"ANCOVA, prespecified covariates:     {fmt_ci(b, lo, hi)}")
if ARGS.outcome == "weight":                        # unit check against Table 3
    lb_u = smf.ols("weight_lb_val_6mo ~ treat_now", data=raw).fit(cov_type="HC1")   # all with 6-mo weight, as Table 3
    print(f"  Unadjusted in POUNDS: {fmt_ci(lb_u.params['treat_now'], *lb_u.conf_int().loc['treat_now'])}"
          " -> the published unadjusted weight difference is in lb, not kg.\n"
          "  The published adjusted +1.95 kg uses 5 imputed baselines; see outputs/replication/.")
RESULTS["replication"] = {"ancova_est": b, "ancova_lo": lo, "ancova_hi": hi,
                          "paper_model": list(pm),
                          "unadjusted_est": unadj.params["treat_now"],
                          "paper_adj": list(pa), "paper_unadj": list(pu)}

# %% [markdown]
# ## §3 Baseline balance (analytic sample)

# %%
header("§3 Baseline balance in the analytic sample (standardized mean differences)")
balance_vars = ["age", "male", "hispanic", "race_nonhisp_black", "rural_site",
                "hba1c_val_0mo", "days_a1c_to_consent", "bmi_val_0mo", "sys_val_0mo",
                "dias_val_0mo", "hdl_val_0mo", "ldl_val_0mo", "trigly_val_0mo",
                "drug_insulin_0mo", "drug_metformin_0mo", "drug_glp1_0mo",
                "op_n_0mo", "ed_or_ip_n_0mo"]
rows = []
for v in balance_vars:
    t, c = df.loc[W == 1, v].dropna(), df.loc[W == 0, v].dropna()
    sd = np.sqrt((t.var() + c.var()) / 2)
    rows.append({"variable": v, "mean_treat": t.mean(), "mean_control": c.mean(),
                 "smd": (t.mean() - c.mean()) / sd if sd > 0 else 0.0,
                 "pct_missing": df[v].isna().mean() * 100})
balance = pd.DataFrame(rows)
print(balance.round(3).to_string(index=False))
print("|SMD| > 0.1 flags chance imbalance worth adjusting for (all are in X).")
balance.to_csv(OUT / "table_balance.csv", index=False)

# %% [markdown]
# ## §4 Candidate effect modifiers (baseline only)
# Post-randomization variables (the survey indices, 6-mo medication use,
# engagement) are excluded: they are measured after assignment and are
# themselves affected by treatment. Baseline LDL (the single lipid marker,
# ~20% missing) enters X only with --ldl (patients with LDL measured; run id
# "<outcome>_ldl_cc"). Variable lists are in §1. Nothing is imputed.

# %%
header("§4 Candidate effect modifiers")


def build_X(data, use_lipids):
    """Baseline covariates as observed. The sample is complete case, so nothing is missing."""
    cols = CONTINUOUS + BINARY + (LIPIDS if use_lipids else [])
    Xr = data[cols].astype(float)
    assert not Xr.isna().any().any(), "missing covariate in the complete-case sample"
    return Xr


X_df = build_X(df, ARGS.ldl)
X = X_df.to_numpy()
FEATURES = list(X_df.columns)
print(f"{len(FEATURES)} features: {FEATURES}")
print(f"Events-per-variable analogue: N/p = {len(df) / len(FEATURES):.1f} patients per feature")

# %% [markdown]
# ## §5 Average treatment effect
# Propensity model: `StratumPropensity` returns the known randomization
# probability within each stratum (site x HbA1c category). Using an ML
# classifier here (as v1 did) only adds noise in an RCT.

# %%
header("§5 Average treatment effect")
STRATUM_COLS = [FEATURES.index("rural_site"), FEATURES.index("hba1c_cat_control")]
# Settings shared by every forest; passed to worker processes as JSON.
CFG = {"n_trees": N_TREES, "n_trees_nuis": N_TREES_NUIS, "k_folds": K_FOLDS,
       "min_leaf_cf": MIN_LEAF_CF, "stratum_cols": STRATUM_COLS}

cf_full = make_forest(SEED, CFG)
cf_full.fit(Y, W, X=X)
cf_ate = cf_full.ate_inference(X=X)
cf_lo, cf_hi = cf_ate.conf_int_mean(alpha=0.05)
print(f"Causal-forest ATE: {fmt_ci(cf_ate.mean_point, cf_lo, cf_hi)} {UNIT}")
RESULTS["ate_forest"] = {"est": float(cf_ate.mean_point), "lo": float(cf_lo), "hi": float(cf_hi)}
# (The AIPW ATE is reported in §6, from the cross-fitted doubly-robust scores.)

# %% [markdown]
# ## §6 Cross-fitted CATEs and doubly-robust scores
# For each of N_REPEATS random K-fold partitions (stratified on stratum x arm):
#   * fit a causal forest on K-1 folds and predict CATE tau_hat(x) for the
#     held-out fold -> every patient's CATE comes from a model that never saw
#     that patient's outcome;
#   * fit arm-specific outcome models mu_0, mu_1 on K-1 folds and form the
#     AIPW score  Gamma_i = mu1 - mu0 + W(Y - mu1)/e - (1-W)(Y - mu0)/(1-e)
#     with e = known stratum propensity. mean(Gamma) is the AIPW ATE and
#     Gamma_i is an unbiased (noisy) proxy for patient i's treatment effect.
# Results are aggregated across repeats by the median (Chernozhukov et al.).
# Each repeat runs in its own worker process (code/hte_lib.py) and is rerun
# if it crashes; see hte_lib's docstring for why.

# %%
header("§6 Cross-fitting (out-of-fold CATE and AIPW scores)")
MAX_WORKERS = max(1, min(N_REPEATS, (os.cpu_count() or 2) // 2))
MAX_ATTEMPTS = 4


def crossfit(Yv, Wv, Xv, e, strata, n_repeats, seed0, cfg=None):
    """Returns arrays (n_repeats x n) of out-of-fold CATE and AIPW scores."""
    cfg = cfg or CFG
    worker = str(ROOT / "code" / "hte_lib.py")
    with tempfile.TemporaryDirectory() as tmp:
        inp = Path(tmp) / "in.npz"
        np.savez(inp, Y=Yv, W=Wv, X=Xv, e=e, strata=strata)
        pending, running, attempts, out = list(range(n_repeats)), {}, {}, {}
        while pending or running:
            while pending and len(running) < MAX_WORKERS:
                r = pending.pop(0)
                attempts[r] = attempts.get(r, 0) + 1
                job = json.dumps({**cfg, "repeat": r, "seed0": seed0})
                running[r] = subprocess.Popen([sys.executable, worker, str(inp),
                                               str(Path(tmp) / f"out_{r}.npz"), job],
                                              stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            for r, p in list(running.items()):
                if p.poll() is None:
                    continue
                err = p.stderr.read().decode(errors="replace")
                del running[r]
                if p.returncode == 0:
                    out[r] = np.load(Path(tmp) / f"out_{r}.npz")
                    out[r] = (out[r]["tau"], out[r]["gamma"])
                    print(f"  repeat {r + 1}/{n_repeats} done"
                          + (f" (attempt {attempts[r]})" if attempts[r] > 1 else ""))
                elif attempts[r] < MAX_ATTEMPTS:
                    print(f"  repeat {r + 1}: worker exited with code {p.returncode}; rerunning")
                    pending.append(r)
                else:
                    raise RuntimeError(f"repeat {r} failed {MAX_ATTEMPTS} times:\n{err[-2000:]}")
            time.sleep(0.2)
    return (np.vstack([out[r][0] for r in range(n_repeats)]),
            np.vstack([out[r][1] for r in range(n_repeats)]))


STRATA = df["stratum"].to_numpy()
TAU, GAMMA = crossfit(Y, W, X, E_HAT, STRATA, N_REPEATS, SEED)
tau_bar = TAU.mean(axis=0)          # repeat-averaged out-of-fold CATE (descriptive)
gamma_bar = GAMMA.mean(axis=0)      # repeat-averaged AIPW score

aipw_est = np.median(GAMMA.mean(axis=1))
aipw_se = np.median(GAMMA.std(axis=1, ddof=1)) / np.sqrt(len(Y))
print(f"\nAIPW ATE (median over repeats): "
      f"{fmt_ci(aipw_est, aipw_est - 1.96 * aipw_se, aipw_est + 1.96 * aipw_se)} {UNIT}")
print(f"Out-of-fold CATE: mean {tau_bar.mean():+.3f}, SD {tau_bar.std():.3f}, "
      f"range [{tau_bar.min():+.2f}, {tau_bar.max():+.2f}]")
print(f"Share with predicted benefit (tau<0): {(tau_bar < 0).mean():.1%}  "
      "(descriptive only; not evidence of heterogeneity)")
RESULTS["ate_aipw"] = {"est": aipw_est, "se": aipw_se,
                       "lo": aipw_est - 1.96 * aipw_se, "hi": aipw_est + 1.96 * aipw_se}
RESULTS["cate_oof"] = {"mean": tau_bar.mean(), "sd": tau_bar.std(),
                       "min": tau_bar.min(), "max": tau_bar.max()}

# %% [markdown]
# ## §7 Omnibus tests: is there detectable heterogeneity at all?
# * **BLP calibration** (Chernozhukov, Demirer, Duflo & Fernández-Val 2018;
#   grf::test_calibration): regress Gamma on [1, tau_hat - mean(tau_hat)].
#   beta2 > 0 means the forest's out-of-fold ranking carries real signal;
#   beta2 ~ 1 means well calibrated. One-sided test H1: beta2 > 0.
# * **GATES**: mean Gamma within terciles of predicted CATE.
# * **RATE / AUTOC** (Yadlowsky et al. 2025, JASA): does prioritizing patients
#   by predicted benefit select patients who actually benefit more?
#   Benefit = -Gamma (HbA1c reduction); priority = -tau_hat.

# %%
header("§7 Omnibus heterogeneity tests (out-of-fold)")


def blp_test(gamma, tau, weights=None):
    Xb = sm.add_constant(tau - np.average(tau, weights=weights))
    m = (sm.WLS(gamma, Xb, weights=weights) if weights is not None else sm.OLS(gamma, Xb)).fit(cov_type="HC1")
    b2, se2 = m.params[1], m.bse[1]
    return {"beta1": m.params[0], "beta1_se": m.bse[0], "beta2": b2, "beta2_se": se2,
            "p_two_sided": m.pvalues[1], "p_one_sided": 1 - norm.cdf(b2 / se2)}


def gates(gamma, tau, n_groups=3):
    """Group ATEs by quantile of predicted CATE (G1 = most predicted benefit)."""
    g = pd.qcut(pd.Series(tau).rank(method="first"), n_groups, labels=False).to_numpy()
    out = []
    for k in range(n_groups):
        x = gamma[g == k]
        out.append({"group": k + 1, "n": len(x), "est": x.mean(), "se": x.std(ddof=1) / np.sqrt(len(x))})
    out = pd.DataFrame(out)
    diff = out.est.iloc[0] - out.est.iloc[-1]
    diff_se = np.sqrt(out.se.iloc[0] ** 2 + out.se.iloc[-1] ** 2)
    return out, diff, diff_se


def toc_curve(benefit, priority):
    """Targeting operator characteristic: TOC(q) = mean benefit in top-q - overall mean."""
    order = np.argsort(-priority, kind="stable")
    b = benefit[order]
    return np.cumsum(b) / np.arange(1, len(b) + 1) - b.mean()


def rate_autoc(benefit, priority, n_boot, seed):
    est = toc_curve(benefit, priority).mean()
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, len(benefit), len(benefit))
        boots[i] = toc_curve(benefit[idx], priority[idx]).mean()
    return est, boots.std(ddof=1)


blp_rows, gates_rows, rate_rows = [], [], []
for r in range(TAU.shape[0]):
    blp_rows.append(blp_test(GAMMA[r], TAU[r]))
    g, d, dse = gates(GAMMA[r], TAU[r])
    g["repeat"] = r; gates_rows.append(g)
    a, ase = rate_autoc(-GAMMA[r], -TAU[r], N_BOOT, SEED + r)
    rate_rows.append({"autoc": a, "se": ase, "gates_diff": d, "gates_diff_se": dse})
blp_df, rate_df = pd.DataFrame(blp_rows), pd.DataFrame(rate_rows)
gates_df = pd.concat(gates_rows)

BLP = blp_df.median()
print("BLP calibration test (median across repeats):")
print(f"  beta1 (mean forest prediction): {BLP.beta1:+.3f} (SE {BLP.beta1_se:.3f})")
print(f"  beta2 (differential prediction): {BLP.beta2:+.3f} (SE {BLP.beta2_se:.3f}); "
      f"one-sided p = {BLP.p_one_sided:.3f}; two-sided p = {BLP.p_two_sided:.3f}")
print(f"  per-repeat beta2: {np.round(blp_df.beta2.values, 2).tolist()}")
print("  Interpretation: beta2 <= 0 means the out-of-fold CATE ranking does NOT\n"
      "  track who actually responds (no detectable heterogeneity). A significantly\n"
      "  NEGATIVE beta2 indicates the forest is fitting noise, not signal.")

G = gates_df.groupby("group")[["n", "est", "se"]].median()
print("\nGATES by tercile of predicted CATE (median across repeats; G1 = most predicted benefit):")
for k, row in G.iterrows():
    print(f"  G{k}: n={int(row.n)}  {fmt_ci(row.est, row.est - 1.96 * row.se, row.est + 1.96 * row.se)}")
gd, gdse = rate_df.gates_diff.median(), rate_df.gates_diff_se.median()
print(f"  G1 - G3: {fmt_ci(gd, gd - 1.96 * gdse, gd + 1.96 * gdse)} "
      f"(negative = predicted responders really do better)")

AUTOC, AUTOC_SE = rate_df.autoc.median(), rate_df.se.median()
print(f"\nRATE (AUTOC), benefit = reduction in {OC['desc']}: {AUTOC:+.3f} (SE {AUTOC_SE:.3f}), "
      f"one-sided p (AUTOC>0) = {1 - norm.cdf(AUTOC / AUTOC_SE):.3f}")

RESULTS["blp"] = BLP.to_dict()
RESULTS["blp_beta2_by_repeat"] = blp_df.beta2.tolist()
RESULTS["gates"] = G.reset_index().to_dict(orient="records")
RESULTS["gates_g1_minus_g3"] = {"est": gd, "se": gdse}
RESULTS["rate_autoc"] = {"est": AUTOC, "se": AUTOC_SE, "p_one_sided": 1 - norm.cdf(AUTOC / AUTOC_SE)}
HETEROGENEITY_DETECTED = (BLP.p_one_sided < 0.05) and (RESULTS["rate_autoc"]["p_one_sided"] < 0.05)

# %% [markdown]
# ## §8 Which variables? (descriptive + linear effect-modification)
# * Forest split-importance (what v1 reported): describes the forest, is not a
#   test, and favors continuous variables with many split points.
# * BLP of the AIPW scores on standardized covariates (grf
#   `best_linear_projection`): each coefficient = change in treatment effect
#   per 1-SD (or 0->1) change in the covariate, adjusting for the others.
#   Holm and Benjamini-Hochberg adjustment across covariates.

# %%
header("§8 Effect modifiers: forest importance and BLP on covariates")
importance = pd.Series(cf_full.feature_importances_, index=FEATURES).sort_values(ascending=False)
print("Causal-forest split importance (top 10):")
print(importance.head(10).round(3).to_string())
importance.rename("importance").to_csv(OUT / "table_variable_importance.csv")

clinical = list(FEATURES)
Z = X_df[clinical].copy()
cont_in = [c for c in clinical if c in CONTINUOUS + LIPIDS]
Z[cont_in] = StandardScaler().fit_transform(Z[cont_in])
# As in docs/02: baseline weight (~BMI) is left out of the joint linear projection only; the
# ANCOVA interaction test below uses every modifier (Z_all). Until 2026-11-09, drug_any_db_0mo
# and hba1c_cat_control were also dropped from both as "near-collinear"; they are not (VIF 1.4,
# 2.2; full rank; verification/multiagent_2026-10-29/C_compare), so the PI chose to follow the plan.
Z_all = Z.copy()
Z = Z.drop(columns=[c for c in ["weight_val_0mo"] if c in Z])
blp_cov = sm.OLS(gamma_bar, sm.add_constant(Z)).fit(cov_type="HC1")
bc = pd.DataFrame({"coef": blp_cov.params, "se": blp_cov.bse, "p": blp_cov.pvalues}).drop("const")
bc["p_holm"] = multipletests(bc.p, method="holm")[1]
bc["p_bh"] = multipletests(bc.p, method="fdr_bh")[1]
bc["label"] = [LABELS.get(i, i) for i in bc.index]
bc = bc.sort_values("p")
print("\nBLP of AIPW scores on covariates (continuous per 1 SD; negative = larger benefit):")
print(bc[["label", "coef", "se", "p", "p_holm", "p_bh"]].round(3).to_string())
wald_p = blp_cov.f_test(np.eye(len(blp_cov.params))[1:]).pvalue
print(f"Joint Wald test that all covariate coefficients = 0: p = {float(wald_p):.3f}")
bc.to_csv(OUT / "table_blp_covariates.csv")
RESULTS["blp_covariates_joint_p"] = float(wald_p)
RESULTS["blp_covariates_min_p_holm"] = float(bc.p_holm.min())

# Classical parametric cross-check: ANCOVA with all treatment x covariate interactions.
Z = Z_all
inter = smf.ols("Y ~ treat_now * (" + " + ".join(Z.columns) + ")",
                data=pd.concat([df[["Y", "treat_now"]], Z], axis=1)).fit(cov_type="HC1")
inter_terms = [t for t in inter.params.index if t.startswith("treat_now:")]
R_mat = np.array([[1.0 if p == t else 0.0 for p in inter.params.index] for t in inter_terms])
inter_p = float(inter.f_test(R_mat).pvalue)
print(f"ANCOVA joint test of all treatment x covariate interactions: p = {inter_p:.3f}")
RESULTS["ancova_interactions_joint_p"] = inter_p

# Randomization check of the same interaction test (Ding, Feller & Miratrix
# 2016, JRSS-B). The HC1 Wald test with ~20 interactions in a few hundred
# patients can over-reject, especially with skewed or binary outcomes. Under
# H0 "constant additive effect tau", Y(0) = Y - tau*W; re-randomizing W within
# strata gives the exact null distribution of the Wald statistic. tau is
# plugged in at the ANCOVA estimate (DFM maximize over a CI; with a null
# average effect the plug-in is the usual approximation).
N_PERM = 200 if FAST else 2000
Xmain = sm.add_constant(Z.to_numpy(float))
tau0 = float(smf.ols("Y ~ treat_now + " + " + ".join(Z.columns),
                     data=pd.concat([df[["Y", "treat_now"]], Z], axis=1)).fit().params["treat_now"])
Y0 = Y - tau0 * W


def inter_wald(w, y):
    D = np.column_stack([Xmain, w, Xmain[:, 1:] * w[:, None]])
    m = sm.OLS(y, D).fit(cov_type="HC1")
    k = Xmain.shape[1] + 1                          # first interaction column
    return float(m.wald_test(np.eye(D.shape[1])[k:], scalar=True).statistic)


rng_p = np.random.default_rng(SEED + 99)
stat_obs = inter_wald(W.astype(float), Y)
perm = np.empty(N_PERM)
for i in range(N_PERM):
    w_p = W.copy()
    for s in np.unique(STRATA):                     # re-randomize within strata
        idx = np.where(STRATA == s)[0]
        w_p[idx] = rng_p.permutation(W[idx])
    perm[i] = inter_wald(w_p.astype(float), Y0 + tau0 * w_p)
inter_perm_p = (1 + np.sum(perm >= stat_obs)) / (N_PERM + 1)
print(f"  randomization p for the same statistic ({N_PERM} stratified re-randomizations, "
      f"constant effect {tau0:+.2f}): p = {inter_perm_p:.3f}")
RESULTS["ancova_interactions_randomization_p"] = float(inter_perm_p)

# %% [markdown]
# ## §9 Pre-specified subgroups
# Family A = Doyle's PRE-SPECIFIED subgroups (protocol SAP 9.1.1: baseline HbA1c
# and site; Supplement 2 eTables 17-19: also sex, household size, ZIP-level
# income/distance/education, prior-year inpatient and ED use, diabetes
# knowledge, insulin and metformin, predicted engagement) that the public data
# can reproduce. Household size, ZIP variables, dispensed insulin and predicted
# engagement are not in the public file; diabetes knowledge is only measured
# at follow-up. Family B = subgroups we ADDED (not in Doyle): age, race/
# ethnicity, staleness of the baseline A1c. Holm adjustment within each family.
# Effects are AIPW (mean Gamma) within subgroup; the interaction p compares the
# two levels. `paper_model_est` = the paper's model fitted within the subgroup,
# for comparison with eTables 18-19 (which use a somewhat larger sample).

# %%
header("§9 Pre-specified subgroups (AIPW within subgroup)")
df["age_ge65"] = (df["age"] >= 65).astype(int)
df["nonhisp_white"] = ((df["hispanic"] == 0) & (df["race"] == "White")).astype(int)
df["a1c_test_gt90d"] = (df["days_a1c_to_consent"] > 90).astype(int)
PRESPEC = [  # (variable, name, label level 0, label level 1, family)
    ("hba1c_cat_control", "Baseline HbA1c", "<9.5%", ">=9.5%", "Doyle pre-specified"),
    ("rural_site", "Site", "Urban", "Rural", "Doyle pre-specified"),
    ("male", "Sex", "Female", "Male", "Doyle pre-specified"),
    ("ip_i_0mo", "Prior-year inpatient admission", "No", "Yes", "Doyle pre-specified"),
    ("ed_i_0mo", "Prior-year ED visit", "No", "Yes", "Doyle pre-specified"),
    ("drug_insulin_0mo", "Insulin prescribed", "No", "Yes", "Doyle pre-specified"),
    ("drug_metformin_0mo", "Metformin prescribed", "No", "Yes", "Doyle pre-specified"),
    ("age_ge65", "Age", "<65 y", ">=65 y", "Added (not in Doyle)"),
    ("nonhisp_white", "Race/ethnicity", "Hispanic or non-White", "Non-Hispanic White", "Added (not in Doyle)"),
    ("a1c_test_gt90d", "Baseline A1c test before consent", "<=90 d", ">90 d", "Added (not in Doyle)"),
]


def subgroup_effect(mask, gamma=gamma_bar):
    g = gamma[mask]
    se = g.std(ddof=1) / np.sqrt(len(g))
    dm = Y[mask & (W == 1)].mean() - Y[mask & (W == 0)].mean()
    return g.mean(), se, dm


pre_rows = []
for var, name, lab0, lab1, fam in PRESPEC:
    v = df[var].to_numpy()
    res = {}
    for level, lab in [(0, lab0), (1, lab1)]:
        m = v == level
        est, se, dm = subgroup_effect(m)
        res[level] = (est, se)
        pm_sub = smf.ols(f_paper, data=df[m]).fit(cov_type="HC1").params["treat_now"]
        pre_rows.append({"family": fam, "subgroup": name, "level": lab,
                         "n": int(m.sum()), "n_treat": int((m & (W == 1)).sum()),
                         "est": est, "lo": est - 1.96 * se, "hi": est + 1.96 * se,
                         "diff_in_means": dm, "paper_model_est": pm_sub})
    d = res[1][0] - res[0][0]; dse = np.hypot(res[1][1], res[0][1])
    p = 2 * (1 - norm.cdf(abs(d / dse)))
    pre_rows[-1]["p_interaction"] = p; pre_rows[-2]["p_interaction"] = p
pre = pd.DataFrame(pre_rows)
p_unique = pre.drop_duplicates("subgroup")[["family", "subgroup", "p_interaction"]].copy()
p_unique["p_holm"] = p_unique.groupby("family")["p_interaction"].transform(
    lambda x: multipletests(x, method="holm")[1])
pre = pre.merge(p_unique[["subgroup", "p_holm"]], on="subgroup")
print(pre[["family", "subgroup", "level", "n", "n_treat", "est", "lo", "hi", "paper_model_est",
           "p_interaction", "p_holm"]].round(3).to_string(index=False))
PRE_DOYLE = pre[pre.family == "Doyle pre-specified"]

# %% [markdown]
# ## §10 ACCORD-style data-driven subgroups
# (a) **Tertiles** of the top continuous modifiers by forest importance, as in
#     Kianmehr 2022. Cut-points do not use outcomes, but the variables were
#     chosen by the forest on the same data -> p-values are optimistic.
# (b) **Honest causal tree** (Athey & Imbens 2016 logic; Jiao 2024 Fig. 1-3):
#     a 50/50 split stratified on stratum x arm. The DISCOVERY half (with its
#     own cross-fitted AIPW scores) grows a depth-2 tree; the ESTIMATION half
#     (its own AIPW scores) estimates leaf effects. Single-variable splits of
#     the top N_TOP_SPLITS (=5) modifiers by forest importance are found and
#     estimated the same way (Jiao 2024: one tree per top-5 variable).
# (c) **Split stability**: repeat the discovery step N_STABILITY times and
#     record which variable is chosen at the root (Kianmehr 2022 repeated the
#     honest tree 4000 times and kept the most frequently selected variable).
#     A real modifier should be chosen consistently; noise scatters.

# %%
header("§10a Tertiles of the top continuous modifiers (exploratory)")
top_cont = [v for v in importance.index if v in CONTINUOUS + LIPIDS][:3]
tert_rows = []
for v in top_cont:
    cuts = df[v].quantile([1 / 3, 2 / 3]).to_numpy()
    t = np.digitize(df[v], cuts, right=True)
    for k in range(3):
        m = t == k
        est, se, dm = subgroup_effect(m)
        x = df.loc[m, v]
        tert_rows.append({"family": "tertile", "subgroup": LABELS[v], "level": f"T{k + 1}: {x.min():g}-{x.max():g}",
                          "n": int(m.sum()), "n_treat": int((m & (W == 1)).sum()), "est": est,
                          "lo": est - 1.96 * se, "hi": est + 1.96 * se, "diff_in_means": dm})
    trend = sm.OLS(gamma_bar, sm.add_constant(t.astype(float))).fit(cov_type="HC1")
    for row in tert_rows[-3:]:
        row["p_interaction"] = trend.pvalues[1]
    print(f"{LABELS[v]}: linear trend across tertiles p = {trend.pvalues[1]:.3f}")
tert = pd.DataFrame(tert_rows)
print(tert[["subgroup", "level", "n", "est", "lo", "hi"]].round(3).to_string(index=False))

header("§10b Honest (sample-split) causal tree")
disc, est_idx = train_test_split(np.arange(len(Y)), test_size=0.5, random_state=SEED,
                                 stratify=STRATA * 2 + W)


def within_half_gamma(idx, seed):
    """Cross-fitted AIPW scores using only the patients in idx."""
    g = np.zeros(len(idx))
    skf = StratifiedKFold(n_splits=K_FOLDS, shuffle=True, random_state=seed)
    for tr, te in skf.split(idx, (STRATA * 2 + W)[idx]):
        g[te] = aipw_scores(Y, W, X, E_HAT, idx[tr], idx[te], seed, N_TREES_NUIS)
    return g


g_disc, g_est = within_half_gamma(disc, SEED + 1), within_half_gamma(est_idx, SEED + 2)
cand = [c for c in clinical if c not in ("drug_any_db_0mo",)]
Xc = X_df[cand]
tree = DecisionTreeRegressor(max_depth=2, min_samples_leaf=MIN_LEAF_TREE,
                             random_state=SEED).fit(Xc.iloc[disc], g_disc)


def leaf_rules(tree, names):
    """Human-readable rule for each leaf id."""
    t, rules = tree.tree_, {}
    def walk(node, path):
        if t.feature[node] == -2:
            rules[node] = " & ".join(path) or "all"
            return
        f, thr = names[t.feature[node]], t.threshold[node]
        walk(t.children_left[node], path + [f"{LABELS.get(f, f)} <= {thr:.4g}"])
        walk(t.children_right[node], path + [f"{LABELS.get(f, f)} > {thr:.4g}"])
    walk(0, [])
    return rules


rules = leaf_rules(tree, cand)
leaf_est = tree.apply(Xc.iloc[est_idx])
leaf_disc = tree.apply(Xc.iloc[disc])
tree_rows = []
for leaf, rule in rules.items():
    m = leaf_est == leaf
    if m.sum() < 5:
        continue
    g = g_est[m]; se = g.std(ddof=1) / np.sqrt(m.sum())
    tree_rows.append({"family": "honest_tree", "subgroup": "Depth-2 tree leaf", "level": rule,
                      "n": int(m.sum()), "n_treat": int(W[est_idx][m].sum()),
                      "est_discovery": g_disc[leaf_disc == leaf].mean(),
                      "est": g.mean(), "lo": g.mean() - 1.96 * se, "hi": g.mean() + 1.96 * se})
tree_tab = pd.DataFrame(tree_rows)
leaf_dum = pd.get_dummies(leaf_est, drop_first=True).astype(float)
tree_p = float(sm.OLS(g_est, sm.add_constant(leaf_dum)).fit(cov_type="HC1").f_pvalue) if leaf_dum.shape[1] else np.nan
tree_tab["p_interaction"] = tree_p
print("Leaves: discovery-half estimate (optimistic) vs estimation-half estimate (honest):")
print(tree_tab[["level", "n", "est_discovery", "est", "lo", "hi"]].round(3).to_string(index=False))
print(f"Honest test of any difference between leaves (estimation half): p = {tree_p:.3f}")
print("Shrinkage from discovery to estimation half shows how much of the 'signal' was noise.")

# Single-variable ACCORD-style splits for the top N_TOP_SPLITS modifiers (honest).
split_rows = []
for v in [x for x in importance.index if x in cand][:N_TOP_SPLITS]:
    st = DecisionTreeRegressor(max_depth=1, min_samples_leaf=MIN_LEAF_TREE,
                               random_state=SEED).fit(Xc[[v]].iloc[disc], g_disc)
    if st.tree_.feature[0] == -2:
        continue
    thr = st.tree_.threshold[0]
    xv = Xc[v].to_numpy()[est_idx]
    ests = {}
    for side, m in [("<=", xv <= thr), (">", xv > thr)]:
        g = g_est[m]; se = g.std(ddof=1) / np.sqrt(m.sum())
        ests[side] = (g.mean(), se)
        split_rows.append({"family": "honest_split", "subgroup": LABELS.get(v, v), "level": f"{side} {thr:.4g}",
                           "n": int(m.sum()), "n_treat": int(W[est_idx][m].sum()),
                           "est": g.mean(), "lo": g.mean() - 1.96 * se, "hi": g.mean() + 1.96 * se})
    d = ests["<="][0] - ests[">"][0]; dse = np.hypot(ests["<="][1], ests[">"][1])
    p = 2 * (1 - norm.cdf(abs(d / dse)))
    split_rows[-1]["p_interaction"] = p; split_rows[-2]["p_interaction"] = p
splits = pd.DataFrame(split_rows, columns=["family", "subgroup", "level", "n", "n_treat",
                                           "est", "lo", "hi", "p_interaction"])
print("\nHonest single-variable splits (threshold from discovery half, estimate from estimation half):")
print(splits[["subgroup", "level", "n", "est", "lo", "hi", "p_interaction"]].round(3).to_string(index=False)
      if len(splits) else "  no admissible split")

header("§10c Split stability of the root variable")
root_vars, root_thr = [], []
for s in range(N_STABILITY):
    d_idx, _ = train_test_split(np.arange(len(Y)), test_size=0.5, random_state=SEED + 1000 + s,
                                stratify=STRATA * 2 + W)
    t1 = DecisionTreeRegressor(max_depth=1, min_samples_leaf=MIN_LEAF_TREE,
                               random_state=s).fit(Xc.iloc[d_idx], gamma_bar[d_idx])
    if t1.tree_.feature[0] != -2:
        root_vars.append(cand[t1.tree_.feature[0]]); root_thr.append(t1.tree_.threshold[0])
stab = pd.Series(root_vars).map(lambda x: LABELS.get(x, x)).value_counts(normalize=True)
print(f"Root-split variable across {N_STABILITY} random half-samples (share of splits):")
print(stab.head(8).round(3).to_string())
stab.rename("share_of_root_splits").to_csv(OUT / "table_split_stability.csv")
RESULTS["split_stability_top"] = stab.head(5).to_dict()

# %% [markdown]
# ## §11 Sensitivity analyses
# S1 Change outcomes: Y = 6-mo LEVEL (paper's scale). Because the baseline
#    value is in X, CATE(level) = CATE(change) in the population; only the
#    nuisance-model efficiency differs. A different answer = instability.
#    ED/IP: Y = number of ED/IP visits (count) instead of any use.
# S2 LDL runs only: baseline LDL removed from X (same patients). In the main runs
#    adding LDL would need imputation; the "<outcome>_ldl_cc" run covers it.
# S3 Attrition / exclusion: inverse-probability-of-inclusion weights from a
#    logistic model in all 465 randomized-and-retained patients, using the
#    covariates observed for all of them. Skipped when nobody is excluded.
# S4 Restrict to patients whose baseline value of the outcome was measured
#    <=180 days before consent (not applicable to ED/IP).
# PC Positive control: inject known heterogeneity (0.25, 0.5, 1 SD of Y) and
#    check it is detected.

# %%
header("§11 Sensitivity analyses")
sens = []


def omnibus(Yv, Xv, mask=None, label="", n_rep=N_REPEATS_SENS, weights=None):
    m = np.ones(len(Yv), bool) if mask is None else mask
    T, Gm = crossfit(Yv[m], W[m], Xv[m], E_HAT[m], STRATA[m], n_rep, SEED + 7)
    b = pd.DataFrame([blp_test(Gm[r], T[r], None if weights is None else weights[m])
                      for r in range(n_rep)]).median()
    a = pd.DataFrame([rate_autoc(-Gm[r], -T[r], N_BOOT, SEED + r) for r in range(n_rep)],
                     columns=["autoc", "se"]).median()
    ate = np.median([np.average(Gm[r], weights=None if weights is None else weights[m]) for r in range(n_rep)])
    row = {"analysis": label, "n": int(m.sum()), "ate_aipw": ate, "blp_beta2": b.beta2,
           "blp_beta2_se": b.beta2_se, "blp_p_one_sided": b.p_one_sided,
           "autoc": a.autoc, "autoc_se": a.se, "autoc_p_one_sided": 1 - norm.cdf(a.autoc / a.se)}
    print(f"  {label}: N={row['n']}, ATE {ate:+.3f}, BLP beta2 {b.beta2:+.2f} "
          f"(one-sided p {b.p_one_sided:.3f}), AUTOC {a.autoc:+.3f} (p {row['autoc_p_one_sided']:.3f})")
    return row


sens.append({"analysis": "Main", "n": len(Y), "ate_aipw": aipw_est, "blp_beta2": BLP.beta2,
             "blp_beta2_se": BLP.beta2_se, "blp_p_one_sided": BLP.p_one_sided,
             "autoc": AUTOC, "autoc_se": AUTOC_SE,
             "autoc_p_one_sided": RESULTS["rate_autoc"]["p_one_sided"]})
if OC["kind"] == "change":
    print("S1 outcome = 6-mo level")
    sens.append(omnibus(df["Y_level"].to_numpy(), X, label="S1 Outcome = 6-mo level"))
else:
    print("S1 outcome = number of ED/IP visits in 6 months")
    sens.append(omnibus(df["Y_level"].to_numpy(float), X, label="S1 Outcome = no. of ED/IP visits"))

if ARGS.ldl:
    print("S2 LDL removed (same patients)")
    sens.append(omnibus(Y, build_X(df, False).to_numpy(), label="S2 LDL removed"))
else:
    print("S2 not run: adding LDL needs patients with LDL measured (see the _ldl_cc run)")

if raw["in_sample"].all():
    print("S3 skipped: every randomized-and-retained patient is in the analytic sample")
    RESULTS["attrition_model_auc"] = None
else:
    print("S3 inverse probability of inclusion weights (missing outcome / incomplete covariates)")
    # Inclusion model on the covariates observed for all 465 (no imputation). The
    # ones that define exclusion (BMI, BP, LDL, baseline weight) are partly missing
    # and so cannot enter; exclusion on them is assumed ignorable given the rest.
    always = [c for c in CONTINUOUS + BINARY if raw[c].notna().all()]
    R_all = raw[always].astype(float).copy()
    R_all["treat_now"] = raw["treat_now"]
    R_s = StandardScaler().fit_transform(R_all)
    obs_model = LogisticRegression(C=1.0, max_iter=2000).fit(R_s, raw["in_sample"])
    p_obs = obs_model.predict_proba(R_s)[:, 1]
    auc = roc_auc_score(raw["in_sample"], p_obs)
    ipw = (1 / p_obs[raw["in_sample"].to_numpy() == 1])
    ipw = ipw / ipw.mean()
    print(f"  inclusion model AUC = {auc:.2f} (0.5 = exclusion unrelated to baseline X); "
          f"weights range {ipw.min():.2f}-{ipw.max():.2f}")
    sens.append(omnibus(Y, X, label="S3 IPW for excluded patients", weights=ipw))
    RESULTS["attrition_model_auc"] = auc

if OC["date"] is not None:
    print("S4 baseline value measured <=180 days before consent")
    recent = df["days_base_to_consent"].to_numpy() <= 180
    sens.append(omnibus(Y, X, mask=recent, label="S4 Baseline measured <=180 d before consent"))

# Positive control (power check): add a KNOWN heterogeneous effect to the real
# outcomes -- treated patients aged >= median get an extra -delta/2, younger
# treated patients +delta/2 (mean effect unchanged). delta is set in SDs of Y
# so that outcomes are comparable. If the pipeline cannot detect a contrast of
# size delta between two halves of the sample, a null result above says little
# about smaller real heterogeneity.
print("Positive control: injected age-based heterogeneity of known size")
older = (df["age"] >= df["age"].median()).to_numpy()
SD_Y = Y.std(ddof=1)
for k in PC_SD:
    delta = k * SD_Y
    Y_pc = Y + W * np.where(older, -delta / 2, delta / 2)
    sens.append(omnibus(Y_pc, X, label=f"PC injected contrast {delta:.2f} {UNIT} = {k:g} SD (older vs younger)"))
RESULTS["sd_y"] = SD_Y
sens_df = pd.DataFrame(sens)
sens_df.to_csv(OUT / "table_sensitivity.csv", index=False)

# %% [markdown]
# ## §12 Figures

# %%
header("§12 Figures")

# Fig 1: forest split importance
fig, ax = plt.subplots(figsize=(7, 5))
imp = importance.head(12)[::-1]
ax.barh([LABELS.get(i, i) for i in imp.index], imp.values, color=ACCENT, height=0.6)
ax.set_xlabel("Share of causal-forest splits (descriptive, not a test)")
ax.set_title("Fig 1. Causal-forest variable importance", loc="left")
ax.grid(axis="y", visible=False)
fig.tight_layout(); fig.savefig(OUT / "fig1_variable_importance.png", dpi=200); plt.close(fig)

# Fig 2: out-of-fold CATE distribution
fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(tau_bar, bins=30, color=ACCENT, edgecolor="#fcfcfb", linewidth=1)
ax.axvline(0, color=INK2, ls="--", lw=1)
ax.axvline(aipw_est, color=ACCENT2, lw=2, label=f"AIPW ATE = {aipw_est:+.2f}")
ax.set_xlabel(f"Out-of-fold predicted CATE, {OC['desc']} at 6 mo ({UNIT}; negative = lower with program)")
ax.set_ylabel("Patients")
ax.set_title("Fig 2. Spread of predicted effects (spread alone does not imply heterogeneity)",
             loc="left", fontsize=10)
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(OUT / "fig2_cate_distribution.png", dpi=200); plt.close(fig)

# Fig 3: GATES
fig, ax = plt.subplots(figsize=(6, 4))
ax.errorbar(G.index, G.est, yerr=1.96 * G.se, fmt="o", color=ACCENT, ms=8, capsize=4, lw=2)
ax.axhline(aipw_est, color=ACCENT2, lw=1.5, label="Overall AIPW ATE")
ax.axhline(0, color=INK2, ls="--", lw=1)
ax.set_xticks([1, 2, 3], ["G1\n(most predicted\nbenefit)", "G2", "G3\n(least predicted\nbenefit)"])
ax.set_ylabel(f"Observed effect (AIPW), {UNIT}")
ax.set_title("Fig 3. Group ATEs by tercile of out-of-fold CATE", loc="left")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(OUT / "fig3_gates.png", dpi=200); plt.close(fig)

# Fig 4: TOC curve (repeat 0) with bootstrap band
benefit, prio = -GAMMA[0], -TAU[0]
qs = np.arange(1, len(Y) + 1) / len(Y)
toc = toc_curve(benefit, prio)
rng = np.random.default_rng(SEED)
boot_toc = np.array([toc_curve(benefit[i], prio[i]) for i in
                     (rng.integers(0, len(Y), len(Y)) for _ in range(min(N_BOOT, 500)))])
lo_b, hi_b = np.percentile(boot_toc, [2.5, 97.5], axis=0)
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.fill_between(qs, lo_b, hi_b, color=ACCENT, alpha=0.18, lw=0)
ax.plot(qs, toc, color=ACCENT, lw=2)
ax.axhline(0, color=INK2, ls="--", lw=1)
ax.set_xlabel("Fraction treated, prioritized by predicted benefit")
ax.set_ylabel(f"Extra reduction vs. treating at random ({UNIT})")
ax.set_title(f"Fig 4. Targeting curve (TOC); AUTOC = {AUTOC:+.2f} (SE {AUTOC_SE:.2f})", loc="left")
ax.set_xlim(0, 1)
fig.tight_layout(); fig.savefig(OUT / "fig4_toc_curve.png", dpi=200); plt.close(fig)

# Fig 5: subgroup forest plot (pre-specified, tertiles, honest splits)
plot_df = pd.concat([pre.assign(section=pre.family),
                     tert.assign(section="Tertiles (exploratory)"),
                     splits.assign(section="Honest single splits")], ignore_index=True)
fig, ax = plt.subplots(figsize=(8.5, 0.32 * len(plot_df) + 1.6))
ypos, ylabels, heads, y = [], [], [], 0
for sec, chunk in plot_df.groupby("section", sort=False):
    heads.append(len(ylabels)); ylabels.append(sec); ypos.append(y); y += 1
    for _, r in chunk.iterrows():
        ax.errorbar(r.est, y, xerr=[[r.est - r.lo], [r.hi - r.est]], fmt="o", color=ACCENT,
                    ms=5, capsize=3, lw=1.5)
        ylabels.append(f"{r.subgroup}: {r.level} (n={r.n})"); ypos.append(y); y += 1
    y += 0.5
ax.axvline(0, color=INK2, ls="--", lw=1)
ax.axvline(aipw_est, color=ACCENT2, lw=1.2)
ax.set_yticks(ypos, ylabels, fontsize=8)
for i in heads:                          # section headers in bold ink
    ax.get_yticklabels()[i].set_fontweight("bold"); ax.get_yticklabels()[i].set_color(INK)
ax.invert_yaxis(); ax.grid(axis="y", visible=False)
ax.set_xlabel(f"{EFFECT_LAB} (AIPW)\n{OC['better']}")
ax.set_title("Fig 5. Subgroup effects (orange line = overall ATE)", loc="left")
fig.tight_layout(); fig.savefig(OUT / "fig5_subgroups.png", dpi=200); plt.close(fig)

# Fig 6: honest tree leaves, discovery vs estimation half
if len(tree_tab):
    fig, ax = plt.subplots(figsize=(8, 0.6 * len(tree_tab) + 1.5))
    yy = np.arange(len(tree_tab))
    ax.scatter(tree_tab.est_discovery, yy - 0.15, color=INK2, marker="s", s=40,
               label="Discovery half (used to grow the tree)")
    ax.errorbar(tree_tab.est, yy + 0.15, xerr=[tree_tab.est - tree_tab.lo, tree_tab.hi - tree_tab.est],
                fmt="o", color=ACCENT, ms=7, capsize=3, label="Estimation half (honest, 95% CI)")
    ax.axvline(0, color=INK2, ls="--", lw=1)
    ax.set_yticks(yy, [f"{r.level}\n(n={r.n})" for r in tree_tab.itertuples()], fontsize=8)
    ax.invert_yaxis(); ax.grid(axis="y", visible=False)
    ax.set_xlabel(EFFECT_LAB)
    ax.set_title(f"Fig 6. Honest depth-2 causal tree (leaf-difference p = {tree_p:.2f})", loc="left")
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=2)
    fig.tight_layout(); fig.savefig(OUT / "fig6_honest_tree.png", dpi=200); plt.close(fig)
print(f"Figures written to {OUT.relative_to(ROOT)}/")

# %% [markdown]
# ## §13 Write results

# %%
header("§13 Results files")
subgroups_all = pd.concat([pre, tert, splits, tree_tab], ignore_index=True)
subgroups_all.to_csv(OUT / "table_subgroups.csv", index=False)
pd.DataFrame({"tau_oof": tau_bar, "gamma_aipw": gamma_bar, "treat_now": W, "Y": Y,
              "subject_id": df["subject_id"]}).to_csv(OUT / "patient_level_cate.csv", index=False)
RESULTS["heterogeneity_detected"] = bool(HETEROGENEITY_DETECTED)
RESULTS["prespecified_min_p_holm"] = float(PRE_DOYLE.p_holm.min())
RESULTS["added_subgroups_min_p_holm"] = float(pre[pre.family != "Doyle pre-specified"].p_holm.min())
RESULTS["honest_tree_leaf_p"] = tree_p
RESULTS["honest_splits_min_p"] = float(splits.p_interaction.min()) if len(splits) else None
RESULTS["n"], RESULTS["n_treat"] = int(len(Y)), int(W.sum())
RESULTS["sensitivity"] = sens_df.to_dict(orient="records")


def to_builtin(o):
    if isinstance(o, dict):
        return {str(k): to_builtin(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [to_builtin(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


(OUT / "results.json").write_text(json.dumps(to_builtin(RESULTS), indent=2))

verdict = ("Formal tests DETECT heterogeneity (BLP and RATE both p<0.05); subgroup findings "
           "below are candidates for replication." if HETEROGENEITY_DETECTED else
           "Formal out-of-fold tests do NOT detect treatment-effect heterogeneity. Subgroup "
           "estimates below should be read as consistent with a uniform (null) effect; any "
           "nominally significant subgroup is expected under multiple testing.")
md = f"""# HTE results: food-as-medicine trial, {OC['desc']} at 6 months (auto-generated)

Run `{RUN_ID}` ({'FAST smoke test' if FAST else 'full'}): LDL in main model: {ARGS.ldl};
complete case (no imputation); seed {SEED}; {N_TREES} trees; {K_FOLDS}-fold cross-fitting
x {N_REPEATS} repeats. Sign: {OC['better']}.

## Bottom line
{verdict}

## Average effect (N={len(Y)})
| Estimator | {EFFECT_LAB} (95% CI) |
|---|---|
| Paper, adjusted (longitudinal FE; N={OC['paper_n']}{'; ' + OC['paper_note'] if 'paper_note' in OC else ''}) | {fmt_ci(*OC['paper_adj'])} |
| Paper model replicated (first differences) | {fmt_ci(*pm)} |
| ANCOVA, prespecified covariates | {fmt_ci(b, lo, hi)} |
| AIPW (cross-fitted, known propensity) | {fmt_ci(aipw_est, aipw_est - 1.96 * aipw_se, aipw_est + 1.96 * aipw_se)} |
| Causal forest | {fmt_ci(cf_ate.mean_point, cf_lo, cf_hi)} |

## Omnibus heterogeneity tests (out-of-fold)
| Test | Estimate | p |
|---|---|---|
| BLP calibration beta2 (H1: >0) | {BLP.beta2:+.2f} (SE {BLP.beta2_se:.2f}) | {BLP.p_one_sided:.3f} (one-sided) |
| GATES G1 - G3 | {gd:+.2f} (SE {gdse:.2f}) | {2 * (1 - norm.cdf(abs(gd / gdse))):.3f} |
| RATE AUTOC (H1: >0) | {AUTOC:+.3f} (SE {AUTOC_SE:.3f}) | {RESULTS['rate_autoc']['p_one_sided']:.3f} (one-sided) |
| Joint BLP on covariates | - | {float(wald_p):.3f} |
| ANCOVA joint interactions (HC1 Wald) | - | {inter_p:.3f} |
| Same, randomization p (constant-effect null) | - | {inter_perm_p:.3f} |

## Pre-specified subgroups (Doyle family: smallest Holm-adjusted interaction p = {PRE_DOYLE.p_holm.min():.3f})
{pre[['family', 'subgroup', 'level', 'n', 'est', 'lo', 'hi', 'paper_model_est', 'p_interaction', 'p_holm']].round(3).to_markdown(index=False)}

## Honest depth-2 tree (estimation half; leaf-difference p = {tree_p:.3f})
{tree_tab[['level', 'n', 'est_discovery', 'est', 'lo', 'hi']].round(3).to_markdown(index=False) if len(tree_tab) else 'No split found.'}

## Sensitivity analyses
{sens_df.round(3).to_markdown(index=False)}

## Root-split stability ({N_STABILITY} half-samples)
{stab.head(6).round(3).to_markdown()}

Figures: fig1-fig6 in this folder. Full console log: console_log.txt.
"""
(OUT / "results_summary.md").write_text(md)
print(f"Wrote results_summary.md, results.json, table_*.csv, patient_level_cate.csv to {OUT.relative_to(ROOT)}/")
print("\n" + verdict)
