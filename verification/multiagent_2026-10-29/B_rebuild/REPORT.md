# Agent B: blind independent rebuild (29 Oct 2026) — 2 of 4 outcomes

Saved by the main session from Agent B's final message (the agent's tool was blocked from writing
this file). Content is the agent's, lightly condensed.

**Blinding:** complied. Used only `docs/02_analysis_plan.md`, `data/raw/` (CSV, README.pdf data
dictionary, `01_published_tables.do`, SOURCE.md) and `requirements.txt`. Did not open `code/`,
`outputs/`, `notebooks/`, `manuscript/`, `docs/01`, `docs/03-08`, `WORKFLOW.md`, the project README
or other `verification/` subfolders (saw only the folder names A_rerun and E_red_team in a
listing).

**Status:** HbA1c and weight finished. SBP, ED/IP and the positive control not run: the run script
was stopped by Claude Code for low system memory (~2.2 of 15.7 GB free) and not restarted. To
finish: `sh run_all.sh`, or `python -W ignore b_analysis.py sbp`, `... ed_ip`, `... hba1c --pc`,
then `python b_summarize.py` (~12–15 min per outcome at n_jobs=3; positive control ~45 min per
outcome).

Files: `b_data.py`, `b_analysis.py`, `b_summarize.py`, `run_all.sh`; `results.json`,
`average_effects_and_tests.csv`, `subgroups.csv`, `weight_x_prior_ed.csv`, `out/hba1c.json`,
`out/weight.json`, logs.

## Main numbers
Analytic N match the plan: HbA1c 341 (165 treated), weight 392 (190), SBP 371 (182), ED/IP 456
(225).

| | HbA1c (pp) | Weight (kg) |
|---|---|---|
| AIPW (95% CI), p | −0.085 (−0.479 to 0.308), 0.67 | +0.888 (−0.282 to 2.059), 0.14 |
| Median-over-repeats AIPW | −0.091 (−0.488 to 0.306) | +0.873 (−0.304 to 2.049) |
| AIPW with full-465 stratum propensity | −0.095 | +0.880 |
| Unadjusted difference | −0.177 | +0.959 |
| BLP β2 (median), 1-sided p | −4.39, 0.997 | −0.65, 0.72 |
| AUTOC (median), 1-sided p | −0.53, 0.998 | −0.15, 0.61 |
| Heterogeneity detected | No | No |
| Joint interaction HC1 Wald | χ² 15.2, df 21, p 0.81 | χ² 30.2, df 22, p 0.11 |
| Interaction randomization p (2000) | 0.913 | 0.298 |
| Joint test of AIPW scores on covariates | 0.83 | 0.099 |

Subgroups with Holm p < 0.10: only weight × prior-year ED visit: no ED (n 222, 104 treated)
+2.54 kg (1.08 to 4.00); ED (n 170, 86 treated) −1.26 (−3.13 to 0.61); difference −3.80,
interaction p 0.002, Holm 0.012. OLS with stratum dummies: −3.68 kg (−6.07 to −1.30), HC1
p 0.0027; within-stratum permutation p 0.002 (2000 permutations) under both the sharp null and
the constant-effect null. Next: weight × insulin p 0.025, Holm 0.15. HbA1c: all Holm p ≥ 0.86
(smallest raw: prior inpatient 0.12).

## Data quirks
1. Triglyceride placeholder 9999999 → missing.
2. Posted kg column = lb/2.2222; recomputed lb/2.20462; baseline kg is the 22nd modifier for weight.
3. Stratum flag `hba1c_cat_control` disagrees with baseline HbA1c ≥ 9.5 for 4 patients; the flag
   (used for randomization) was used.
4. `drug_any_db_0mo` and `drug_any_db_summary_0mo` disagree for 25 patients.
5. Missing baselines: BMI 6, BP 8, either 9.
6. Treated share per stratum 0.47–0.50 in all 465; 0.449–0.519 in the analytic samples (plan's
   0.44–0.53 matches the latter).

## Decisions where the plan was ambiguous
1. Propensity = treated share within stratum of the analytic sample (full-465 version as
   sensitivity: ±0.01).
2. Modifiers: days = `consent_date − hba1c_date_0mo`; outpatient = `op_n_0mo`; ED/IP =
   `ed_or_ip_n_0mo`; any IP = `ip_i_0mo`; any ED = `ed_i_0mo` (assumed "_0mo" = prior year; the
   dictionary says "at baseline"); any diabetes drug = `drug_any_db_0mo`.
3. AIPW outcome models: RF 500 trees, min leaf 5, max_features 1/3, per arm; same 5×5 folds as
   the forest; primary = patient AIPW scores averaged over repeats, SE = sd/√n.
4. Causal forest: same RF outcome model; treatment model returns the known stratum share
   (residuals verified = W − e); inner cv 5, 4000 trees, min leaf 10, max_samples 0.45, honest.
5. BLP: AIPW scores on [1, centred OOF CATE], HC1, 1-sided normal p; plain median over repeats
   (no ×2).
6. RATE/AUTOC: rank by predicted CATE, highest first (grf convention); AUTOC = mean over j of
   (top-j mean − overall mean); 1000 bootstraps with ranking fixed; median over repeats.
7. Detected only if both median 1-sided p < 0.05.
8. Interaction test: Y ~ W + Xc + W·Xc (Xc centred), strata through main effects only; HC1 Wald;
   randomization p from 2000 within-stratum permutations under a constant plug-in effect;
   p = (1+count)/(1+B).
9. Subgroups: means of full-sample AIPW scores within level (not refitted); 2-sided z test of the
   difference; Holm across 7 per outcome.
10. Weight × prior ED: OLS with treatment, ED, product, 3 stratum dummies, HC1, t p; permutation
    under both nulls with |t|.
11. Positive control (coded, not run): older = age ≥ analytic-sample median; inject δ × SD into
    treated older patients, δ = 0.25/0.5/1.
12. Seeds: base 20261029 (folds +1000·rep; AIPW forests +10·rep+k; causal forests +100·rep+k;
    bootstrap +rep; rerandomization +77; permutations +99).

Note from the agent: the HbA1c BLP β2 is strongly negative (−4.4) and GATES run opposite to the
predicted ordering, the pattern of out-of-fold predictions that are noise, not a real negative
signal.
