# Phase 2: code review

Files reviewed:
- `code/hte_food_as_medicine.py` (read at about 17:2x; the file was being edited concurrently and its mtime is now 17:30)
- `code/hte_lib.py`
- `code/run_all.py`
- `docs/04_results_across_outcomes.md`
- `outputs/v3/*/results.json`

Line numbers refer to the current on-disk version unless marked "(as read)". I fixed nothing.

## Summary
- I found no leakage, no post-randomization covariates, no wrong propensity, and no AIPW sign or unit errors in the HTE machinery. The out-of-fold design is sound and my blind re-implementation reproduces its key numbers.
- The material problem is an **interpretive error about the published weight result**, which persists in the doc and the v3 outputs.
- Secondary problems:
  - a "survives correction for 4 outcomes" claim that rests on the most favourable of several equally valid p-values;
  - mixed-vintage output files;
  - a few minor statistical shortcuts.

## Issues

### 1. MODERATE (substantive claim wrong): "the adjusted +1.95 is in pounds"
**Where the claim appears:**
- `docs/04_results_across_outcomes.md:29, 41-49, 143-144`
- `outputs/v3/summary_across_outcomes.md:24`
- `outputs/v3/weight/results_summary.md:13` ("published as kg; reproduces in lb")
- The pre-edit code (as read): `OUTCOMES["weight"]["paper_note"]`, the §2 comment, and the §2 `lb_m` check (Y_lb ~ treat + weight_lb_0 + d_time), which printed "the published weight differences match pounds".

**Evidence that the adjusted column is kg:**
- The do-file runs on `weight_val`, which is kg (×0.45 of lb).
- It multiply-imputes the 5 missing baselines with `mi impute mvn` in long format. The imputation model does not condition on the patient's own other-period weight.
- My emulation of that procedure (`v2_weight_mi.py`, `v2b_weight_mi_12mo.py`, 100 seeds) gives:

| Emulation | Estimate | Published |
|---|---|---|
| kg, 6 mo | **1.78 ± 1.75** | 1.95 (0.07 to 3.83) |
| kg, 12 mo | **1.75 ± 1.96** | 1.88 (−0.14 to 3.89) |
| lb, 6 mo | 3.97 ± 3.84 | — |
| lb, 12 mo | 3.91 ± 4.37 | — |

- The point estimate and the CI width both match kg.
- The analyst's lb argument rests on a coincidence: the complete-case lb estimate is 1.85, but its CI (−0.71 to 4.42) does not match the published CI width.

**What is right:** only the *unadjusted* −8.14 (−20.57 to 4.30) is in lb, and that reproduces exactly.

**Consequence:**
- The doc's conclusion that the kg effect is ≈+0.8–0.9 kg and n.s. is correct, but only as a complete-case finding.
- The paper's 1.95 kg is a kg number. It is inflated by 5 imputed baselines, which add about ±30 kg of noise to Δweight. It is not a unit error.
- "Raise with authors: weight differences are in pounds" (doc §2, §7c) should become: "unadjusted column in lb; adjusted column driven by imputation of 5 baselines without using the patient's other weight".
- The same MI explanation fits SBP (doc table: −0.74 vs −0.38, "paper imputes 6 baselines"). The doc applies it to SBP but not to weight, which is inconsistent.

**Status:** the analyst appears to have reached this independently during my review. The current code lines 212-214, 329-332 and 360-362 and `run_all.py:113-114` now say "adjusted kg with 5 imputed baselines". The doc and all `outputs/v3` summaries are **stale** and still state the lb claim.

### 2. MODERATE: the ×4-outcome survival of weight × ED is fragile and method-dependent
- **Where:** `docs/04_results_across_outcomes.md:17, 86, 95`.
- **The claim:** "survives Holm, and a further correction for 4 outcomes" (0.048).
- **Why it is fragile:** this uses the pipeline's AIPW interaction p (0.00178, `outputs/v3/weight/table_subgroups.csv`). Other equally valid estimators give:

| Method | Raw p | Holm (×7) | Further ×4 |
|---|---|---|---|
| Pipeline AIPW | 0.00178 | 0.012 | 0.048 |
| My AIPW | 0.0025 | 0.018 | 0.07 |
| OLS HC1 | 0.0026 | 0.018 | 0.07 |
| Within-stratum permutation (10k) | 0.0031 | 0.022 | 0.087 |
| Pipeline's own console re-run (17:31) | 0.002 | 0.012 | — |

- The subgroup point estimates are robust: about +2.5 vs −1.2 kg, interaction ≈ −3.6 kg in every method.
- But "survives ×4" holds only at the knife edge for one estimator.
- The family definition is also generous. 7 subgroups × 4 outcomes × 2 covariate sets were examined, plus 3 added subgroups and data-driven splits.
- Recommendation: report it as nominal/Holm-within-outcome only, and add the randomization p.

### 3. MINOR-MODERATE: outputs/v3 mixes runs of different code versions
- **Where:** `outputs/v3/weight/`.
- `results.json`, `results_summary.md`, `table_subgroups.csv` and `table_sensitivity.csv` are timestamped 14:32.
- `console_log.txt`, `table_balance.csv`, `table_blp_covariates.csv`, `table_variable_importance.csv` and `table_split_stability.csv` are timestamped 17:30-17:31, from a newer run that was in progress or did not finish.
- They disagree slightly. For example, prior-ED "No" is 2.485 (1.046 to 3.924) in the CSV and 2.508 (1.057 to 3.959) in the console log.
- The doc's numbers match the 14:32 files.
- Rerun `run_all.py` to completion before citing anything.

### 4. MINOR: weight is not exactly kg
- **Where:** `code/run_all.py:114` ("Our analysis is in kg (lb / 2.20462)"), and the `LABELS`/`UNIT` "kg" used throughout.
- In the data `weight_val = weight_lb_val × 0.45` exactly (ratio 2.2222 for all 459 rows).
- So the analysis unit is 0.8% smaller than true kg. That is immaterial numerically (2.49 → 2.51 kg), but the stated conversion is false.
- The published Table 3 means use ×0.4536 while its SDs use ×0.45. This is an inconsistency in the paper and does not affect this pipeline.

### 5. MINOR: medians of p-values without the Chernozhukov adjustment
- **Where:** `hte_food_as_medicine.py:578` (`BLP = blp_df.median()`), `:596` (AUTOC median), and `omnibus()` at `:909-919`.
- Chernozhukov–Demirer–Duflo–Fernández-Val recommend reporting 2 × the median p (and adjusted-level CIs) when aggregating over sample splits.
- Using the plain median is anti-conservative.
- It does not change any conclusion, because all omnibus tests are null. It would matter if a positive were claimed.

### 6. MINOR: the S3 "IPW for excluded patients" is only partially weighted
- **Where:** `omnibus()` at `:909-916`, called at `:961`.
- The weights enter the BLP regression and the ATE average.
- They do not enter the cross-fitted forests or the AIPW nuisance fits (`crossfit` is unweighted), or the RATE/AUTOC bootstrap.
- The doc's "S3 … changes no conclusion" should be read with that in mind.

### 7. MINOR: interpreting a negative BLP slope
- **Where:** doc `:64-66`; code `:584-586` (print text).
- "A significantly NEGATIVE beta2 indicates the forest is fitting noise" is a reasonable heuristic, not a theorem.
- With proper K-fold out-of-fold CATEs the slope should be about 0 under no heterogeneity. A t of about −2.7 (−4.2/1.5) is unusual and could be chance in this one dataset; repeats share the data, so their agreement is not independent evidence.
- My grf OOB `test_calibration` gives −14.5. That is the familiar OOB artifact and is much larger than the K-fold value.
- Either way the conclusion "no detectable heterogeneity" stands.

### 8. MINOR: the RATE p-value is sensitive to how priorities are cross-fitted
- **Where:** `rate_autoc` at `:558`, used at `:572`.
- The pipeline uses 5-fold out-of-fold priorities scored against same-repeat AIPW scores and reports p=0.998 for HbA1c.
- With a split-half train/evaluate RATE (grf), I get one-sided p of 0.43–0.82.
- Both are null. "p≈0.998" should not be read as strong evidence *against* heterogeneity. It only reflects the negative-slope phenomenon in item 7.

### 9. MINOR: the "unadjusted" weight comparison uses the analytic sample
- **Where:** `:350-352` (`unadj` uses `df`, the N=395 analytic sample).
- It prints "Unadjusted difference in 6-mo level: −4.78 kg   paper: −8.14", which compares N=395 against the paper's N=400, and kg against lb.
- The separate `lb_u` line, which uses `raw`, is the right comparison.
- Mislabelled comparison only.

## Checked and found correct
- **Analysis restrictions:** 6-mo only, `treat_now`, strata = 2·rural + hba1c_cat (`:280`). Samples: HbA1c 349/170, weight 395/192 (both weights required, `:286-294`), ED/IP 465.
- **Propensity:** the known stratum share is used everywhere, both for `E_HAT` (`:310`) and inside the forest (`StratumPropensity`, `hte_lib.py:21-46`).
- **AIPW score formula** (`hte_lib.py:62-73`): correct signs. Arm-specific RF nuisances are fit on training folds only. There is no leakage between folds (`crossfit_repeat`, `hte_lib.py:76-86`), and the forest is refit per fold.
- **Covariates:** baseline only (`:228-238`). No 6-mo medications, survey indices or engagement variables are used. `days_a1c_to_consent` uses baseline dates.
- **Paper-model replication:** the FD form (`:341`), and the binary form with y0 as a regressor (`:346`), equal the FE xtreg treatment coefficient. Verified independently: HbA1c −0.102, ED/IP −7.129.
- **Subgroup AIPW and Holm:** the subgroup effect is the mean of the AIPW scores within the subgroup, and the interaction is the difference of independent subgroup means (`:723-745`). Holm is applied within the Doyle family of 7 (`:748`). These reproduce my numbers.
- **BLP:** `blp_test` regresses Γ on [1, τ̂ − mean τ̂] with HC1, which is the score-based BLP. My 5-fold cross-fit BLP with grf gives a slope of −4.27, matching the reported −4.23.
- **Randomization test of the joint interactions** (`:662-688`): treatment is re-randomized within strata under a constant-effect null. This is correct as described.
