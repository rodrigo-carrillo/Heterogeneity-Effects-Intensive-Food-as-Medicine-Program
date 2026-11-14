# Agent C: original pipeline vs B's blind rebuild (2 Nov 2026)

Saved by the main session from Agent C's final message (the agent's tool was blocked from writing
this file). Content is the agent's, condensed. Files: `comparison.csv` (265 rows,
`c_build_comparison.py`), `c_tests.csv`, `c_collinearity.csv` (`c_reconcile_tests.py`),
`c_aipw_hba1c.csv` (`c_reconcile_aipw.py`), `c_ldl_interaction.csv` (`c_ldl_interaction.py`).

**Verdict:** B's rebuild confirms the original analysis and every manuscript conclusion: HbA1c
null (−0.05 vs −0.09); no outcome meets the heterogeneity criterion (all 1-sided BLP/AUTOC
p ≥ .72 in both); weight × prior ED nominal in both (Holm .013 vs .012; interaction p .0019 vs
.0017; OLS identical, −3.683, SE 1.218; permutation p .002 in both); positive control 1 SD
detected for all 4 outcomes in both, 0.5 SD for none. No computational bug. Same samples and
outcome values (N 341/392/371/456; treated 165/190/182/225; all 56 subgroup sizes).

**Tally (265 rows):** 197 agree; 65 different numbers, same conclusion; 2 conclusion differs
(nominal only, neither feeds a decision: SBP interaction asymptotic p .076 vs .021; SBP ×
metformin raw p .032 vs .052, Holm .23 vs .36); 1 not run by B (trim > 10 kg check).

## Condensed comparison (original → B)
| | HbA1c | Weight | SBP | ED/IP |
|---|---|---|---|---|
| AIPW | −0.053 → −0.085 | 0.90 → 0.89 | −0.15 → +0.25 | −6.04 → −5.83 |
| BLP β2, 1-sided p | −3.48, .999 → −4.39, .997 | −0.74, .80 → −0.65, .72 | −1.65, .84 → −3.13, .98 | −0.96, .75 → −2.25, .97 |
| AUTOC, p | −0.58, .998 → −0.53, .998 | −0.13, .59 → −0.15, .61 | −2.11, .83 → −3.43, .96 | −4.41, .85 → −7.26, .96 |
| Detected | No/No | No/No | No/No | No/No |
| Interaction asymptotic p | .63 → .81 | .059 → .11 | .076 → .021 | .079 → .077 |
| Interaction randomization p | .80 → .91 | .16 → .30 | .17 → .077 | .19 → .21 |
| AIPW scores on covariates p | .71 → .83 | .061 → .099 | .047 → .027 | .019 → .0066 |
| Smallest Holm p (7 variables) | 1.0 → .86 (prior IP) | .013 → .012 (prior ED) | .23 → .36 | .78 → .63 |
| Positive control 0.25/0.5/1 SD | N/N/Y both | N/N/Y both | N/N/Y both | N/N/Y both |

CI-excludes-0 status identical in every subgroup cell.

## Causes of differences
1. **Interaction-test covariate set** (SBP randomization .173 vs .077, and others). Original
   (`hte_food_as_medicine.py` §8, matrix `Z`) drops `drug_any_db_0mo` and `hba1c_cat_control`
   (and baseline weight for weight) as "near-collinear"; they are not (VIF 1.4 and 2.2–2.3; full
   rank). B keeps all modifiers, as the plan says ("joint test of all treatment × covariate
   interactions"). Reproduced both sides exactly; adding back `hba1c_cat_control` alone moves SBP
   .173 → .076. Original uses F, B χ² (SBP with B's set: F .026, χ² .021). **The original's
   omission is an unplanned deviation with a wrong justification.** No conclusion changes (all
   randomization p ≥ .05; not part of the decision rule). LDL runs robust (SBP .017 vs .015;
   ED/IP .045 vs .046).
2. **AIPW scores on covariates:** same covariate/F-vs-χ² issue plus different AIPW scores; B also
   does not drop baseline weight for weight (the plan requires it here; the original follows the
   plan). Not reported in the manuscript.
3. **AIPW outcome model `max_features`:** original sklearn default 1.0, B 1/3. Reproduced 2×2:
   HbA1c −0.053 (orig) / −0.087 (orig folds, 1/3) / −0.037 (B folds, 1.0) / −0.085 (B). ~0.17 SE;
   plan silent; both defensible.
4. **Forest outputs:** forest randomness, outcome model, propensity handling (B fixed stratum
   table; original refits shares within inner folds); B's AUTOC uses the mirror sign convention.
   All valid.
5. **Positive control:** original ±δ/2 for treated older/younger, 3 repeats; B +δ for older
   treated, 5 repeats; same contrast. HbA1c 0.5 SD: original BLP p .064, B .96 (neither detected).

## To fix
Cause 1: code comment, Table note ("collinear terms omitted") and eMethods 2 ("2 collinear
modifiers") are wrong. Either rerun the test with all modifiers (follows the plan; SBP
randomization p → .08) or describe the omission as an analyst choice departing from the plan.
Numbers that would change under B's choices (no conclusion changes): HbA1c AIPW −0.09 (−0.48 to
0.31); interaction P .91/.30/.08/.21; smallest Holm P HbA1c .86 (prior inpatient), SBP .36;
weight × ED CIs 1.08 to 4.00 and −3.13 to 0.60.
