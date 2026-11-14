# HTE results: food-as-medicine trial, SBP change at 6 months (auto-generated)

Run `sbp` (full): LDL in main model: False;
complete case (no imputation); seed 20260621; 4000 trees; 5-fold cross-fitting
x 5 repeats. Sign: negative = program lowered SBP more.

## Bottom line
Formal out-of-fold tests do NOT detect treatment-effect heterogeneity. Subgroup estimates below should be read as consistent with a uniform (null) effect; any nominally significant subgroup is expected under multiple testing.

## Average effect (N=371)
| Estimator | Effect on SBP change, mm Hg (95% CI) |
|---|---|
| Paper, adjusted (longitudinal FE; N=378) | -0.38 (-4.16 to +3.40) |
| Paper model replicated (first differences) | -0.59 (-4.31 to +3.14) |
| ANCOVA, prespecified covariates | -0.76 (-4.65 to +3.12) |
| AIPW (cross-fitted, known propensity) | -0.15 (-4.05 to +3.76) |
| Causal forest | +0.22 (-4.80 to +5.23) |

## Omnibus heterogeneity tests (out-of-fold)
| Test | Estimate | p |
|---|---|---|
| BLP calibration beta2 (H1: >0) | -1.65 (SE 1.66) | 0.835 (one-sided) |
| GATES G1 - G3 | +5.84 (SE 4.89) | 0.233 |
| RATE AUTOC (H1: >0) | -2.110 (SE 2.209) | 0.830 (one-sided) |
| Joint BLP on covariates | - | 0.035 |
| ANCOVA joint interactions (HC1 Wald) | - | 0.026 |
| Same, randomization p (constant-effect null) | - | 0.080 |

## Pre-specified subgroups (Doyle family: smallest Holm-adjusted interaction p = 0.227)
| family               | subgroup                         | level                 |   n |    est |      lo |     hi |   paper_model_est |   p_interaction |   p_holm |
|:---------------------|:---------------------------------|:----------------------|----:|-------:|--------:|-------:|------------------:|----------------:|---------:|
| Doyle pre-specified  | Baseline HbA1c                   | <9.5%                 | 167 | -2.122 |  -8.087 |  3.843 |            -2.037 |           0.383 |    1     |
| Doyle pre-specified  | Baseline HbA1c                   | >=9.5%                | 204 |  1.378 |  -3.752 |  6.508 |             0.509 |           0.383 |    1     |
| Doyle pre-specified  | Site                             | Urban                 | 253 |  1.807 |  -2.819 |  6.433 |             1.188 |           0.147 |    0.733 |
| Doyle pre-specified  | Site                             | Rural                 | 118 | -4.494 | -11.635 |  2.646 |            -4.217 |           0.147 |    0.733 |
| Doyle pre-specified  | Sex                              | Female                | 198 | -0.942 |  -6.379 |  4.496 |            -1.115 |           0.688 |    1     |
| Doyle pre-specified  | Sex                              | Male                  | 173 |  0.654 |  -4.925 |  6.234 |            -0.889 |           0.688 |    1     |
| Doyle pre-specified  | Prior-year inpatient admission   | No                    | 280 | -0.564 |  -5.085 |  3.956 |            -0.881 |           0.742 |    1     |
| Doyle pre-specified  | Prior-year inpatient admission   | Yes                   |  91 |  0.932 |  -6.754 |  8.618 |            -0.149 |           0.742 |    1     |
| Doyle pre-specified  | Prior-year ED visit              | No                    | 210 |  1.563 |  -3.493 |  6.619 |             0.955 |           0.315 |    1     |
| Doyle pre-specified  | Prior-year ED visit              | Yes                   | 161 | -2.494 |  -8.576 |  3.589 |            -2.594 |           0.315 |    1     |
| Doyle pre-specified  | Insulin prescribed               | No                    | 184 |  3.666 |  -1.491 |  8.824 |             2.26  |           0.053 |    0.315 |
| Doyle pre-specified  | Insulin prescribed               | Yes                   | 187 | -3.999 |  -9.785 |  1.786 |            -3.622 |           0.053 |    0.315 |
| Doyle pre-specified  | Metformin prescribed             | No                    | 178 | -4.611 | -10.367 |  1.145 |            -4.627 |           0.032 |    0.227 |
| Doyle pre-specified  | Metformin prescribed             | Yes                   | 193 |  3.873 |  -1.35  |  9.096 |             2.529 |           0.032 |    0.227 |
| Added (not in Doyle) | Age                              | <65 y                 | 287 |  0.042 |  -4.248 |  4.332 |            -0.585 |           0.836 |    1     |
| Added (not in Doyle) | Age                              | >=65 y                |  84 | -1.015 | -10.055 |  8.024 |             0.506 |           0.836 |    1     |
| Added (not in Doyle) | Race/ethnicity                   | Hispanic or non-White |  65 |  3.8   |  -5.457 | 13.057 |             2.422 |           0.352 |    1     |
| Added (not in Doyle) | Race/ethnicity                   | Non-Hispanic White    | 306 | -1.047 |  -5.339 |  3.246 |            -1.113 |           0.352 |    1     |
| Added (not in Doyle) | Baseline A1c test before consent | <=90 d                | 215 |  0.146 |  -5.241 |  5.532 |            -0.804 |           0.836 |    1     |
| Added (not in Doyle) | Baseline A1c test before consent | >90 d                 | 156 | -0.67  |  -6.223 |  4.882 |            -0.837 |           0.836 |    1     |

## Honest depth-2 tree (estimation half; leaf-difference p = 0.841)
| level                                                              |   n |   est_discovery |    est |      lo |     hi |
|:-------------------------------------------------------------------|----:|----------------:|-------:|--------:|-------:|
| Prior-year outpatient visits <= 4.5 & BMI <= 36.41                 |  33 |          -3.421 |  4.402 |  -9.398 | 18.201 |
| Prior-year outpatient visits <= 4.5 & BMI > 36.41                  |  30 |         -22.366 | -2.977 | -17.509 | 11.554 |
| Prior-year outpatient visits > 4.5 & Other glucose-lowering <= 0.5 |  83 |           9.013 |  0.311 |  -8.355 |  8.977 |
| Prior-year outpatient visits > 4.5 & Other glucose-lowering > 0.5  |  40 |          -5.806 | -3.265 | -15.268 |  8.739 |

## Sensitivity analyses
| analysis                                                     |   n |   ate_aipw |   blp_beta2 |   blp_beta2_se |   blp_p_one_sided |   autoc |   autoc_se |   autoc_p_one_sided |
|:-------------------------------------------------------------|----:|-----------:|------------:|---------------:|------------------:|--------:|-----------:|--------------------:|
| Main                                                         | 371 |     -0.148 |      -1.646 |          1.66  |             0.835 |  -2.11  |      2.209 |               0.83  |
| S1 Outcome = 6-mo level                                      | 371 |     -0.882 |      -0.721 |          1.565 |             0.668 |   1.16  |      2.162 |               0.296 |
| S3 IPW for excluded patients                                 | 371 |     -0.594 |      -1.239 |          1.58  |             0.781 |  -0.125 |      2.22  |               0.522 |
| S4 Baseline measured <=180 d before consent                  | 356 |     -0.123 |      -1.401 |          1.34  |             0.852 |  -2.062 |      2.046 |               0.843 |
| PC injected contrast 5.23 mm Hg = 0.25 SD (older vs younger) | 371 |     -0.435 |      -0.703 |          1.372 |             0.683 |  -0.279 |      2.046 |               0.554 |
| PC injected contrast 10.45 mm Hg = 0.5 SD (older vs younger) | 371 |     -0.575 |       0.335 |          1.172 |             0.4   |   1.205 |      1.968 |               0.27  |
| PC injected contrast 20.91 mm Hg = 1 SD (older vs younger)   | 371 |     -0.609 |       1.726 |          0.552 |             0     |   4.356 |      2.045 |               0.017 |

## Root-split stability (200 half-samples)
|                              |   proportion |
|:-----------------------------|-------------:|
| Diastolic BP                 |        0.205 |
| Systolic BP                  |        0.18  |
| Metformin                    |        0.135 |
| Baseline HbA1c (%)           |        0.085 |
| Prior-year outpatient visits |        0.08  |
| Sulfonylurea                 |        0.065 |

Figures: fig1-fig6 in this folder. Full console log: console_log.txt.
