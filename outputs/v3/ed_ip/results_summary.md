# HTE results: food-as-medicine trial, any ED/IP use at 6 months (auto-generated)

Run `ed_ip` (full): LDL in main model: False;
complete case (no imputation); seed 20260621; 4000 trees; 5-fold cross-fitting
x 5 repeats. Sign: negative = fewer patients with ED/IP use.

## Bottom line
Formal out-of-fold tests do NOT detect treatment-effect heterogeneity. Subgroup estimates below should be read as consistent with a uniform (null) effect; any nominally significant subgroup is expected under multiple testing.

## Average effect (N=456)
| Estimator | Effect on any ED/IP use, pp (95% CI) |
|---|---|
| Paper, adjusted (longitudinal FE; N=465) | -7.13 (-15.30 to +1.04) |
| Paper model replicated (first differences) | -7.08 (-15.32 to +1.15) |
| ANCOVA, prespecified covariates | -8.07 (-16.51 to +0.38) |
| AIPW (cross-fitted, known propensity) | -6.04 (-14.12 to +2.04) |
| Causal forest | -5.22 (-16.55 to +6.11) |

## Omnibus heterogeneity tests (out-of-fold)
| Test | Estimate | p |
|---|---|---|
| BLP calibration beta2 (H1: >0) | -0.96 (SE 1.39) | 0.749 (one-sided) |
| GATES G1 - G3 | -1.29 (SE 10.39) | 0.901 |
| RATE AUTOC (H1: >0) | -4.411 (SE 4.313) | 0.847 (one-sided) |
| Joint BLP on covariates | - | 0.016 |
| ANCOVA joint interactions (HC1 Wald) | - | 0.085 |
| Same, randomization p (constant-effect null) | - | 0.213 |

## Pre-specified subgroups (Doyle family: smallest Holm-adjusted interaction p = 0.782)
| family               | subgroup                         | level                 |   n |     est |      lo |     hi |   paper_model_est |   p_interaction |   p_holm |
|:---------------------|:---------------------------------|:----------------------|----:|--------:|--------:|-------:|------------------:|----------------:|---------:|
| Doyle pre-specified  | Baseline HbA1c                   | <9.5%                 | 199 |  -1.112 | -13.615 | 11.391 |            -1.258 |           0.293 |    1     |
| Doyle pre-specified  | Baseline HbA1c                   | >=9.5%                | 257 |  -9.864 | -20.364 |  0.636 |           -11.81  |           0.293 |    1     |
| Doyle pre-specified  | Site                             | Urban                 | 327 |  -4.707 | -14.155 |  4.741 |            -4.204 |           0.608 |    1     |
| Doyle pre-specified  | Site                             | Rural                 | 129 |  -9.434 | -24.847 |  5.979 |           -14.195 |           0.608 |    1     |
| Doyle pre-specified  | Sex                              | Female                | 249 | -10.534 | -21.572 |  0.505 |           -12.647 |           0.229 |    1     |
| Doyle pre-specified  | Sex                              | Male                  | 207 |  -0.644 | -12.386 | 11.097 |            -0.481 |           0.229 |    1     |
| Doyle pre-specified  | Prior-year inpatient admission   | No                    | 349 |  -6.828 | -15.752 |  2.096 |            -6.791 |           0.747 |    1     |
| Doyle pre-specified  | Prior-year inpatient admission   | Yes                   | 107 |  -3.49  | -21.722 | 14.743 |            -8.016 |           0.747 |    1     |
| Doyle pre-specified  | Prior-year ED visit              | No                    | 266 |  -2.94  | -12.557 |  6.677 |            -2.86  |           0.387 |    1     |
| Doyle pre-specified  | Prior-year ED visit              | Yes                   | 190 | -10.391 | -24.255 |  3.474 |           -10.594 |           0.387 |    1     |
| Doyle pre-specified  | Insulin prescribed               | No                    | 229 |  -0.912 | -12.664 | 10.839 |             0.321 |           0.209 |    1     |
| Doyle pre-specified  | Insulin prescribed               | Yes                   | 227 | -11.222 | -22.205 | -0.238 |           -14.319 |           0.209 |    1     |
| Doyle pre-specified  | Metformin prescribed             | No                    | 228 | -12.565 | -23.695 | -1.435 |           -13.815 |           0.112 |    0.782 |
| Doyle pre-specified  | Metformin prescribed             | Yes                   | 228 |   0.476 | -11.12  | 12.073 |             0.156 |           0.112 |    0.782 |
| Added (not in Doyle) | Age                              | <65 y                 | 359 |  -4.791 | -13.674 |  4.091 |            -5.05  |           0.579 |    0.579 |
| Added (not in Doyle) | Age                              | >=65 y                |  97 | -10.683 | -29.496 |  8.131 |           -13.713 |           0.579 |    0.579 |
| Added (not in Doyle) | Race/ethnicity                   | Hispanic or non-White |  86 |  11.042 |  -7.397 | 29.481 |            12.085 |           0.044 |    0.132 |
| Added (not in Doyle) | Race/ethnicity                   | Non-Hispanic White    | 370 | -10.016 | -18.927 | -1.105 |           -11.485 |           0.044 |    0.132 |
| Added (not in Doyle) | Baseline A1c test before consent | <=90 d                | 263 | -10.987 | -21.642 | -0.333 |           -12.038 |           0.159 |    0.317 |
| Added (not in Doyle) | Baseline A1c test before consent | >90 d                 | 193 |   0.691 | -11.565 | 12.947 |            -0.442 |           0.159 |    0.317 |

## Honest depth-2 tree (estimation half; leaf-difference p = 0.926)
| level                                                             |   n |   est_discovery |     est |      lo |     hi |
|:------------------------------------------------------------------|----:|----------------:|--------:|--------:|-------:|
| Prior-year outpatient visits <= 7.5 & Baseline HbA1c (%) <= 10.35 |  84 |         -20.217 | -12.955 | -31.914 |  6.004 |
| Prior-year outpatient visits <= 7.5 & Baseline HbA1c (%) > 10.35  |  53 |           4.348 | -11.881 | -36.031 | 12.269 |
| Prior-year outpatient visits > 7.5 & Baseline HbA1c (%) <= 9.55   |  44 |          33.571 |  -2.024 | -30.387 | 26.339 |
| Prior-year outpatient visits > 7.5 & Baseline HbA1c (%) > 9.55    |  47 |          -7.753 |  -6.967 | -33.958 | 20.024 |

## Sensitivity analyses
| analysis                                                   |   n |   ate_aipw |   blp_beta2 |   blp_beta2_se |   blp_p_one_sided |   autoc |   autoc_se |   autoc_p_one_sided |
|:-----------------------------------------------------------|----:|-----------:|------------:|---------------:|------------------:|--------:|-----------:|--------------------:|
| Main                                                       | 456 |     -6.042 |      -0.958 |          1.393 |             0.749 |  -4.411 |      4.313 |               0.847 |
| S1 Outcome = no. of ED/IP visits                           | 456 |     -0.094 |      -2.299 |          1.549 |             0.959 |  -0.262 |      0.17  |               0.939 |
| S3 IPW for excluded patients                               | 456 |     -6.544 |      -1.45  |          1.369 |             0.833 |  -5.878 |      4.181 |               0.92  |
| PC injected contrast 11.80 pp = 0.25 SD (older vs younger) | 456 |     -6.662 |      -0.742 |          1.229 |             0.699 |  -4.195 |      4.233 |               0.839 |
| PC injected contrast 23.60 pp = 0.5 SD (older vs younger)  | 456 |     -6.727 |       0.684 |          0.881 |             0.232 |  -0.513 |      4.283 |               0.548 |
| PC injected contrast 47.19 pp = 1 SD (older vs younger)    | 456 |     -6.593 |       1.383 |          0.39  |             0     |  10.258 |      4.176 |               0.007 |

## Root-split stability (200 half-samples)
|                                    |   proportion |
|:-----------------------------------|-------------:|
| Systolic BP                        |        0.315 |
| Baseline HbA1c (%)                 |        0.175 |
| Days, baseline A1c test to consent |        0.17  |
| Prior-year outpatient visits       |        0.085 |
| Age (y)                            |        0.085 |
| BMI                                |        0.06  |

Figures: fig1-fig6 in this folder. Full console log: console_log.txt.
