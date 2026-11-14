# HTE results: food-as-medicine trial, SBP change at 6 months (auto-generated)

Run `sbp_ldl_cc` (full): LDL in main model: True;
complete case (no imputation); seed 20260621; 4000 trees; 5-fold cross-fitting
x 5 repeats. Sign: negative = program lowered SBP more.

## Bottom line
Formal out-of-fold tests do NOT detect treatment-effect heterogeneity. Subgroup estimates below should be read as consistent with a uniform (null) effect; any nominally significant subgroup is expected under multiple testing.

## Average effect (N=308)
| Estimator | Effect on SBP change, mm Hg (95% CI) |
|---|---|
| Paper, adjusted (longitudinal FE; N=378) | -0.38 (-4.16 to +3.40) |
| Paper model replicated (first differences) | +0.39 (-3.71 to +4.49) |
| ANCOVA, prespecified covariates | +0.96 (-3.38 to +5.30) |
| AIPW (cross-fitted, known propensity) | +1.10 (-3.21 to +5.42) |
| Causal forest | +1.47 (-3.75 to +6.68) |

## Omnibus heterogeneity tests (out-of-fold)
| Test | Estimate | p |
|---|---|---|
| BLP calibration beta2 (H1: >0) | -3.13 (SE 1.73) | 0.958 (one-sided) |
| GATES G1 - G3 | +7.41 (SE 5.35) | 0.166 |
| RATE AUTOC (H1: >0) | -2.979 (SE 1.941) | 0.938 (one-sided) |
| Joint BLP on covariates | - | 0.014 |
| ANCOVA joint interactions (HC1 Wald) | - | 0.002 |
| Same, randomization p (constant-effect null) | - | 0.015 |

## Pre-specified subgroups (Doyle family: smallest Holm-adjusted interaction p = 0.609)
| family               | subgroup                         | level                 |   n |    est |      lo |     hi |   paper_model_est |   p_interaction |   p_holm |
|:---------------------|:---------------------------------|:----------------------|----:|-------:|--------:|-------:|------------------:|----------------:|---------:|
| Doyle pre-specified  | Baseline HbA1c                   | <9.5%                 | 148 | -0.805 |  -6.946 |  5.337 |            -1.203 |           0.431 |    1     |
| Doyle pre-specified  | Baseline HbA1c                   | >=9.5%                | 160 |  2.642 |  -3.355 |  8.638 |             1.664 |           0.431 |    1     |
| Doyle pre-specified  | Site                             | Urban                 | 210 |  3.617 |  -1.41  |  8.644 |             2.902 |           0.087 |    0.609 |
| Doyle pre-specified  | Site                             | Rural                 |  98 | -4.653 | -12.678 |  3.373 |            -4.628 |           0.087 |    0.609 |
| Doyle pre-specified  | Sex                              | Female                | 170 |  0.612 |  -5.352 |  6.575 |             0.383 |           0.849 |    1     |
| Doyle pre-specified  | Sex                              | Male                  | 138 |  1.447 |  -4.712 |  7.605 |            -0.204 |           0.849 |    1     |
| Doyle pre-specified  | Prior-year inpatient admission   | No                    | 237 | -0.323 |  -5.285 |  4.639 |            -0.918 |           0.256 |    1     |
| Doyle pre-specified  | Prior-year inpatient admission   | Yes                   |  71 |  5.355 |  -3.097 | 13.808 |             4.722 |           0.256 |    1     |
| Doyle pre-specified  | Prior-year ED visit              | No                    | 176 |  2.03  |  -3.615 |  7.675 |             1.765 |           0.583 |    1     |
| Doyle pre-specified  | Prior-year ED visit              | Yes                   | 132 | -0.406 |  -7.018 |  6.205 |            -0.342 |           0.583 |    1     |
| Doyle pre-specified  | Insulin prescribed               | No                    | 160 |  4.521 |  -1.085 | 10.126 |             2.583 |           0.093 |    0.609 |
| Doyle pre-specified  | Insulin prescribed               | Yes                   | 148 | -2.836 |  -9.351 |  3.68  |            -2.774 |           0.093 |    0.609 |
| Doyle pre-specified  | Metformin prescribed             | No                    | 142 | -2.398 |  -8.781 |  3.984 |            -2.287 |           0.153 |    0.763 |
| Doyle pre-specified  | Metformin prescribed             | Yes                   | 166 |  3.881 |  -1.889 |  9.65  |             1.939 |           0.153 |    0.763 |
| Added (not in Doyle) | Age                              | <65 y                 | 238 |  0.52  |  -4.274 |  5.315 |            -0.237 |           0.707 |    1     |
| Added (not in Doyle) | Age                              | >=65 y                |  70 |  2.569 |  -6.992 | 12.129 |             3.354 |           0.707 |    1     |
| Added (not in Doyle) | Race/ethnicity                   | Hispanic or non-White |  56 |  5.996 |  -3.623 | 15.615 |             3.915 |           0.264 |    0.792 |
| Added (not in Doyle) | Race/ethnicity                   | Non-Hispanic White    | 252 | -0.128 |  -4.911 |  4.656 |            -0.3   |           0.264 |    0.792 |
| Added (not in Doyle) | Baseline A1c test before consent | <=90 d                | 184 |  2.053 |  -3.719 |  7.826 |             1.054 |           0.545 |    1     |
| Added (not in Doyle) | Baseline A1c test before consent | >90 d                 | 124 | -0.598 |  -6.948 |  5.751 |            -0.755 |           0.545 |    1     |

## Honest depth-2 tree (estimation half; leaf-difference p = 0.568)
| level                            |   n |   est_discovery |    est |     lo |     hi |
|:---------------------------------|----:|----------------:|-------:|-------:|-------:|
| Age (y) <= 67 & Metformin <= 0.5 |  63 |          -4.434 |  2.724 | -6.474 | 11.922 |
| Age (y) <= 67 & Metformin > 0.5  |  71 |           8.897 | -0.743 | -9.141 |  7.655 |
| Age (y) > 67                     |  20 |         -13.764 | 10.099 | -9.002 | 29.201 |

## Sensitivity analyses
| analysis                                                     |   n |   ate_aipw |   blp_beta2 |   blp_beta2_se |   blp_p_one_sided |   autoc |   autoc_se |   autoc_p_one_sided |
|:-------------------------------------------------------------|----:|-----------:|------------:|---------------:|------------------:|--------:|-----------:|--------------------:|
| Main                                                         | 308 |      1.103 |      -3.127 |          1.729 |             0.958 |  -2.979 |      1.941 |               0.938 |
| S1 Outcome = 6-mo level                                      | 308 |      0.814 |      -4.011 |          1.382 |             0.997 |  -5.647 |      1.843 |               0.999 |
| S2 LDL removed                                               | 308 |      0.941 |      -2.534 |          1.209 |             0.953 |  -4.345 |      1.955 |               0.987 |
| S3 IPW for excluded patients                                 | 308 |      1.135 |      -2.598 |          1.259 |             0.936 |  -3.628 |      1.985 |               0.966 |
| S4 Baseline measured <=180 d before consent                  | 296 |      1.284 |      -3.143 |          1.652 |             0.945 |  -2.569 |      2.24  |               0.874 |
| PC injected contrast 5.27 mm Hg = 0.25 SD (older vs younger) | 308 |      0.841 |      -1.729 |          1.171 |             0.865 |  -1.962 |      2.001 |               0.837 |
| PC injected contrast 10.55 mm Hg = 0.5 SD (older vs younger) | 308 |      0.952 |      -0.31  |          1.059 |             0.597 |  -0.456 |      1.921 |               0.594 |
| PC injected contrast 21.09 mm Hg = 1 SD (older vs younger)   | 308 |      1.18  |       1.81  |          0.589 |             0.001 |   5.676 |      1.955 |               0.002 |

## Root-split stability (200 half-samples)
|                                    |   proportion |
|:-----------------------------------|-------------:|
| BMI                                |        0.23  |
| Diastolic BP                       |        0.16  |
| Systolic BP                        |        0.13  |
| Baseline HbA1c (%)                 |        0.105 |
| Days, baseline A1c test to consent |        0.09  |
| LDL                                |        0.08  |

Figures: fig1-fig6 in this folder. Full console log: console_log.txt.
