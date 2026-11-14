# HTE results: food-as-medicine trial, HbA1c change at 6 months (auto-generated)

Run `hba1c_ldl_cc` (full): LDL in main model: True;
complete case (no imputation); seed 20260621; 4000 trees; 5-fold cross-fitting
x 5 repeats. Sign: negative = program lowered HbA1c more.

## Bottom line
Formal out-of-fold tests do NOT detect treatment-effect heterogeneity. Subgroup estimates below should be read as consistent with a uniform (null) effect; any nominally significant subgroup is expected under multiple testing.

## Average effect (N=285)
| Estimator | Effect on HbA1c change, pp (95% CI) |
|---|---|
| Paper, adjusted (longitudinal FE; N=349) | -0.10 (-0.46 to +0.25) |
| Paper model replicated (first differences) | -0.10 (-0.49 to +0.29) |
| ANCOVA, prespecified covariates | -0.06 (-0.46 to +0.34) |
| AIPW (cross-fitted, known propensity) | -0.06 (-0.48 to +0.35) |
| Causal forest | -0.06 (-0.55 to +0.44) |

## Omnibus heterogeneity tests (out-of-fold)
| Test | Estimate | p |
|---|---|---|
| BLP calibration beta2 (H1: >0) | -3.41 (SE 1.72) | 0.979 (one-sided) |
| GATES G1 - G3 | +1.04 (SE 0.53) | 0.049 |
| RATE AUTOC (H1: >0) | -0.339 (SE 0.178) | 0.972 (one-sided) |
| Joint BLP on covariates | - | 0.473 |
| ANCOVA joint interactions (HC1 Wald) | - | 0.438 |
| Same, randomization p (constant-effect null) | - | 0.714 |

## Pre-specified subgroups (Doyle family: smallest Holm-adjusted interaction p = 0.653)
| family               | subgroup                         | level                 |   n |    est |     lo |    hi |   paper_model_est |   p_interaction |   p_holm |
|:---------------------|:---------------------------------|:----------------------|----:|-------:|-------:|------:|------------------:|----------------:|---------:|
| Doyle pre-specified  | Baseline HbA1c                   | <9.5%                 | 140 | -0.216 | -0.647 | 0.214 |            -0.235 |           0.51  |    1     |
| Doyle pre-specified  | Baseline HbA1c                   | >=9.5%                | 145 |  0.059 | -0.636 | 0.754 |            -0.041 |           0.51  |    1     |
| Doyle pre-specified  | Site                             | Urban                 | 200 | -0.246 | -0.761 | 0.269 |            -0.269 |           0.183 |    1     |
| Doyle pre-specified  | Site                             | Rural                 |  85 |  0.323 | -0.336 | 0.983 |             0.306 |           0.183 |    1     |
| Doyle pre-specified  | Sex                              | Female                | 156 |  0.241 | -0.332 | 0.814 |             0.282 |           0.093 |    0.653 |
| Doyle pre-specified  | Sex                              | Male                  | 129 | -0.46  | -1.045 | 0.125 |            -0.436 |           0.093 |    0.653 |
| Doyle pre-specified  | Prior-year inpatient admission   | No                    | 218 |  0.091 | -0.368 | 0.551 |             0.109 |           0.17  |    1     |
| Doyle pre-specified  | Prior-year inpatient admission   | Yes                   |  67 | -0.621 | -1.527 | 0.286 |            -0.634 |           0.17  |    1     |
| Doyle pre-specified  | Prior-year ED visit              | No                    | 164 |  0.165 | -0.356 | 0.686 |             0.222 |           0.187 |    1     |
| Doyle pre-specified  | Prior-year ED visit              | Yes                   | 121 | -0.403 | -1.065 | 0.26  |            -0.387 |           0.187 |    1     |
| Doyle pre-specified  | Insulin prescribed               | No                    | 143 | -0.094 | -0.655 | 0.467 |            -0.007 |           0.932 |    1     |
| Doyle pre-specified  | Insulin prescribed               | Yes                   | 142 | -0.058 | -0.663 | 0.547 |            -0.156 |           0.932 |    1     |
| Doyle pre-specified  | Metformin prescribed             | No                    | 134 | -0.04  | -0.608 | 0.528 |            -0.188 |           0.87  |    1     |
| Doyle pre-specified  | Metformin prescribed             | Yes                   | 151 | -0.108 | -0.702 | 0.485 |            -0.207 |           0.87  |    1     |
| Added (not in Doyle) | Age                              | <65 y                 | 222 |  0.086 | -0.405 | 0.577 |            -0.028 |           0.086 |    0.257 |
| Added (not in Doyle) | Age                              | >=65 y                |  63 | -0.647 | -1.325 | 0.03  |            -0.324 |           0.086 |    0.257 |
| Added (not in Doyle) | Race/ethnicity                   | Hispanic or non-White |  57 |  0.641 | -0.45  | 1.732 |             0.318 |           0.134 |    0.269 |
| Added (not in Doyle) | Race/ethnicity                   | Non-Hispanic White    | 228 | -0.256 | -0.69  | 0.179 |            -0.286 |           0.134 |    0.269 |
| Added (not in Doyle) | Baseline A1c test before consent | <=90 d                | 161 |  0.015 | -0.551 | 0.581 |            -0.068 |           0.619 |    0.619 |
| Added (not in Doyle) | Baseline A1c test before consent | >90 d                 | 124 | -0.194 | -0.792 | 0.404 |            -0.191 |           0.619 |    0.619 |

## Honest depth-2 tree (estimation half; leaf-difference p = 0.404)
| level                                                   |   n |   est_discovery |    est |     lo |    hi |
|:--------------------------------------------------------|----:|----------------:|-------:|-------:|------:|
| Baseline HbA1c (%) <= 9.75 & Baseline HbA1c (%) <= 8.65 |  36 |           0.407 |  0.451 | -0.393 | 1.295 |
| Baseline HbA1c (%) <= 9.75 & Baseline HbA1c (%) > 8.65  |  42 |          -1.799 | -0.416 | -1.373 | 0.541 |
| Baseline HbA1c (%) > 9.75                               |  65 |           1.033 | -0.069 | -1.241 | 1.102 |

## Sensitivity analyses
| analysis                                                  |   n |   ate_aipw |   blp_beta2 |   blp_beta2_se |   blp_p_one_sided |   autoc |   autoc_se |   autoc_p_one_sided |
|:----------------------------------------------------------|----:|-----------:|------------:|---------------:|------------------:|--------:|-----------:|--------------------:|
| Main                                                      | 285 |     -0.062 |      -3.415 |          1.722 |             0.979 |  -0.339 |      0.178 |               0.972 |
| S1 Outcome = 6-mo level                                   | 285 |     -0.092 |      -2.07  |          1.806 |             0.874 |  -0.194 |      0.171 |               0.871 |
| S2 LDL removed                                            | 285 |     -0.045 |      -2.021 |          1.548 |             0.904 |  -0.263 |      0.183 |               0.925 |
| S3 IPW for excluded patients                              | 285 |     -0.044 |      -2.193 |          1.734 |             0.916 |  -0.308 |      0.17  |               0.965 |
| S4 Baseline measured <=180 d before consent               | 232 |     -0.051 |      -2.836 |          1.656 |             0.955 |  -0.331 |      0.195 |               0.955 |
| PC injected contrast 0.55 pp = 0.25 SD (older vs younger) | 285 |     -0.07  |      -0.176 |          1.163 |             0.56  |  -0.07  |      0.178 |               0.654 |
| PC injected contrast 1.10 pp = 0.5 SD (older vs younger)  | 285 |     -0.095 |       1.559 |          0.796 |             0.024 |   0.271 |      0.166 |               0.052 |
| PC injected contrast 2.19 pp = 1 SD (older vs younger)    | 285 |     -0.154 |       1.835 |          0.344 |             0     |   0.814 |      0.17  |               0     |

## Root-split stability (200 half-samples)
|                                    |   proportion |
|:-----------------------------------|-------------:|
| Baseline HbA1c (%)                 |        0.28  |
| Age (y)                            |        0.14  |
| Systolic BP                        |        0.105 |
| LDL                                |        0.085 |
| Days, baseline A1c test to consent |        0.065 |
| BMI                                |        0.065 |

Figures: fig1-fig6 in this folder. Full console log: console_log.txt.
