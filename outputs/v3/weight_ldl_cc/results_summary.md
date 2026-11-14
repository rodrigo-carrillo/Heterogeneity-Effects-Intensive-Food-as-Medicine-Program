# HTE results: food-as-medicine trial, weight change at 6 months (auto-generated)

Run `weight_ldl_cc` (full): LDL in main model: True;
complete case (no imputation); seed 20260621; 4000 trees; 5-fold cross-fitting
x 5 repeats. Sign: negative = less weight gain with program.

## Bottom line
Formal out-of-fold tests do NOT detect treatment-effect heterogeneity. Subgroup estimates below should be read as consistent with a uniform (null) effect; any nominally significant subgroup is expected under multiple testing.

## Average effect (N=321)
| Estimator | Effect on weight change, kg (95% CI) |
|---|---|
| Paper, adjusted (longitudinal FE; N=400; unadjusted published in lb; adjusted kg with 5 imputed baselines) | +1.95 (+0.07 to +3.83) |
| Paper model replicated (first differences) | +0.64 (-0.60 to +1.88) |
| ANCOVA, prespecified covariates | +0.62 (-0.66 to +1.91) |
| AIPW (cross-fitted, known propensity) | +0.75 (-0.49 to +2.00) |
| Causal forest | +0.89 (-0.87 to +2.65) |

## Omnibus heterogeneity tests (out-of-fold)
| Test | Estimate | p |
|---|---|---|
| BLP calibration beta2 (H1: >0) | -0.91 (SE 1.23) | 0.781 (one-sided) |
| GATES G1 - G3 | +0.67 (SE 1.52) | 0.657 |
| RATE AUTOC (H1: >0) | +0.249 (SE 0.673) | 0.356 (one-sided) |
| Joint BLP on covariates | - | 0.005 |
| ANCOVA joint interactions (HC1 Wald) | - | 0.012 |
| Same, randomization p (constant-effect null) | - | 0.067 |

## Pre-specified subgroups (Doyle family: smallest Holm-adjusted interaction p = 0.009)
| family               | subgroup                         | level                 |   n |    est |     lo |    hi |   paper_model_est |   p_interaction |   p_holm |
|:---------------------|:---------------------------------|:----------------------|----:|-------:|-------:|------:|------------------:|----------------:|---------:|
| Doyle pre-specified  | Baseline HbA1c                   | <9.5%                 | 151 |  1.419 | -0.367 | 3.205 |             1.027 |           0.32  |    1     |
| Doyle pre-specified  | Baseline HbA1c                   | >=9.5%                | 170 |  0.16  | -1.561 | 1.88  |             0.181 |           0.32  |    1     |
| Doyle pre-specified  | Site                             | Urban                 | 223 |  1.157 | -0.351 | 2.664 |             0.853 |           0.325 |    1     |
| Doyle pre-specified  | Site                             | Rural                 |  98 | -0.169 | -2.339 | 2.001 |            -0.031 |           0.325 |    1     |
| Doyle pre-specified  | Sex                              | Female                | 179 |  1.799 |  0.002 | 3.597 |             2.035 |           0.056 |    0.334 |
| Doyle pre-specified  | Sex                              | Male                  | 142 | -0.568 | -2.196 | 1.059 |            -0.984 |           0.056 |    0.334 |
| Doyle pre-specified  | Prior-year inpatient admission   | No                    | 247 |  1.159 | -0.213 | 2.531 |             1.044 |           0.269 |    1     |
| Doyle pre-specified  | Prior-year inpatient admission   | Yes                   |  74 | -0.607 | -3.418 | 2.205 |             0.093 |           0.269 |    1     |
| Doyle pre-specified  | Prior-year ED visit              | No                    | 182 |  2.531 |  0.996 | 4.066 |             2.297 |           0.001 |    0.009 |
| Doyle pre-specified  | Prior-year ED visit              | Yes                   | 139 | -1.577 | -3.555 | 0.4   |            -1.364 |           0.001 |    0.009 |
| Doyle pre-specified  | Insulin prescribed               | No                    | 167 |  1.475 | -0.065 | 3.015 |             1.366 |           0.237 |    1     |
| Doyle pre-specified  | Insulin prescribed               | Yes                   | 154 | -0.032 | -2     | 1.936 |            -0.222 |           0.237 |    1     |
| Doyle pre-specified  | Metformin prescribed             | No                    | 147 |  1.385 | -0.556 | 3.326 |             1.432 |           0.362 |    1     |
| Doyle pre-specified  | Metformin prescribed             | Yes                   | 174 |  0.217 | -1.377 | 1.811 |             0.079 |           0.362 |    1     |
| Added (not in Doyle) | Age                              | <65 y                 | 252 |  0.631 | -0.812 | 2.074 |             0.358 |           0.689 |    1     |
| Added (not in Doyle) | Age                              | >=65 y                |  69 |  1.195 | -1.157 | 3.547 |             1.381 |           0.689 |    1     |
| Added (not in Doyle) | Race/ethnicity                   | Hispanic or non-White |  62 |  1.142 | -1.791 | 4.076 |             0.554 |           0.77  |    1     |
| Added (not in Doyle) | Race/ethnicity                   | Non-Hispanic White    | 259 |  0.659 | -0.71  | 2.027 |             0.617 |           0.77  |    1     |
| Added (not in Doyle) | Baseline A1c test before consent | <=90 d                | 190 |  1.099 | -0.538 | 2.735 |             0.986 |           0.507 |    1     |
| Added (not in Doyle) | Baseline A1c test before consent | >90 d                 | 131 |  0.249 | -1.649 | 2.148 |             0.235 |           0.507 |    1     |

## Honest depth-2 tree (estimation half; leaf-difference p = 0.515)
| level                                                       |   n |   est_discovery |   est |     lo |    hi |
|:------------------------------------------------------------|----:|----------------:|------:|-------:|------:|
| Systolic BP <= 140.5 & Prior-year outpatient visits <= 12.5 | 109 |           0.534 | 1.081 | -1.24  | 3.401 |
| Systolic BP <= 140.5 & Prior-year outpatient visits > 12.5  |  24 |          -4.741 | 1.105 | -3.438 | 5.648 |
| Systolic BP > 140.5                                         |  28 |           4.707 | 3.385 |  0.044 | 6.727 |

## Sensitivity analyses
| analysis                                                  |   n |   ate_aipw |   blp_beta2 |   blp_beta2_se |   blp_p_one_sided |   autoc |   autoc_se |   autoc_p_one_sided |
|:----------------------------------------------------------|----:|-----------:|------------:|---------------:|------------------:|--------:|-----------:|--------------------:|
| Main                                                      | 321 |      0.751 |      -0.906 |          1.227 |             0.781 |   0.249 |      0.673 |               0.356 |
| S1 Outcome = 6-mo level                                   | 321 |     -0.403 |      -0.035 |          1.719 |             0.508 |   0.166 |      1.484 |               0.456 |
| S2 LDL removed                                            | 321 |      0.74  |       0.964 |          1.218 |             0.214 |   0.644 |      0.679 |               0.171 |
| S3 IPW for excluded patients                              | 321 |      0.939 |       0.187 |          1.289 |             0.442 |   0.46  |      0.656 |               0.242 |
| S4 Baseline measured <=180 d before consent               | 310 |      0.881 |      -1.933 |          0.946 |             0.989 |  -1.419 |      0.703 |               0.978 |
| PC injected contrast 1.44 kg = 0.25 SD (older vs younger) | 321 |      0.695 |       0.867 |          1.185 |             0.246 |   0.36  |      0.641 |               0.287 |
| PC injected contrast 2.88 kg = 0.5 SD (older vs younger)  | 321 |      0.702 |       1.093 |          1.079 |             0.155 |   0.461 |      0.661 |               0.242 |
| PC injected contrast 5.76 kg = 1 SD (older vs younger)    | 321 |      0.679 |       2.048 |          0.538 |             0     |   1.729 |      0.65  |               0.004 |

## Root-split stability (200 half-samples)
|                              |   proportion |
|:-----------------------------|-------------:|
| Systolic BP                  |        0.265 |
| Any prior-year ED            |        0.26  |
| Prior-year ED/IP visits      |        0.095 |
| Prior-year outpatient visits |        0.09  |
| LDL                          |        0.07  |
| Age (y)                      |        0.05  |

Figures: fig1-fig6 in this folder. Full console log: console_log.txt.
