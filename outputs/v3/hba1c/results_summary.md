# HTE results: food-as-medicine trial, HbA1c change at 6 months (auto-generated)

Run `hba1c` (full): LDL in main model: False;
complete case (no imputation); seed 20260621; 4000 trees; 5-fold cross-fitting
x 5 repeats. Sign: negative = program lowered HbA1c more.

## Bottom line
Formal out-of-fold tests do NOT detect treatment-effect heterogeneity. Subgroup estimates below should be read as consistent with a uniform (null) effect; any nominally significant subgroup is expected under multiple testing.

## Average effect (N=341)
| Estimator | Effect on HbA1c change, pp (95% CI) |
|---|---|
| Paper, adjusted (longitudinal FE; N=349) | -0.10 (-0.46 to +0.25) |
| Paper model replicated (first differences) | -0.09 (-0.45 to +0.27) |
| ANCOVA, prespecified covariates | -0.01 (-0.38 to +0.36) |
| AIPW (cross-fitted, known propensity) | -0.05 (-0.44 to +0.33) |
| Causal forest | -0.03 (-0.53 to +0.46) |

## Omnibus heterogeneity tests (out-of-fold)
| Test | Estimate | p |
|---|---|---|
| BLP calibration beta2 (H1: >0) | -3.48 (SE 1.35) | 0.999 (one-sided) |
| GATES G1 - G3 | +1.55 (SE 0.48) | 0.001 |
| RATE AUTOC (H1: >0) | -0.583 (SE 0.197) | 0.998 (one-sided) |
| Joint BLP on covariates | - | 0.802 |
| ANCOVA joint interactions (HC1 Wald) | - | 0.805 |
| Same, randomization p (constant-effect null) | - | 0.907 |

## Pre-specified subgroups (Doyle family: smallest Holm-adjusted interaction p = 1.000)
| family               | subgroup                         | level                 |   n |    est |     lo |    hi |   paper_model_est |   p_interaction |   p_holm |
|:---------------------|:---------------------------------|:----------------------|----:|-------:|-------:|------:|------------------:|----------------:|---------:|
| Doyle pre-specified  | Baseline HbA1c                   | <9.5%                 | 158 | -0.047 | -0.47  | 0.376 |            -0.076 |           0.978 |    1     |
| Doyle pre-specified  | Baseline HbA1c                   | >=9.5%                | 183 | -0.036 | -0.648 | 0.576 |            -0.146 |           0.978 |    1     |
| Doyle pre-specified  | Site                             | Urban                 | 238 | -0.134 | -0.608 | 0.34  |            -0.172 |           0.446 |    1     |
| Doyle pre-specified  | Site                             | Rural                 | 103 |  0.174 | -0.461 | 0.808 |             0.18  |           0.446 |    1     |
| Doyle pre-specified  | Sex                              | Female                | 182 |  0.207 | -0.332 | 0.747 |             0.156 |           0.17  |    1     |
| Doyle pre-specified  | Sex                              | Male                  | 159 | -0.325 | -0.862 | 0.212 |            -0.36  |           0.17  |    1     |
| Doyle pre-specified  | Prior-year inpatient admission   | No                    | 256 |  0.123 | -0.309 | 0.556 |             0.075 |           0.157 |    1     |
| Doyle pre-specified  | Prior-year inpatient admission   | Yes                   |  85 | -0.536 | -1.339 | 0.267 |            -0.622 |           0.157 |    1     |
| Doyle pre-specified  | Prior-year ED visit              | No                    | 194 |  0.106 | -0.366 | 0.577 |             0.126 |           0.398 |    1     |
| Doyle pre-specified  | Prior-year ED visit              | Yes                   | 147 | -0.234 | -0.866 | 0.397 |            -0.266 |           0.398 |    1     |
| Doyle pre-specified  | Insulin prescribed               | No                    | 168 |  0.046 | -0.488 | 0.579 |             0.051 |           0.66  |    1     |
| Doyle pre-specified  | Insulin prescribed               | Yes                   | 173 | -0.125 | -0.673 | 0.422 |            -0.219 |           0.66  |    1     |
| Doyle pre-specified  | Metformin prescribed             | No                    | 168 |  0.018 | -0.501 | 0.538 |            -0.058 |           0.764 |    1     |
| Doyle pre-specified  | Metformin prescribed             | Yes                   | 173 | -0.099 | -0.659 | 0.462 |            -0.236 |           0.764 |    1     |
| Added (not in Doyle) | Age                              | <65 y                 | 265 |  0.083 | -0.37  | 0.537 |            -0.072 |           0.172 |    0.516 |
| Added (not in Doyle) | Age                              | >=65 y                |  76 | -0.475 | -1.135 | 0.185 |            -0.067 |           0.172 |    0.516 |
| Added (not in Doyle) | Race/ethnicity                   | Hispanic or non-White |  61 |  0.535 | -0.485 | 1.554 |             0.139 |           0.211 |    0.516 |
| Added (not in Doyle) | Race/ethnicity                   | Non-Hispanic White    | 280 | -0.166 | -0.575 | 0.242 |            -0.205 |           0.211 |    0.516 |
| Added (not in Doyle) | Baseline A1c test before consent | <=90 d                | 196 |  0.003 | -0.524 | 0.53  |            -0.128 |           0.789 |    0.789 |
| Added (not in Doyle) | Baseline A1c test before consent | >90 d                 | 145 | -0.101 | -0.651 | 0.449 |            -0.132 |           0.789 |    0.789 |

## Honest depth-2 tree (estimation half; leaf-difference p = 0.181)
| level                          |   n |   est_discovery |    est |     lo |    hi |
|:-------------------------------|----:|----------------:|-------:|-------:|------:|
| Age (y) <= 48                  |  43 |           1.273 | -0.722 | -1.998 | 0.554 |
| Age (y) > 48 & Age (y) <= 53.5 |  29 |          -1.645 |  0.861 | -0.304 | 2.026 |
| Age (y) > 48 & Age (y) > 53.5  |  99 |          -0.05  | -0.102 | -0.801 | 0.598 |

## Sensitivity analyses
| analysis                                                  |   n |   ate_aipw |   blp_beta2 |   blp_beta2_se |   blp_p_one_sided |   autoc |   autoc_se |   autoc_p_one_sided |
|:----------------------------------------------------------|----:|-----------:|------------:|---------------:|------------------:|--------:|-----------:|--------------------:|
| Main                                                      | 341 |     -0.053 |      -3.483 |          1.351 |             0.999 |  -0.583 |      0.197 |               0.998 |
| S1 Outcome = 6-mo level                                   | 341 |     -0.046 |      -4.684 |          1.771 |             0.998 |  -0.5   |      0.175 |               0.998 |
| S3 IPW for excluded patients                              | 341 |      0.03  |      -3.876 |          1.763 |             0.994 |  -0.496 |      0.19  |               0.996 |
| S4 Baseline measured <=180 d before consent               | 276 |      0.033 |      -2.855 |          1.906 |             0.933 |  -0.17  |      0.233 |               0.768 |
| PC injected contrast 0.57 pp = 0.25 SD (older vs younger) | 341 |     -0.021 |      -1.191 |          1.158 |             0.848 |  -0.218 |      0.171 |               0.899 |
| PC injected contrast 1.13 pp = 0.5 SD (older vs younger)  | 341 |     -0.044 |       1.109 |          0.73  |             0.064 |   0.171 |      0.154 |               0.134 |
| PC injected contrast 2.26 pp = 1 SD (older vs younger)    | 341 |     -0.084 |       1.506 |          0.289 |             0     |   0.698 |      0.174 |               0     |

## Root-split stability (200 half-samples)
|                              |   proportion |
|:-----------------------------|-------------:|
| Age (y)                      |        0.275 |
| Baseline HbA1c (%)           |        0.155 |
| Systolic BP                  |        0.11  |
| BMI                          |        0.075 |
| Any prior-year inpatient     |        0.075 |
| Prior-year outpatient visits |        0.07  |

Figures: fig1-fig6 in this folder. Full console log: console_log.txt.
