# HTE results: food-as-medicine trial, any ED/IP use at 6 months (auto-generated)

Run `ed_ip_ldl_cc` (full): LDL in main model: True;
complete case (no imputation); seed 20260621; 4000 trees; 5-fold cross-fitting
x 5 repeats. Sign: negative = fewer patients with ED/IP use.

## Bottom line
Formal out-of-fold tests do NOT detect treatment-effect heterogeneity. Subgroup estimates below should be read as consistent with a uniform (null) effect; any nominally significant subgroup is expected under multiple testing.

## Average effect (N=366)
| Estimator | Effect on any ED/IP use, pp (95% CI) |
|---|---|
| Paper, adjusted (longitudinal FE; N=465) | -7.13 (-15.30 to +1.04) |
| Paper model replicated (first differences) | -4.60 (-13.86 to +4.67) |
| ANCOVA, prespecified covariates | -5.43 (-14.98 to +4.13) |
| AIPW (cross-fitted, known propensity) | -2.41 (-11.71 to +6.90) |
| Causal forest | -2.44 (-14.40 to +9.52) |

## Omnibus heterogeneity tests (out-of-fold)
| Test | Estimate | p |
|---|---|---|
| BLP calibration beta2 (H1: >0) | -0.65 (SE 1.05) | 0.735 (one-sided) |
| GATES G1 - G3 | +7.11 (SE 11.61) | 0.541 |
| RATE AUTOC (H1: >0) | -4.436 (SE 4.615) | 0.832 (one-sided) |
| Joint BLP on covariates | - | 0.000 |
| ANCOVA joint interactions (HC1 Wald) | - | 0.004 |
| Same, randomization p (constant-effect null) | - | 0.046 |

## Pre-specified subgroups (Doyle family: smallest Holm-adjusted interaction p = 0.470)
| family               | subgroup                         | level                 |   n |     est |      lo |     hi |   paper_model_est |   p_interaction |   p_holm |
|:---------------------|:---------------------------------|:----------------------|----:|--------:|--------:|-------:|------------------:|----------------:|---------:|
| Doyle pre-specified  | Baseline HbA1c                   | <9.5%                 | 171 |   2.849 | -11.064 | 16.763 |             1.272 |           0.321 |    1     |
| Doyle pre-specified  | Baseline HbA1c                   | >=9.5%                | 195 |  -6.591 | -18.982 |  5.801 |            -9.732 |           0.321 |    1     |
| Doyle pre-specified  | Site                             | Urban                 | 261 |  -0.864 | -11.721 |  9.992 |            -2.004 |           0.666 |    1     |
| Doyle pre-specified  | Site                             | Rural                 | 105 |  -5.451 | -23.259 | 12.358 |           -11.042 |           0.666 |    1     |
| Doyle pre-specified  | Sex                              | Female                | 207 |  -3.578 | -16.02  |  8.864 |            -6.377 |           0.735 |    1     |
| Doyle pre-specified  | Sex                              | Male                  | 159 |  -0.36  | -14.272 | 13.552 |            -1.861 |           0.735 |    1     |
| Doyle pre-specified  | Prior-year inpatient admission   | No                    | 287 |  -4.88  | -15.102 |  5.343 |            -6.22  |           0.303 |    1     |
| Doyle pre-specified  | Prior-year inpatient admission   | Yes                   |  79 |   7.627 | -13.888 | 29.143 |             1.03  |           0.303 |    1     |
| Doyle pre-specified  | Prior-year ED visit              | No                    | 217 |  -0.398 | -11.504 | 10.707 |            -1.131 |           0.66  |    1     |
| Doyle pre-specified  | Prior-year ED visit              | Yes                   | 149 |  -4.775 | -20.823 | 11.273 |            -7.21  |           0.66  |    1     |
| Doyle pre-specified  | Insulin prescribed               | No                    | 190 |   6.107 |  -6.99  | 19.204 |             5.396 |           0.067 |    0.47  |
| Doyle pre-specified  | Insulin prescribed               | Yes                   | 176 | -11.127 | -24.124 |  1.871 |           -15.46  |           0.067 |    0.47  |
| Doyle pre-specified  | Metformin prescribed             | No                    | 175 | -10.037 | -23.154 |  3.081 |           -12.992 |           0.11  |    0.662 |
| Doyle pre-specified  | Metformin prescribed             | Yes                   | 191 |   5.018 |  -7.998 | 18.034 |             4.273 |           0.11  |    0.662 |
| Added (not in Doyle) | Age                              | <65 y                 | 288 |  -0.319 | -10.509 |  9.872 |            -2.572 |           0.478 |    0.5   |
| Added (not in Doyle) | Age                              | >=65 y                |  78 |  -9.053 | -30.896 | 12.789 |           -12.404 |           0.478 |    0.5   |
| Added (not in Doyle) | Race/ethnicity                   | Hispanic or non-White |  75 |  18.381 |  -1.679 | 38.442 |            18.295 |           0.025 |    0.074 |
| Added (not in Doyle) | Race/ethnicity                   | Non-Hispanic White    | 291 |  -7.479 | -17.853 |  2.894 |           -10.357 |           0.025 |    0.074 |
| Added (not in Doyle) | Baseline A1c test before consent | <=90 d                | 215 |  -6.734 | -18.818 |  5.351 |            -9.495 |           0.25  |    0.5   |
| Added (not in Doyle) | Baseline A1c test before consent | >90 d                 | 151 |   4.303 | -10.112 | 18.719 |             2.098 |           0.25  |    0.5   |

## Honest depth-2 tree (estimation half; leaf-difference p = 0.204)
| level                                     |   n |   est_discovery |     est |      lo |     hi |
|:------------------------------------------|----:|----------------:|--------:|--------:|-------:|
| Systolic BP <= 121 & Diastolic BP <= 65.5 |  22 |          -3.461 |  -9.268 | -44.363 | 25.828 |
| Systolic BP <= 121 & Diastolic BP > 65.5  |  48 |          37.236 |  10.11  | -15.783 | 36.002 |
| Systolic BP > 121 & LDL <= 99.5           |  61 |         -42.222 |   2.972 | -20.829 | 26.774 |
| Systolic BP > 121 & LDL > 99.5            |  52 |           3.766 | -28.841 | -56.458 | -1.225 |

## Sensitivity analyses
| analysis                                                   |   n |   ate_aipw |   blp_beta2 |   blp_beta2_se |   blp_p_one_sided |   autoc |   autoc_se |   autoc_p_one_sided |
|:-----------------------------------------------------------|----:|-----------:|------------:|---------------:|------------------:|--------:|-----------:|--------------------:|
| Main                                                       | 366 |     -2.408 |      -0.65  |          1.047 |             0.735 |  -4.436 |      4.615 |               0.832 |
| S1 Outcome = no. of ED/IP visits                           | 366 |     -0.032 |      -2.005 |          1.551 |             0.87  |  -0.147 |      0.157 |               0.825 |
| S2 LDL removed                                             | 366 |     -1.945 |      -0.048 |          1.146 |             0.517 |  -3.189 |      4.634 |               0.754 |
| S3 IPW for excluded patients                               | 366 |     -2.298 |      -0.222 |          1.247 |             0.571 |  -5.67  |      4.492 |               0.897 |
| PC injected contrast 11.89 pp = 0.25 SD (older vs younger) | 366 |     -1.715 |      -0.103 |          1.17  |             0.535 |  -4.942 |      4.539 |               0.862 |
| PC injected contrast 23.79 pp = 0.5 SD (older vs younger)  | 366 |     -1.717 |       0.718 |          1.016 |             0.24  |  -2.131 |      4.787 |               0.672 |
| PC injected contrast 47.58 pp = 1 SD (older vs younger)    | 366 |     -2.053 |       1.871 |          0.502 |             0     |   9.153 |      5.239 |               0.04  |

## Root-split stability (200 half-samples)
|                                    |   proportion |
|:-----------------------------------|-------------:|
| Systolic BP                        |        0.485 |
| Baseline HbA1c (%)                 |        0.12  |
| Age (y)                            |        0.085 |
| Days, baseline A1c test to consent |        0.065 |
| Prior-year outpatient visits       |        0.06  |
| Diastolic BP                       |        0.04  |

Figures: fig1-fig6 in this folder. Full console log: console_log.txt.
