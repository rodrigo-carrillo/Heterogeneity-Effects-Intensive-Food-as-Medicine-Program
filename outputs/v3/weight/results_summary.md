# HTE results: food-as-medicine trial, weight change at 6 months (auto-generated)

Run `weight` (full): LDL in main model: False;
complete case (no imputation); seed 20260621; 4000 trees; 5-fold cross-fitting
x 5 repeats. Sign: negative = less weight gain with program.

## Bottom line
Formal out-of-fold tests do NOT detect treatment-effect heterogeneity. Subgroup estimates below should be read as consistent with a uniform (null) effect; any nominally significant subgroup is expected under multiple testing.

## Average effect (N=392)
| Estimator | Effect on weight change, kg (95% CI) |
|---|---|
| Paper, adjusted (longitudinal FE; N=400; unadjusted published in lb; adjusted kg with 5 imputed baselines) | +1.95 (+0.07 to +3.83) |
| Paper model replicated (first differences) | +0.82 (-0.35 to +1.99) |
| ANCOVA, prespecified covariates | +0.70 (-0.52 to +1.92) |
| AIPW (cross-fitted, known propensity) | +0.90 (-0.28 to +2.09) |
| Causal forest | +0.88 (-0.83 to +2.59) |

## Omnibus heterogeneity tests (out-of-fold)
| Test | Estimate | p |
|---|---|---|
| BLP calibration beta2 (H1: >0) | -0.74 (SE 1.03) | 0.797 (one-sided) |
| GATES G1 - G3 | -0.55 (SE 1.40) | 0.695 |
| RATE AUTOC (H1: >0) | -0.128 (SE 0.577) | 0.588 (one-sided) |
| Joint BLP on covariates | - | 0.092 |
| ANCOVA joint interactions (HC1 Wald) | - | 0.123 |
| Same, randomization p (constant-effect null) | - | 0.286 |

## Pre-specified subgroups (Doyle family: smallest Holm-adjusted interaction p = 0.013)
| family               | subgroup                         | level                 |   n |    est |     lo |    hi |   paper_model_est |   p_interaction |   p_holm |
|:---------------------|:---------------------------------|:----------------------|----:|-------:|-------:|------:|------------------:|----------------:|---------:|
| Doyle pre-specified  | Baseline HbA1c                   | <9.5%                 | 172 |  1.831 |  0.108 | 3.553 |             1.343 |           0.167 |    0.67  |
| Doyle pre-specified  | Baseline HbA1c                   | >=9.5%                | 220 |  0.169 | -1.443 | 1.781 |             0.203 |           0.167 |    0.67  |
| Doyle pre-specified  | Site                             | Urban                 | 273 |  1.176 | -0.195 | 2.548 |             0.921 |           0.501 |    1     |
| Doyle pre-specified  | Site                             | Rural                 | 119 |  0.26  | -2.03  | 2.55  |             0.526 |           0.501 |    1     |
| Doyle pre-specified  | Sex                              | Female                | 211 |  1.772 |  0.079 | 3.464 |             1.784 |           0.113 |    0.566 |
| Doyle pre-specified  | Sex                              | Male                  | 181 | -0.12  | -1.738 | 1.498 |            -0.259 |           0.113 |    0.566 |
| Doyle pre-specified  | Prior-year inpatient admission   | No                    | 296 |  0.952 | -0.375 | 2.28  |             0.97  |           0.88  |    1     |
| Doyle pre-specified  | Prior-year inpatient admission   | Yes                   |  96 |  0.731 | -1.825 | 3.287 |             0.354 |           0.88  |    1     |
| Doyle pre-specified  | Prior-year ED visit              | No                    | 222 |  2.54  |  1.069 | 4.011 |             2.481 |           0.002 |    0.013 |
| Doyle pre-specified  | Prior-year ED visit              | Yes                   | 170 | -1.246 | -3.131 | 0.638 |            -1.239 |           0.002 |    0.013 |
| Doyle pre-specified  | Insulin prescribed               | No                    | 197 |  2.28  |  0.726 | 3.833 |             1.917 |           0.02  |    0.123 |
| Doyle pre-specified  | Insulin prescribed               | Yes                   | 195 | -0.498 | -2.259 | 1.264 |            -0.344 |           0.02  |    0.123 |
| Doyle pre-specified  | Metformin prescribed             | No                    | 189 |  1.532 | -0.227 | 3.29  |             1.601 |           0.311 |    0.933 |
| Doyle pre-specified  | Metformin prescribed             | Yes                   | 203 |  0.308 | -1.277 | 1.894 |             0.126 |           0.311 |    0.933 |
| Added (not in Doyle) | Age                              | <65 y                 | 307 |  0.671 | -0.686 | 2.029 |             0.528 |           0.453 |    0.952 |
| Added (not in Doyle) | Age                              | >=65 y                |  85 |  1.717 | -0.652 | 4.085 |             1.962 |           0.453 |    0.952 |
| Added (not in Doyle) | Race/ethnicity                   | Hispanic or non-White |  71 |  1.849 | -0.746 | 4.445 |             1.155 |           0.435 |    0.952 |
| Added (not in Doyle) | Race/ethnicity                   | Non-Hispanic White    | 321 |  0.688 | -0.635 | 2.011 |             0.726 |           0.435 |    0.952 |
| Added (not in Doyle) | Baseline A1c test before consent | <=90 d                | 226 |  1.416 | -0.123 | 2.955 |             1.388 |           0.317 |    0.952 |
| Added (not in Doyle) | Baseline A1c test before consent | >90 d                 | 166 |  0.193 | -1.645 | 2.031 |             0.072 |           0.317 |    0.952 |

## Honest depth-2 tree (estimation half; leaf-difference p = 0.264)
| level                                                                            |   n |   est_discovery |    est |     lo |    hi |
|:---------------------------------------------------------------------------------|----:|----------------:|-------:|-------:|------:|
| Prior-year outpatient visits <= 3.5                                              |  58 |          -1.632 |  2.767 | -0.883 | 6.418 |
| Prior-year outpatient visits > 3.5 & Days, baseline A1c test to consent <= 100.5 |  86 |           0.595 |  0.524 | -2.416 | 3.464 |
| Prior-year outpatient visits > 3.5 & Days, baseline A1c test to consent > 100.5  |  52 |           4.323 | -1.063 | -3.857 | 1.731 |

## Sensitivity analyses
| analysis                                                  |   n |   ate_aipw |   blp_beta2 |   blp_beta2_se |   blp_p_one_sided |   autoc |   autoc_se |   autoc_p_one_sided |
|:----------------------------------------------------------|----:|-----------:|------------:|---------------:|------------------:|--------:|-----------:|--------------------:|
| Main                                                      | 392 |      0.902 |      -0.74  |          1.031 |             0.797 |  -0.128 |      0.577 |               0.588 |
| S1 Outcome = 6-mo level                                   | 392 |      0.378 |       0.196 |          1.475 |             0.447 |   0.037 |      0.891 |               0.483 |
| S3 IPW for excluded patients                              | 392 |      1.074 |       0.883 |          1.266 |             0.243 |   0.575 |      0.738 |               0.218 |
| S4 Baseline measured <=180 d before consent               | 375 |      1.088 |      -0.609 |          1.072 |             0.725 |  -0.289 |      0.61  |               0.682 |
| PC injected contrast 1.50 kg = 0.25 SD (older vs younger) | 392 |      0.902 |       1.483 |          1.318 |             0.134 |   0.665 |      0.696 |               0.17  |
| PC injected contrast 3.01 kg = 0.5 SD (older vs younger)  | 392 |      0.885 |       1.838 |          1.179 |             0.06  |   0.575 |      0.7   |               0.206 |
| PC injected contrast 6.02 kg = 1 SD (older vs younger)    | 392 |      0.775 |       2.336 |          0.549 |             0     |   1.428 |      0.68  |               0.018 |

## Root-split stability (200 half-samples)
|                              |   proportion |
|:-----------------------------|-------------:|
| Systolic BP                  |        0.27  |
| Any prior-year ED            |        0.27  |
| Prior-year ED/IP visits      |        0.11  |
| Insulin                      |        0.085 |
| Baseline HbA1c (%)           |        0.065 |
| Prior-year outpatient visits |        0.05  |

Figures: fig1-fig6 in this folder. Full console log: console_log.txt.
