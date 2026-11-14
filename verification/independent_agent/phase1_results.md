# Phase 1: blind independent verification

I wrote all of this before opening anything under code/, outputs/, notebooks/ or docs/.
Scripts are in this folder: `v1_replication.py`, `v2_weight_mi.py`, `v2b_weight_mi_12mo.py`, `v3_subgroups.py` (log: `v3_subgroups.log`, CSVs `subgroups_*.csv`), `v4_grf_hba1c.R` (log), and `v4b_grf_crossfit_blp.R` (log).
Python came from the project venv. R grf was installed to `~/Rlibs_verifier`.

| # | Claim | My result | Claimed | Verdict |
|---|---|---|---|---|
| 1a | HbA1c adjusted, 6 mo | FE (within) −0.102 (−0.458 to 0.254), p=.573; FD −0.102 (−0.459 to 0.254), n=349 | −0.10 (−0.46 to 0.25) | VERIFIED |
| 1b | Any ED/IP, 6 mo | FE −7.129 (−15.296 to 1.038); FD −7.129 (−15.301 to 1.043), n=465 | −7.13 (−15.30 to 1.04) | VERIFIED |
| 2 | Weight units | Means/SDs: kg. Unadjusted: **lb** (exact). Adjusted: **kg**, with the gap from complete case explained by the authors' MI (see below) | "adjusted diff is in lb, so ≈ +0.8–0.9 kg and n.s." | **DISCREPANCY (partly)** |
| 3 | Weight × prior ED | No ED +2.46 kg (1.00 to 3.92), n=225; ED −1.19 (−3.06 to 0.68), n=170; interaction p=.0025 (AIPW), .0026 (OLS), .0031 (permutation, 10k); Holm .018 / .018 / .022 | +2.49 (1.05 to 3.92); −1.28 (−3.15 to 0.59); p≈.002; Holm .012 | VERIFIED (small numeric differences) |
| 4 | HbA1c: no HTE | Cross-fitted BLP slope median −4.27 (10 seeds, range −5.4 to −3.4), one-sided p≈.997. grf OOB `test_calibration` slope −14.5, p=1.00. Cross-fitted RATE/AUTOC one-sided p .43–.82. Subgroup min Holm p .87 (OLS), .44 (AIPW) | slope ≈ −4.2, p≈.999; RATE p≈.998; min Holm .88 | VERIFIED (conclusion). RATE p of .998 not reproduced with cross-fitting |
| 5 | ed_i_0mo definition and counts | 191 overall; 94 control, 97 treated (lab sample 148: 81/67) | 191; 94/97 | VERIFIED |

## Claim 1: replication
- I fitted the FE model directly. I stacked the 2 periods, within-transformed by patient, and used the regressors post, post×treat, post×baseline and (for HbA1c) days from consent. SEs are cluster-robust with Stata's xtreg small-sample factor (G/(G−1))·((N−1)/(N−K)) and t on G−1 df.
- I also fitted the first-difference form: Δy ~ 1 + treat + y0 [+ Δtime], with HC1 SEs.
- The two forms give the same point estimate to every digit. CIs differ only in the 3rd decimal, because the df/small-sample factors differ. The algebraic equivalence holds: with T=2, demeaning gives ±Δ/2, and the post dummy becomes the intercept.
- As an extra check I ran outpatient visits (Table 4): 0.995 (0.275 to 1.715), published 1.00 (0.28 to 1.72).

## Claim 2: weight units
- **The data's unit.** `weight_val` is in kg. Implied height √(weight/BMI) has mean 1.68 m (SD 0.10). The data also has `weight_lb_val`, and the ratio is exactly 2.2222 everywhere, so `weight_val = lb × 0.45` (not 0.4536).
- **Means/SDs.** Published means equal lb mean × 0.45359 (108.39, 104.70), while published SDs equal lb SD × 0.45 (29.85, 27.06). Both are kg, but the two conversions are inconsistent (a cosmetic 0.8% difference).
- **Unadjusted difference.** It reproduces exactly in **pounds**: `weight_lb_val_6mo` gives −8.14 (−20.57 to 4.30). The kg version would be −3.66 (−9.26 to 1.93). This part of the claim is correct.
- **Adjusted difference.** Run in the do-file's kg on complete cases (n=395), I get 0.835 (−0.319 to 1.988), p=.16. In lb it is 1.855 (−0.708 to 4.417). The published 1.95 (0.07 to 3.83) matches neither directly.
  - The 5 patients who have 6-mo weight but no baseline weight explain the gap. The do-file imputes their baseline with `mi impute mvn` in long (mlong) format. The predictors are only treat/post/strata/demographics, **not** the patient's own other-period weight, so the imputed baseline is close to a population draw (SD ≈ 28 kg). That puts very noisy Δweight values into the FE regression.
  - I emulated this procedure (`v2_weight_mi.py`: 10 imputations, Rubin's rules, 100 seeds):
    - In kg: median estimate 1.78 (seed 5–95% range 1.63–1.89), CI half-width ≈1.75, p<.05 in 70% of seeds.
    - In lb: 3.97 ± 3.84.
  - 12 mo works the same way. kg-MI gives 1.75 ± 1.96, i.e. about (−0.21 to 3.71), against published 1.88 (−0.14 to 3.89). lb-MI gives 3.91 ± 4.37, which is far off.
- **Conclusion.**
  - The published *adjusted* weight effects match **kg** under the authors' (poor) MI. They are not pounds. Mislabelled units cannot be the explanation for the adjusted column; only the *unadjusted* column is in lb.
  - Still, the substantive point that the kg effect is ≈ +0.8–0.9 kg and not significant holds for a **complete-case** analysis (0.84, p=.16; AIPW 0.89, SE 0.60). The 1.95 kg "significant" result is driven by 5 imputed baselines, not by a unit error.
  - The claim's reasoning (that the adjusted value is in lb) is wrong, even though its numeric conclusion (≈0.85 kg, n.s.) is right for complete cases.

## Claim 3: weight × prior ED
- **Sample.** Both weights observed: N=395, 192 treated. Without ED n=225 (106 treated); with ED n=170 (86 treated). These match the claim.
- **Outcome.** Δweight_val, in the data's kg (×0.45). True kg would be 0.8% larger, which is immaterial.
- **AIPW.** Cross-fitted 5-fold, averaged over 20 repetitions. Ridge-linear outcome models by arm on the 21 baseline covariates plus baseline weight. Propensity is the sample treated share in each of the 4 strata. Subgroup effects are subgroup means of the scores. Results:
  - No ED: +2.46 (1.00 to 3.92).
  - Prior ED: −1.19 (−3.06 to 0.68).
  - Interaction −3.65, p=.0025.
- **Other methods.**
  - OLS Δw ~ T×ED + strata (HC1): interaction −3.61, p=.0026.
  - The same with all baseline covariates: −3.71, p=.0019.
  - Within-stratum re-randomization test: p=.0031.
  - Stratum-weighted (Hajek) subgroup differences: +2.57 and −1.18.
  - Dropping 4 patients with |Δw|>20 kg: −3.27, p=.004.
- **Holm across the 7 variables.** 0.018 (AIPW), 0.018 (OLS), 0.022 (permutation). The claimed 0.012 implies a raw p of about .0017. That is within run-to-run variation of AIPW nuisance fits, so I treat it as VERIFIED.
- **Note.** Under AIPW, insulin's interaction is also nominal (p=.008, Holm .049). It is not significant by OLS (Holm .22) or permutation (Holm .26), so it is method-sensitive.
- **Caveat.** This is a post-hoc subgroup on a secondary outcome. Whether the 7-variable Holm family is the right family is a judgement call (other outcomes were examined as well).

## Claim 4: HbA1c, no heterogeneity
- grf 2.x causal_forest, 4000 trees, W.hat = known stratum propensity, same 21 covariates, median imputation (19 missing cells).
- ATE −0.06 (SE 0.18).
- **Cross-fitted BLP** (5-fold forests, regression-forest m̂, HC3): median slope −4.27, one-sided p≈.997. This matches the claim almost exactly.
- **grf's own `test_calibration` on OOB predictions:** slope −14.5 (t≈−6). A strongly negative calibration slope is the known OOB/leave-out artifact when there is no signal, because CATE predictions are anti-correlated with the patient's own residual. It should be read as "no evidence of heterogeneity", not as evidence of reversed heterogeneity.
- **RATE/AUTOC.** With priorities trained on the other half and evaluated on the held-out half: one-sided p .43–.82 across seeds. With OOB priorities on the same forest: p≈1.00. The claimed p≈.998 is consistent with the OOB/same-sample version or a strongly negative cross-fit draw, but not with my honest split-sample version. The conclusion is the same either way: no detectable heterogeneity.
- **Subgroups.** Min Holm-adjusted interaction p: .87 by OLS (rural site, raw p=.145), matching the claim's .88. By AIPW it is .44 (rural raw p=.063). Neither is significant.

## Claim 5: ed_i_0mo
- README: "Indicator variable for having any emergency department (ED) claims at baseline".
- The paper (Methods) defines baseline use as "the number of episodes in the year prior to consent". Table 1 labels the row "Prior-year health care use: Any ED visit". The do-file uses `ed_i_0mo` for that Table 1 row.
- Counts: 191 (94 control / 97 treated) in the full sample and 148 (81/67) in the lab sample, all matching Table 1.
- Internal consistency: ed_i=1 ⇒ ed_or_ip_i=1 (0 violations), and ed_or_ip_i = max(ed_i, ip_i).
