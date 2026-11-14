# Agent E: red team (29 Oct 2026)

Saved by the main session from Agent E's final message (the agent's tool was blocked from writing
this file). Content is the agent's, lightly condensed. Not yet independently checked by the main
session.

**Overall verdict:** no coding error, leakage, post-randomization covariate or silent patient drop
found. The primary HbA1c null stands, but only for large heterogeneity. Weight × prior ED is a
real pattern in the data but weaker than the letter presents it: nominal, post hoc,
hypothesis-generating.

**Incomplete:** two background jobs were stopped by Claude Code for low system memory and not
restarted: `e03_placebo_planted.py` (71 of 310 tasks, ~7 per scenario; indicative only) and
`e04_subgroup_fwer.py` (no output). `e05_alt_learners.py` written, not run.

## Findings

| ID | Severity | Finding | Evidence | Effect | Fix |
|---|---|---|---|---|---|
| E1 | MAJOR | Letter implies weight × prior ED was prespecified and contrasts "prespecified decision rules" with Alsan's post hoc analysis. Doyle prespecified subgroups for engagement, HbA1c and utilization, not weight (docs/04 §4); the SAP was drafted 21 Oct in the session in which weight × ED first appeared and approved 28 Oct, when it was named the main finding. | letter_text.py Results ¶2, Discussion ¶1; WORKFLOW "Answered (2026-10-28)" | Credibility framing overstated | Subgroup *variables* prespecified (for HbA1c/utilization); applying them to weight and headlining it was post hoc; soften the Alsan contrast |
| E2 | MAJOR | Weight × prior ED fragile; such hits common by chance in the search run | Null (within-stratum permuted) runs: 2/8 gave Holm P < .05 among Doyle subgroups in ≥1 of 4 outcomes, 1/8 ≤ .013 (main runs only; actual search was wider). Removing 3 most supportive patients → raw P above Holm boundary (.0071); 8 of 392 → above .05. Excluding 39 with \|Δweight\| > 10 kg → P .10 (−1.49 vs −3.68). Most influential patients changed ±18–22 kg (unchecked). Treatment × outpatient visits P .013; × insulin .0075; fully interacted .031. Prior-ED patients: 2× outpatient visits, insulin 59% vs 42%, inpatient 38% vs 14%. Planted effect of observed size (−3.8 kg) Holm-significant in 2/7 (mean −3.2): estimate likely inflated | "1 subgroup did [modify]" reads as established | Report as nominal (P .002, Holm .01, .05 across outcomes, null forest tests); disclose search breadth and sensitivity to extreme changes |
| E3 | MAJOR | "Did not appear to hide benefit" overstates power | Positive control: HbA1c 1 SD (2.26 pp) 7/7, 0.5 SD 3/7; weight 1 SD 5/7, 0.5 SD 0/7; linear modifier with effect SD 0.5 pp 2/7 | Only large heterogeneity (≈2 pp+) ruled out | "No large heterogeneity; modest heterogeneity not excluded"; 2.3 not 2.2 |
| E4 | MINOR | Omnibus tests lean negative under the null (conservative) | 29 null permutations: median calibration z −1.0, RATE z −1.2 (range ≈ −3.1 to +1.6); 0/29 false positives; real HbA1c z −2.6/−3.0 inside null range | No over-rejection; lower power | Drop "both P > .99" or note the shift |
| E5 | MINOR | HbA1c complete-case inclusion differs by arm × prior ED (ED: control 86% vs intervention 68%; no ED 67% vs 74%; interaction P .002; joint P .008). S3 reweighting is additive, cannot correct; docs/02 calls attrition "balanced" (true overall only) | e06, e07 | Arm-specific IPW leaves HbA1c subgroups essentially unchanged (ED −0.63 → −0.65); weight inclusion not differential (P .61) | Limitation; optionally arm-specific S3 |
| E6 | MINOR | Mechanism wording vs cell means: 6-mo change no-ED control −1.98 kg, no-ED intervention +0.64, ED control −0.32, ED intervention −1.51; ~half of +2.54 is loss among no-ED controls | e02 | Interpretive | Present relative to controls; drop "ate more … gained" |
| E7 | MINOR | N = 392 OLS (.0025) and permutation (.0028) values come from no repository script; independent agent's script used N = 395, old kg factor, sharp-null permutation | grep; verification/independent_agent/v3_subgroups.py | Small | Add a script; constant-effect null |
| E8 | MINOR | Honest single splits choose top-5 variables from a forest fit on all patients | hte_food_as_medicine.py:610, 843 | Mild optimism; all null | Importance from discovery half only |
| E9 | NOTE | Omnibus P = median across repeats without ×2 correction (logged in docs/05) | hte:570, 588 | None while all null | Apply ×2 if a positive is claimed |
| E10 | NOTE | "6-month" outcomes measured 150–296 d after consent (control median 183, p75 202): many controls measured after the wait list ended | e01 | Dilutes contrasts; weight × ED stronger restricted to ≤190 d (P .0018) | Optional limitation |
| E11 | NOTE | 1 SD HbA1c 2.26 not ≈2.2; weight has 22 modifiers not 21; ×4 Holm P .054 | summary_across_outcomes.md | None | Fix wording |

## Robustness of weight × prior ED (e02, N = 392; base −3.68 kg, P .0025)
Holds: covariate adjustment (baseline weight, BMI, all covariates, treatment × weight/BMI, timing)
P .0015–.003; 6-month level .002; Huber .003; median regression .006; winsorized .003–.004;
leave-one-out max .004; same sign in every site, HbA1c stratum and sex; within-ED balance
|SMD| ≤ 0.30. Weakens: excluding \|Δ\| > 10 kg .10; treatment × outpatient visits .013; fully
interacted .031; ≥2 ED/IP visits .05; continuous ED/IP count .037; any inpatient stay .78.

## Code checks passed
Outcomes (6-mo minus baseline, sign, lb/2.20462, ED/IP = 100 × indicator); all baseline dates
before consent; complete-case filter drops nothing else; treatment coding; propensity = treated
fraction within randomization stratum; AIPW formula and SE; CATEs and AIPW nuisances out-of-fold on
the same folds (no leakage); calibration test 1-sided; standard RATE bootstrap; interaction test
re-randomizes within strata under a constant-effect null; Holm within family; fixed seeds.
Positive control recovered unbiased (−2.29 vs 2.26 pp; −6.41 vs 6.02 kg); planted ED effect
recovered with correct sign. HbA1c CATEs heavily shrunk (SD 0.09 pp averaged, ~0.16 per repeat)
but not collapsed: the null is not a bug.

Files: e01_data_checks.txt, e02_weight_ed_robustness.txt, e02_weight_ed_specs.csv,
e02_balance_within_ed.csv, e03_results.jsonl, e03_summary_partial.txt, e06_missingness.{txt,csv},
e07_armspecific_ipw.{txt,csv}, scripts/.
