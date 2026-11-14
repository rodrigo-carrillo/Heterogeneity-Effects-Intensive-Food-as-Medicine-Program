# Detailed methods

Methods details and supporting tables for the manuscript *Heterogeneity in the Effects of a
Food-as-Medicine Program: Secondary Analysis of a Randomized Trial*. The article format allows no
supplementary files, so this material is kept here. Section numbers below are referred to in this
file only; "the Table" and "the Figure" are those of the article.

## 1. Data, outcomes, and candidate modifiers

**Source:** replication data posted by the trial investigators on Harvard Dataverse
(doi:[10.7910/DVN/LZUTU0](https://doi.org/10.7910/DVN/LZUTU0), version 2.0), the 465 patients who
completed the trial (of 500 randomized), with electronic health record data. One baseline
triglyceride value recorded as the placeholder 9999999 was set to missing.

**Outcomes** (six months after randomization): change in HbA1c from baseline (percentage points;
primary); change in weight (kg; computed from the posted weight in pounds divided by 2.20462); change
in systolic blood pressure (mm Hg); any emergency department visit or inpatient admission
(percentage points).

**Candidate modifiers** (21; measured before randomization; selected by me for this analysis from the
baseline variables in the public data): age; sex; Hispanic ethnicity; non-Hispanic Black race; rural
site; baseline HbA1c and HbA1c stratum (≥9.5%); days between the baseline HbA1c test and consent; body
mass index; systolic and diastolic blood pressure; prescriptions for insulin, metformin,
glucagon-like peptide 1 receptor agonists, sulfonylureas, other glucose-lowering drugs, and any
glucose-lowering drug; prior-year outpatient visits; prior-year emergency department or inpatient
visits; any prior-year inpatient admission; any prior-year emergency department visit. Baseline
weight was added for the weight outcome. Variables measured after randomization were excluded.

**Complete case:** a patient was analyzed for an outcome if the outcome, its baseline value, and every
modifier were observed (HbA1c, n = 341; weight, n = 392; systolic blood pressure, n = 371; emergency
department visit or hospitalization, n = 456). A sensitivity analysis also required baseline
low-density lipoprotein cholesterol and added it as a modifier; another reweighted the complete cases
to all 465 patients by inverse probability of inclusion. In these sensitivity analyses, no outcome met
the criterion for heterogeneity. In the analysis with low-density lipoprotein cholesterol, the
interaction test (section 2) had P < 0.05 for systolic blood pressure (P = 0.02) and for emergency
department visit or hospitalization (P = 0.046). That is, in these smaller samples (308 and 366
patients), the effects on these two outcomes appeared to vary with baseline characteristics. The test
does not identify which characteristics; low-density lipoprotein cholesterol was not singled out (in
a linear model of the patient-level augmented scores (section 2) on all modifiers, its coefficient had
P = 0.40 for systolic blood pressure and P = 0.36 for emergency department visit or hospitalization,
and no modifier had a Holm-adjusted P < 0.05). This does not change the findings or the analysis plan:
the interaction test is a supporting check, not part of the decision rule (section 2); in the same
analyses, the calibration and rank-weighted tests were null (all P ≥ 0.73) and none of the seven
prespecified subgroup variables had a Holm-adjusted P < 0.05; and in the main analyses of the same
outcomes, which included more patients, the interaction test had P = 0.08 and P = 0.21 (Table in the
article). Under the analysis plan, a signal in a secondary outcome is considered further only if it
is consistent across both covariate sets (with and without low-density lipoprotein cholesterol);
these were not.

Inclusion in the HbA1c sample differed by arm according to prior-year emergency department use. Among
patients with a prior-year visit, 81 of 94 (86%) in the wait-list group and 66 of 97 (68%) in the
intervention group were included; among patients without one, 95 of 141 (67%) and 99 of 133 (74%)
(interaction between arm and prior visit in a logistic model of inclusion, P = 0.002). This is noted
as a limitation.

## 2. Estimation and tests

**Average effects:** augmented inverse probability weighting with cross-fitted random-forest outcome
models; the probability of assignment to the intervention was the proportion assigned within each
randomization stratum (site x HbA1c stratum) of the analyzed sample. In plain terms, inverse
probability weighting compares the groups after weighting each patient by 1 divided by the
probability of the assignment that the patient received; in a randomized trial, this is a difference
in means within randomization strata. The augmented version adds, for each patient, the outcome
predicted from baseline characteristics by a model fitted in each group, and weights only the part of
the outcome that this model does not predict. The estimate remains valid if either the assignment
probabilities or the outcome models are correct (here the assignment probabilities are known from the
design), and the outcome models make it more precise when baseline characteristics predict the
outcome. The same patient-level augmented scores were used for the subgroup effects and the
calibration test.

**Individual effects:** causal forest (EconML 0.17, CausalForestDML; 4000 honest trees, minimum leaf
size 10, subsample fraction 0.45). Estimates were out of fold (cross-fitting). The patients were split
at random into five groups (folds) of similar size and similar composition by randomization stratum
and arm; a forest was trained on four folds (about 80% of the patients) and used to estimate the
effects of the patients in the remaining fold, and this was done for each fold in turn. Every
patient's estimate therefore came from a forest that had not seen that patient's data. The procedure
was repeated five times with different random splits (25 forests in total), giving five out-of-fold
estimates per patient. The tests were computed within each repeat, and their results were combined
across repeats by the median (the median P for the calibration test; the median estimate and SE for
the rank-weighted test).

**Prespecified decision rule:** (1) calibration test, the best linear projection of the doubly robust
score on the centered out-of-fold estimate (1-sided test that the slope is positive); (2) rank-weighted
average treatment effect, area under the targeting operator characteristic curve, with bootstrap SE
(1-sided). Heterogeneity was considered detected only if both P < 0.05. As a parametric check, all
treatment x covariate interactions were tested jointly in a linear regression of the outcome on
treatment, every modifier, and every treatment x modifier product (a Wald test, with robust SEs, that
all products are zero). With about 20 products and a few hundred patients, the usual large-sample P
value of this test can be too small, so its P value was obtained by randomization: treatment labels
were reassigned at random 2000 times within the randomization strata, as in the trial, with each
patient's outcome set to what it would have been under a common treatment effect equal to the
estimated average effect; the P value is the proportion of reassignments with a test statistic at
least as large as the observed one.

**Prespecified subgroup variables:** the seven subgroup variables the trial investigators prespecified
for HbA1c and health care use that are available in the public data (baseline HbA1c <9.5% vs ≥9.5%,
site, sex, any prior-year inpatient admission, any prior-year emergency department visit, insulin
prescribed, metformin prescribed). Within each level, effects were estimated by augmented inverse
probability weighting; the interaction P was Holm-adjusted across the seven variables. Applying them
to weight and systolic blood pressure was my decision, not the trial investigators' (who prespecified
subgroups only for HbA1c and health care use), and is therefore post hoc. Post hoc, after the weight x
prior-year emergency department difference was seen: Holm adjustment across the four outcomes;
ordinary least squares with stratum indicators and robust SEs (HC1, the type that Stata uses with its
robust option); a within-stratum permutation test (10 000 permutations, constant-effect null); and
ordinary least squares after excluding patients whose weight changed by more than 10 kg; after this
exclusion, the interaction was no longer significant (difference in effect, 1.5 kg vs 3.7 kg in the
same model without the exclusion; P = 0.10 vs P = 0.002, neither adjusted for multiplicity).

**Detectable heterogeneity** (positive control; a check prespecified in the analysis plan, similar in
aim to a power calculation but, unlike post hoc power, not based on the observed effects): an
artificial difference in effect between older and younger patients of 0.25, 0.5, or 1 SD of the
outcome was added to the treated patients' outcomes, and the calibration and rank-weighted tests were
rerun, with all candidate modifiers and three cross-fitting repeats (one run per size); only the 1 SD
difference was detected, for every outcome. Age is not an outcome here: it only defines the two halves
of the sample that received different added effects (in the intervention group, half the difference
was subtracted from the outcomes of the older half and added to those of the younger half, so the
average effect was nearly unchanged). The outcomes and the tests were those of the main analysis; any
baseline characteristic could have been used to define the halves.

*Purpose:* the true individual effects are unknown, so a difference of known size was added to the real
data to see whether the calibration and rank-weighted tests would find it. This shows how large
heterogeneity must be for this trial's sample sizes and tests to detect it and confirms that the
analysis detects heterogeneity when it is large enough.

*Interpretation:* the calibration and rank-weighted tests could detect a difference in effect between
two halves of the sample (patients at or above vs below the median age) of 1 SD of the outcome (HbA1c,
2.26 percentage points; weight, 6.02 kg; systolic blood pressure, 20.9 mm Hg; emergency department
visit or hospitalization, 47 percentage points), but not of 0.5 SD (for weight, 3.01 kg). The null
results therefore exclude only large heterogeneity; smaller differences, which may still matter
clinically, cannot be excluded in this study: with this sample size and with all candidate modifiers
combined, these tests have too little power to detect them. A test of one variable at a time can
detect smaller differences along that variable. For weight by prior-year emergency department visit,
the SE of the difference in effect between the two subgroups was 1.22 kg, so this test had 80% power
to detect a difference of 3.4 kg (0.57 SD; 2-sided P < 0.05), or 4.3 kg (0.72 SD) after Holm
adjustment for the seven subgroup variables. Across all outcomes and subgroup variables, these values
ranged from 0.47 to 0.68 SD (0.59 to 0.86 SD after Holm adjustment). These values were calculated post
hoc, from the SEs and not from the observed differences. This also explains why the weight difference
by prior-year emergency department visit (3.8 kg, about 0.6 SD) could appear in a subgroup test, which
examines one variable at a time, while the calibration and rank-weighted tests, which combine all
candidate modifiers, were null. This does not show that the difference is real: the subgroup test was
repeated for seven variables and four outcomes, and after adjustment for all of these comparisons,
P = 0.05.

**Analysis plan:** the plan (outcomes, modifiers, tests, decision rule, multiplicity handling) was
written and fixed before the final analysis (drafted 21 October 2026; approved 28 October 2026), after
earlier exploratory versions of the analysis; it was not registered. Deviation: the plan states that
when the calibration and rank-weighted tests are null, subgroup results are reported as consistent
with a uniform effect; the weight x prior-year emergency department difference is nevertheless
reported, as a nominal, post hoc finding.

## 3. Replication of the published trial results (check of data handling)

Because the analysis code was written with a large language model (section 4), I first checked that
the public data were read and used correctly by reproducing the published trial results. The
investigators' Stata code was translated to Python. The primary results were reproduced exactly:
HbA1c −0.10 percentage points (95% CI, −0.46 to 0.25) and any emergency department visit or
hospitalization −7.13 percentage points (95% CI, −15.30 to 1.04). The investigators' own code was also
run in Stata (versions 19.5 and 19.0, after a one-line change needed on current Stata versions); the
Python translation matched it to the precision of Stata's exported results for every estimate not
involving imputation (absolute differences ≤ 0.00005; relative ≤ 0.002% for health care use and costs).
The published adjusted weight estimate depends on the imputation of five missing baseline weights;
this analysis uses observed baseline weights (complete case), so its average weight effect differs
from the published one. Key estimates are shown in Table S2.

## 4. Use of artificial intelligence in the analysis

**Tool:** Claude Code (Anthropic), a command-line agent, with the Claude Opus 5.5 model, used from June
2026 to November 2026.

**Use:** under my direction, the model wrote and ran the analysis code, reviewed earlier versions of
the analysis (an in-sample bias in a first version was identified and corrected by moving to
out-of-fold estimates), searched for related literature and data, proposed candidate references, and
drafted the text. I set the questions, the analysis plan, and every analytic decision.

**Checks:** (1) replication of the published results in Python and Stata (section 3); (2) an
independent reanalysis of the main findings by a separate model instance with no access to the
analysis code, including a cross-check of the causal forest in R (grf); (3) reproducibility: all
analyses were rerun from scratch on a second computer and operating system, with identical sample
sizes, average effects, and test decisions.

**Sequence of prompts (summary):** (1) June 2026: assess whether the public trial data allow an
analysis of heterogeneous treatment effects, and write a first analysis (superseded); (2) 21 October
2026: review the first analysis (an in-sample bias was found), write the revised analysis plan and
code, add the secondary outcomes, use complete cases, align the subgroup analysis with the trial's
prespecified subgroups, verify the main findings with an independent model instance, and reproduce the
published results in Python and Stata; (3) 22-27 October 2026: rerun all analyses and the Stata check
on a second computer; (4) 28-29 October 2026: search for related
work and other open food-as-medicine data, propose candidate references, draft the manuscript, and
revise it in response to my comments; (5) 29 October to 2 November 2026: independent verification by
separate model instances (a rerun of all analyses, a blind rebuild of the analysis from the plan, a
comparison of the two, an audit of the manuscript against the outputs, and a critical review), and
revision of the manuscript; (6) 9 November 2026: revision in response to my comments, rerun of the
interaction test with all modifiers as specified in the analysis plan, and preparation of the public
repository; (7) 10 November 2026: revision in response to my comments, a calculation of the smallest
difference detectable by the subgroup test, and update of the public repository.

**Prompts:** the prompts from 22 October 2026 onward directly related to the analysis and the
derivation of the manuscript, in sequence, are in [`prompts/prompt_log.md`](prompts/prompt_log.md); the
prompts of the earlier sessions (June 2026 and 21 October 2026) were not retained.

## Table S1. Other randomized food-as-medicine trials with publicly deposited data (search of 28 October 2026)

| Repository record | Clinical outcomes as end points? |
|---|---|
| Hollis-Hansen, Kelseanna, 2025, "Medically Tailored Groceries and Food Resource Coaching", https://doi.org/10.18738/T8/NQIAPM, Texas Data Repository, V4 | No |
| Ranjit, Nalini. "Enhancing the Effectiveness of Home-Delivery Based Produce Prescription Program Implementation Strategies." Ann Arbor, MI: Inter-university Consortium for Political and Social Research (ICPSR) [distributor], August 3, 2026. https://doi.org/10.3886/E238486V2 | No |
| VandenBerg, Alexandra E. "Enhancing Food as Medicine Interventions for Food Insecure Postpartum Women in Central Texas." Ann Arbor, MI: Inter-university Consortium for Political and Social Research (ICPSR) [distributor], August 3, 2026. https://doi.org/10.3886/E237851V1 | No |
| Runge, Carly, 2026, "Improving Diet and Quality of Life through Medically Tailored Meals and Dietitian Counseling in Heart Failure", https://doi.org/10.7910/DVN/Y5AH2R, Harvard Dataverse, V1 | No |

Searched: Harvard Dataverse (which indexes openICPSR, Texas Data Repository, and other Dataverse
installations), ICPSR, Vivli, BioLINCC, Zenodo, Dryad, OSF, and PubMed (Zenodo through web search
only). The Dallas dataset could not be inspected (access page blocked); its content is taken from the
trial's published protocol (PubMed ID 39753253) and registration (NCT06242808).

## Table S2. Reproduction of the published trial estimates (check of data handling)

| Outcome at six months (intervention minus wait list) | Published | Python translation | Investigators' code run in Stata |
|---|---|---|---|
| HbA1c, percentage points (primary outcome) | −0.10 (−0.46 to 0.25) | −0.10 (−0.46 to 0.25) | −0.10 (−0.46 to 0.25) |
| Emergency department visit or hospitalization, percentage points | −7.13 (−15.30 to 1.04) | −7.13 (−15.30 to 1.04) | −7.13 (−15.30 to 1.04) |
| Weight, kg (five missing baseline weights multiply imputed)<sup>a</sup> | 1.95 (0.07 to 3.83) | 1.63 to 1.94<sup>b</sup> | 1.27 to 1.84<sup>c</sup> |

Values are adjusted differences (95% CI) from the trial's models (its Tables 3 and 4). Python and Stata
agreed to within 0.00005 (absolute) for every estimate not involving imputation in the trial's Tables
2 to 4. Of the 118 published cells of the trial's Tables 1 to 4, 93 were reproduced exactly at the
published rounding and 10 multiply imputed cells within the variation expected from imputation; the
other 15 differed slightly (13 lipid or diastolic blood pressure cells, mostly sample sizes; 1 count;
1 rounding).

<sup>a</sup> Multiply imputed estimates vary with the imputation draws. With observed baseline weights
only (395 patients), the investigators' model gave 0.84 kg (95% CI, −0.32 to 2.00).
<sup>b</sup> Range over 10 imputation seeds.
<sup>c</sup> Range over 40 runs that differed only in the sort order of the data (the investigators'
imputation seed was kept).
