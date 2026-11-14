# Heterogeneity in the effects of an intensive food-as-medicine program

Code, aggregate outputs and checks for an exploratory secondary analysis of the randomized clinical
trial by Doyle et al., *Effect of an Intensive Food-as-Medicine Program on Health and Health Care Use:
A Randomized Clinical Trial*, JAMA Intern Med 2024;184(2):154-163 (ClinicalTrials.gov NCT03718832).

The trial found no average effect of the program on 6-month hemoglobin A1c (HbA1c). This analysis asks
whether baseline characteristics modified the program's effects on HbA1c (primary outcome), weight,
systolic blood pressure, and any emergency department visit or hospitalization at 6 months. It uses
causal forests with out-of-fold (cross-fitted) estimates, a prespecified decision rule (calibration
test and rank-weighted average treatment effect, both 1-sided P < .05), and the subgroup variables the
trial investigators prespecified. All analyses are complete case.

This repository accompanies a manuscript on this analysis. The detailed methods and supporting tables
that the article format leaves out are in [`METHODS.md`](METHODS.md).

## Data (not included)
The trial data are publicly posted by the trial investigators on Harvard Dataverse and are **not**
redistributed here:

> Doyle J, Alsan M, Skelley N, Lu Y, Cawley J. Replication Data for: Effect of an Intensive
> Food-as-Medicine Program on Health and Health Care Use: A Randomized Clinical Trial. Harvard
> Dataverse, version 2. doi:[10.7910/DVN/LZUTU0](https://doi.org/10.7910/DVN/LZUTU0)

Download the dataset files and place them in `data/raw/dataverse_files/` (see
[`data/README.md`](data/README.md)). The code expects `data/raw/dataverse_files/analysis_sample.csv`.

## Environment
Python 3.13 with the pinned packages in `requirements.txt` (EconML 0.17, scikit-learn 1.9):
```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```
The Stata check (`verification/stata/`) needs Stata 16.1 or later and the investigators' own
`01_published_tables.do` from the Dataverse record.

## How to run
```bash
python code/run_all.py --fast --jobs 4   # smoke test of all 8 runs -> outputs/v3_fast/ (minutes)
python code/run_all.py --jobs 2          # full analysis: 4 outcomes x 2 covariate sets -> outputs/v3/ (about 1-2 h)
python code/hte_food_as_medicine.py --outcome weight   # one run; flags: --outcome {hba1c,weight,sbp,ed_ip} --ldl --fast
python code/replicate_doyle.py           # reproduction of the published trial results -> outputs/replication/
python code/weight_ed_checks.py          # post hoc checks of weight x prior-year ED visit (after the weight run)
python code/subgroup_detectable_difference.py   # smallest difference detectable by the subgroup tests (post hoc)
python code/manuscript_figure.py --leading-zero  # the manuscript Figure (from outputs/v3/)
python code/checks/v1_leak_check.py      # why an earlier in-sample version overstated heterogeneity
```
Seeds are fixed. Forest-based numbers can differ slightly across platforms (third decimal); the
decisions do not.

## Repository layout
```
code/
  hte_food_as_medicine.py   main analysis, one outcome x covariate set per run
  hte_lib.py                model fitting; cross-fitting repeats run in separate worker processes
  run_all.py                all 8 runs + cross-outcome summary
  replicate_doyle.py        Python translation of the investigators' Stata code (published Tables 1-4)
  weight_ed_checks.py       post hoc checks of the weight x prior-year ED difference
  subgroup_detectable_difference.py  minimum detectable difference of the subgroup tests (post hoc)
  manuscript_figure.py      Figure of the manuscript
  build_notebook.py         regenerates notebooks/hte_food_as_medicine.ipynb from the script
  checks/v1_leak_check.py   in-sample vs out-of-fold calibration slope (first version of the analysis)
notebooks/                  notebook version of the main script (generated; no outputs stored)
outputs/
  v3/<run>/                 results of each run (see below); <run> = hba1c, weight, sbp, ed_ip, and *_ldl_cc
  v3/summary_across_outcomes.{md,csv}
  v3/subgroup_detectable_difference.{md,csv}
  replication/              published vs reproduced trial estimates
verification/               independent checks (see below)
prompts/                    the author's prompts to the language model (analysis, data, manuscript), in sequence
```

Each `outputs/v3/<run>/` folder holds `results_summary.md` (all key numbers), `results.json`,
`console_log.txt`, tables (`table_balance.csv`, `table_blp_covariates.csv`, `table_subgroups.csv`,
`table_sensitivity.csv`, `table_variable_importance.csv`, `table_split_stability.csv`) and figures
`fig1`-`fig6`. The pipeline also writes `patient_level_cate.csv` (one row per patient); those files
contain individual-level trial data and are therefore not shared here; the code regenerates them.
`*_ldl_cc` runs add baseline LDL cholesterol (patients with LDL measured).

## Verification
The analysis code was written with a large language model (Claude Code, Anthropic; Claude Opus 5.5)
under the author's direction. It was checked by:
- `outputs/replication/` and `verification/replication_checks/`: reproduction of the published trial
  results from the public data (Python translation of the investigators' Stata code);
- `verification/stata/`: the investigators' own code run in Stata and compared with the Python
  translation;
- `verification/independent_agent/`: an independent reanalysis by a separate model instance without
  access to the analysis code, including a cross-check of the causal forest in R (`grf`);
- `verification/multiagent_2026-10-29/`: a rerun of all analyses (A), a blind rebuild of the analysis
  from the analysis plan (B), a comparison of the two (C), and a critical review (E).

Some verification reports refer to internal documents (analysis plan, working notes) that are not
part of this repository. Files that list individual patients were left out.

## Use of artificial intelligence
- **Tool:** Claude Code (Anthropic), a command-line agent, with the Claude Opus 5.5 model.
- **When:** June to November 2026.
- **Where and how:** under the author's direction, the model wrote and ran the analysis code; checked
  it (reproduction of the published trial results, an independent reanalysis by a separate model
  instance, and reruns on a second computer; see Verification); searched for related literature and
  data; proposed candidate references, which the author verified in PubMed; and drafted and edited
  the manuscript text. The author set the questions, the analysis plan, and every analytic decision,
  and reviewed all content.

## License
Code and outputs: MIT License (`LICENSE`). The trial data are not included; they are available from
the Harvard Dataverse record above, under its terms.
