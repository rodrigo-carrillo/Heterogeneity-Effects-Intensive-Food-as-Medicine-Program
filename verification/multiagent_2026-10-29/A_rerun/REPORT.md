# Agent A: clean rerun (29 Oct 2026) — INCOMPLETE

Saved by the main session from Agent A's final message (the agent's tool was blocked from writing
this file). Content is the agent's, condensed.

**Verdict:** published-trial replication REPRODUCED (7 files byte-identical, two clean runs). HTE
pipeline NOT VERIFIED: only `hba1c` (main) completed and it matched the committed outputs exactly
(3,207 values, 0 differences > 1e-6; figures byte-identical; console log differs only in line
order). The other 7 runs, the summary, `manuscript/figure1.png` and the determinism check were
not done.

## Setup
- HEAD `2542993`, clean tree. Python 3.13.14, econml 0.17.0, scikit-learn 1.9.1, numpy 2.5.4,
  pandas 3.0.6; Windows 11, 8 cores, 16 GB RAM (shared with agents B and E).
- Clean copy in the session scratchpad (`agentA_clean`); commands `python code/replicate_doyle.py`,
  `python code/run_all.py --jobs 1`.

## What happened
1. Attempt 1: replication 10 s (identical); `hba1c` 19 min 06 s, no worker crashes;
   `hba1c_ldl_cc` failed: `FileNotFoundError` writing
   `outputs\v3\hba1c_ldl_cc\table_variable_importance.csv`, a 260-character path
   (Windows MAX_PATH; `LongPathsEnabled` = 0). Cause: the long scratchpad path, not the analysis.
   **Note:** the pipeline fails on Windows when the project sits deep in a folder tree.
2. Attempt 2 (`subst Q:` short path): replication identical again; the `hba1c` run was stopped by
   Claude Code because the machine was critically low on memory (~2.2 of 15.8 GB free). Not
   restarted (per instructions).

## Numbers confirmed by the clean run
HbA1c: N = 341 (165 treated); AIPW −0.05 (−0.44 to 0.33); BLP p 0.9989, AUTOC p 0.9984;
randomization p 0.802; smallest Holm p 1.0. Every `TABLE_ROWS` cell in `letter_text.py` matches
the committed results (only the HbA1c row confirmed by a fresh rerun).

## Untraceable manuscript numbers
"Linear-model and permutation P = .002 and .003" (weight × prior ED) come from docs/04 §4; no
committed script computes them for N = 392 (`verification/independent_agent/v3_subgroups.py` used
the earlier N = 395). Recomputed with `check_weight_ed_ols_perm.py` on N = 392: OLS interaction
−3.68 (SE 1.22), p = 0.0025; permutation (10 000, within strata) p = 0.0025 vs 0.0028 in docs/04,
within Monte Carlo error (SE ≈ 0.0005).

## To finish
Run from a short path with enough free memory (~19 min per run, ~2.5 h total):
`run_all.py --jobs 1`; `manuscript_figure.py`; a second `hte_food_as_medicine.py --outcome weight`
(determinism); then `compare_outputs.py <project> <clean> diffs.csv`.

Files: `compare_outputs.py` (validated with planted differences of 1e-8, 0.02, 0.1),
`check_weight_ed_ols_perm.py`, `differences_attempt1_hba1c.csv`,
`differences_replication_attempt2.csv`.
