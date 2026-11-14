"""Agent C: number-by-number comparison of the original pipeline (outputs/v3) and B's blind
rebuild (verification/multiagent_2026-10-29/B_rebuild). Reads files only. -> comparison.csv

Verdict rule per row
  estimates : AGREE if |diff| <= 0.25 * SE (SE from the original) and the CI-excludes-0 status is
              the same; DIFFERENT NUMBERS SAME CONCLUSION if only the status is the same;
              CONCLUSION DIFFERS otherwise.
  p-values  : AGREE if same side of 0.05 and (|diff| <= 0.05 or both < 0.01);
              DIFFERENT NUMBERS SAME CONCLUSION if same side of 0.05; CONCLUSION DIFFERS otherwise.
  booleans  : AGREE if equal, else CONCLUSION DIFFERS.
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
B = HERE.parent / "B_rebuild"
OUTS = ["hba1c", "weight", "sbp", "ed_ip"]
rows = []


def add(o, item, orig, reb, kind, se=None, lo=None, note=""):
    if kind == "bool":
        v = "AGREE" if bool(orig) == bool(reb) else "CONCLUSION DIFFERS"
    elif kind == "p":
        same = (orig < 0.05) == (reb < 0.05)
        close = abs(orig - reb) <= 0.05 or (orig < 0.01 and reb < 0.01)
        v = ("AGREE" if close else "DIFFERENT NUMBERS SAME CONCLUSION") if same else "CONCLUSION DIFFERS"
    elif kind == "ci":       # orig/reb = (est, lo, hi)
        so = orig[1] > 0 or orig[2] < 0
        sb = reb[1] > 0 or reb[2] < 0
        se_o = (orig[2] - orig[1]) / 3.92
        close = abs(orig[0] - reb[0]) <= 0.25 * se_o
        v = ("AGREE" if close else "DIFFERENT NUMBERS SAME CONCLUSION") if so == sb else "CONCLUSION DIFFERS"
    elif kind == "est":
        close = se is not None and abs(orig - reb) <= 0.25 * se
        v = "AGREE" if close else "DIFFERENT NUMBERS SAME CONCLUSION"
    else:
        v = "NOT RUN BY B" if kind == "na" else ("AGREE" if orig == reb else "DIFFERS")
    fmt = (lambda t: f"{t[0]:+.3f} ({t[1]:+.3f} to {t[2]:+.3f})") if kind == "ci" else (
        lambda x: f"{x:.4g}" if isinstance(x, float) else str(x))
    rows.append(dict(outcome=o, item=item, original=fmt(orig), rebuild_B=fmt(reb), verdict=v, note=note))


avgB = pd.read_csv(B / "average_effects_and_tests.csv").set_index("outcome")
subB = pd.read_csv(B / "subgroups.csv")
pcB = pd.read_csv(B / "positive_control.csv")
wB = pd.read_csv(B / "weight_x_prior_ed.csv").iloc[0]
VAR = {"Baseline HbA1c": "hba1c_cat_control", "Site": "rural_site", "Sex": "male",
       "Prior-year inpatient admission": "ip_i_0mo", "Prior-year ED visit": "ed_i_0mo",
       "Insulin prescribed": "drug_insulin_0mo", "Metformin prescribed": "drug_metformin_0mo"}

for o in OUTS:
    r = json.loads((ROOT / "outputs" / "v3" / o / "results.json").read_text())
    b = avgB.loc[o]
    add(o, "N", r["n"], int(b.n), "eq")
    add(o, "n treated", r["n_treat"], int(b.n_treat), "eq")
    a = r["ate_aipw"]
    add(o, "AIPW average effect (95% CI)", (a["est"], a["lo"], a["hi"]), (b.ate, b.ate_lo, b.ate_hi), "ci",
        note="B primary = mean of repeat-averaged scores; original = median of repeat means")
    add(o, "BLP beta2 (median)", r["blp"]["beta2"], b.blp_beta2, "est", se=r["blp"]["beta2_se"])
    add(o, "BLP 1-sided p", r["blp"]["p_one_sided"], b.blp_p_one_sided, "p")
    add(o, "AUTOC (median)", r["rate_autoc"]["est"], b.autoc, "est", se=r["rate_autoc"]["se"],
        note="sign conventions differ (orig ranks by -tau, benefit -Gamma; B ranks by tau, benefit Gamma); same H1")
    add(o, "AUTOC 1-sided p", r["rate_autoc"]["p_one_sided"], b.rate_p_one_sided, "p")
    add(o, "Heterogeneity detected", r["heterogeneity_detected"], bool(b.heterogeneity_detected), "bool")
    for k, col in zip((1, 2, 3), ("gates_low", "gates_mid", "gates_high")):
        g = r["gates"][k - 1]
        add(o, f"GATES G{k} (G1 = lowest predicted CATE)", g["est"], b[col], "est", se=g["se"])
    add(o, "Interaction test asymptotic p", r["ancova_interactions_joint_p"], b.interaction_p_asymptotic, "p",
        note="orig: 19 terms, F; B: all 21/22 terms, chi2 (see c_tests.csv)")
    add(o, "Interaction test randomization p", r["ancova_interactions_randomization_p"],
        b.interaction_p_randomization, "p", note="covariate set differs (see c_tests.csv)")
    add(o, "Joint BLP of AIPW scores on covariates p", r["blp_covariates_joint_p"], b.blp_covariates_joint_p, "p",
        note="covariate set + F/chi2 + different AIPW scores")
    so = pd.read_csv(ROOT / "outputs" / "v3" / o / "table_subgroups.csv")
    so = so[so.family == "Doyle pre-specified"]
    for name, var in VAR.items():
        sub = so[so.subgroup == name].reset_index(drop=True)
        sb = subB[(subB.outcome == o) & (subB.variable == var)].sort_values("level").reset_index(drop=True)
        for lev in (0, 1):
            add(o, f"{name} [{sub.level[lev]}] n", int(sub.n[lev]), int(sb.n[lev]), "eq")
            add(o, f"{name} [{sub.level[lev]}] effect (95% CI)", (sub.est[lev], sub.lo[lev], sub.hi[lev]),
                (sb.est[lev], sb.lo[lev], sb.hi[lev]), "ci")
        add(o, f"{name} interaction p", sub.p_interaction[0], sb.p_interaction[0], "p")
        add(o, f"{name} Holm p", sub.p_holm[0], sb.p_holm[0], "p")
    sens = pd.read_csv(ROOT / "outputs" / "v3" / o / "table_sensitivity.csv")
    pcs = sens[sens.analysis.str.startswith("PC")].reset_index(drop=True)
    for i, k in enumerate((0.25, 0.5, 1.0)):
        po = pcs.iloc[i]
        pb = pcB[(pcB.outcome == o) & (pcB.delta_sd == k)].iloc[0]
        add(o, f"PC {k} SD BLP p", po.blp_p_one_sided, pb.blp_p_one_sided, "p",
            note="orig 3 repeats, +-delta/2; B 5 repeats, +delta to older treated")
        add(o, f"PC {k} SD AUTOC p", po.autoc_p_one_sided, pb.rate_p_one_sided, "p")
        add(o, f"PC {k} SD detected", bool(po.blp_p_one_sided < .05 and po.autoc_p_one_sided < .05),
            bool(pb.detected), "bool")

w = json.loads((ROOT / "outputs" / "v3" / "weight" / "weight_ed_checks.json").read_text())
add("weight", "Weight x prior ED: OLS interaction (95% CI)",
    (w["ols_interaction"], w["ols_interaction"] - 1.96 * w["ols_se"], w["ols_interaction"] + 1.96 * w["ols_se"]),
    (wB.coef_interaction, wB.lo, wB.hi), "ci")
add("weight", "Weight x prior ED: OLS p", w["ols_p_hc1"], wB.p_hc1, "p", note="orig normal z; B t(384)")
add("weight", "Weight x prior ED: permutation p, sharp null", w["perm_p_sharp_null"], wB.p_perm_sharp_null, "p",
    note="orig 10000 draws; B 2000")
add("weight", "Weight x prior ED: permutation p, constant-effect null", w["perm_p_constant_effect_null"],
    wB.p_perm_constant_effect, "p")
add("weight", "Weight x prior ED: trimmed (|dY|<=10 kg) p", w["trim10_p_hc1"], "not run", "na",
    note="not run by B (post hoc check, not in plan)")
out = pd.DataFrame(rows)
out.to_csv(HERE / "comparison.csv", index=False)
pd.set_option("display.width", 250, "display.max_colwidth", 60, "display.max_rows", 500)
print(out.drop(columns="note").to_string())
print(out.verdict.value_counts())
