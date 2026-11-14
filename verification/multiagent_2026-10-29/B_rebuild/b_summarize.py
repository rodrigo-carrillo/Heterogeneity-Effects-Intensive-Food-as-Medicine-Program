"""Collect out/*.json into results.json and tidy CSVs."""
import json
from pathlib import Path
import pandas as pd

H = Path(__file__).resolve().parent
OUTCOMES = ["hba1c", "weight", "sbp", "ed_ip"]
res, avg, sub, pc = {}, [], [], []
for o in OUTCOMES:
    f = H / "out" / f"{o}.json"
    if not f.exists():
        continue
    r = json.loads(f.read_text())
    res[o] = r
    om, it = r["omnibus"], r["interaction_test"]
    avg.append(dict(outcome=o, n=r["n"], n_treat=r["n_treat"],
                    ate=r["ate_aipw"]["est"], ate_lo=r["ate_aipw"]["lo"], ate_hi=r["ate_aipw"]["hi"],
                    ate_p=r["ate_aipw"]["p"],
                    ate_median_agg=r["ate_aipw_median_agg"]["est"],
                    ate_median_agg_lo=r["ate_aipw_median_agg"]["lo"], ate_median_agg_hi=r["ate_aipw_median_agg"]["hi"],
                    ate_e_full465=r["ate_aipw_e_full465"]["est"],
                    diff_means=r["ate_unadjusted_diff_means"],
                    blp_beta2=om["blp_beta2_median"], blp_p_one_sided=om["blp_p_median"],
                    autoc=om["autoc_median"], rate_p_one_sided=om["rate_p_median"],
                    heterogeneity_detected=om["decision_detected"],
                    gates_low=om["gates_median"][0], gates_mid=om["gates_median"][1], gates_high=om["gates_median"][2],
                    interaction_wald=it["wald"], interaction_df=it["df"],
                    interaction_p_asymptotic=it["p_asymptotic"], interaction_p_randomization=it["p_randomization"],
                    blp_covariates_joint_p=r["blp_covariates_joint"]["p"]))
    sub += r["subgroups"]
    fp = H / "out" / f"{o}_pc.json"
    if fp.exists():
        rp = json.loads(fp.read_text())
        res[o]["positive_control"] = {k: v for k, v in rp.items() if k in ("pc", "pc_older_def", "pc_older_n", "y_sd")}
        for dlt, om2 in rp["pc"].items():
            pc.append(dict(outcome=o, delta_sd=float(dlt), delta_units=float(dlt) * rp["y_sd"],
                           older_def=rp["pc_older_def"], blp_beta2=om2["blp_beta2_median"],
                           blp_p_one_sided=om2["blp_p_median"], autoc=om2["autoc_median"],
                           rate_p_one_sided=om2["rate_p_median"], detected=om2["decision_detected"]))
(H / "results.json").write_text(json.dumps(res, indent=1))
pd.DataFrame(avg).to_csv(H / "average_effects_and_tests.csv", index=False)
cols = ["outcome", "subgroup", "variable", "level", "n", "n_treat", "est", "se", "lo", "hi",
        "diff_1_minus_0", "p_interaction", "p_holm"]
pd.DataFrame(sub)[cols].to_csv(H / "subgroups.csv", index=False)
if pc:
    pd.DataFrame(pc).to_csv(H / "positive_control.csv", index=False)
if "weight" in res:
    pd.DataFrame([res["weight"]["weight_x_prior_ed"]]).to_csv(H / "weight_x_prior_ed.csv", index=False)
pd.set_option("display.width", 250)
print(pd.DataFrame(avg).round(4).T)
s = pd.DataFrame(sub)
print(s[cols].round(3).to_string())
if pc:
    print(pd.DataFrame(pc).round(4).to_string())
