"""Agent B independent rebuild: data preparation (from docs/02_analysis_plan.md + data dictionary only)."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
CSV = ROOT / "data" / "raw" / "dataverse_files" / "analysis_sample.csv"
LB_PER_KG = 2.20462

CONT = ["age", "hba1c_val_0mo", "days_a1c_to_consent", "bmi_val_0mo", "sys_val_0mo",
        "dias_val_0mo", "op_n_0mo", "ed_or_ip_n_0mo"]
BIN = ["male", "hispanic", "race_nonhisp_black", "rural_site", "hba1c_cat_control",
       "drug_insulin_0mo", "drug_metformin_0mo", "drug_glp1_0mo", "drug_sulf_0mo",
       "drug_other_0mo", "drug_any_db_0mo", "ip_i_0mo", "ed_i_0mo"]
MODIFIERS = CONT + BIN  # 21

SUBGROUPS = {  # Doyle pre-specified, reproducible from public data
    "hba1c_ge_9.5": "hba1c_cat_control",
    "rural_site": "rural_site",
    "male": "male",
    "prior_ip": "ip_i_0mo",
    "prior_ed": "ed_i_0mo",
    "insulin": "drug_insulin_0mo",
    "metformin": "drug_metformin_0mo",
}
OUTCOMES = ["hba1c", "weight", "sbp", "ed_ip"]


def load_raw():
    d = pd.read_csv(CSV).copy()
    quirks = {}
    # placeholder triglyceride (not used in main runs, cleaned per plan)
    quirks["trigly_placeholder_n"] = int((d["trigly_val_0mo"] == 9999999).sum())
    d.loc[d["trigly_val_0mo"] == 9999999, "trigly_val_0mo"] = np.nan
    d = d.copy()
    d["days_a1c_to_consent"] = d["consent_date"] - d["hba1c_date_0mo"]
    d["weight_kg_0mo"] = d["weight_lb_val_0mo"] / LB_PER_KG
    d["weight_kg_6mo"] = d["weight_lb_val_6mo"] / LB_PER_KG
    d["stratum"] = 2 * d["rural_site"] + d["hba1c_cat_control"]
    quirks["kg_column_ratio_lb_over_kg"] = float((d.weight_lb_val_0mo / d.weight_val_0mo).median())
    quirks["hba1c_cat_control_vs_ge9.5_mismatch_n"] = int(
        ((d.hba1c_val_0mo >= 9.5).astype(int) != d.hba1c_cat_control).sum())
    quirks["drug_any_db_vs_summary_mismatch_n"] = int((d.drug_any_db_0mo != d.drug_any_db_summary_0mo).sum())
    # known (design) propensity per stratum from all 465 randomized+analysed patients
    e_tab = d.groupby("stratum")["treat_now"].mean()
    d["e_known"] = d["stratum"].map(e_tab)
    return d, quirks, e_tab.to_dict()


def analytic(d, outcome):
    X_cols = list(MODIFIERS)
    if outcome == "hba1c":
        y = d["hba1c_val_6mo"] - d["hba1c_val_0mo"]
    elif outcome == "weight":
        y = d["weight_kg_6mo"] - d["weight_kg_0mo"]
        X_cols = X_cols + ["weight_kg_0mo"]
    elif outcome == "sbp":
        y = d["sys_val_6mo"] - d["sys_val_0mo"]
    elif outcome == "ed_ip":
        y = 100.0 * d["ed_or_ip_i_6mo"]
    else:
        raise ValueError(outcome)
    df = d.copy()
    df["Y"] = y
    keep = df["Y"].notna() & df[X_cols].notna().all(axis=1)
    df = df.loc[keep].reset_index(drop=True)
    return df, X_cols


if __name__ == "__main__":
    d, q, e = load_raw()
    print(q, e)
    for o in OUTCOMES:
        df, xc = analytic(d, o)
        print(o, len(df), int(df.treat_now.sum()), len(xc), "Y sd", round(df.Y.std(), 3))
