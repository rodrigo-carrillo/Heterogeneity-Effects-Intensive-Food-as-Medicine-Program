"""Shared data prep replicating code/hte_food_as_medicine.py §1 (no project code is modified)."""
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parents[1]
CONT = ["age", "hba1c_val_0mo", "days_a1c_to_consent", "bmi_val_0mo",
        "sys_val_0mo", "dias_val_0mo", "op_n_0mo", "ed_or_ip_n_0mo"]
BIN = ["male", "hispanic", "race_nonhisp_black", "rural_site", "hba1c_cat_control",
       "drug_insulin_0mo", "drug_metformin_0mo", "drug_glp1_0mo", "drug_sulf_0mo",
       "drug_other_0mo", "drug_any_db_0mo", "ip_i_0mo", "ed_i_0mo"]
VAL = {"hba1c": "hba1c_val", "weight": "weight_val", "sbp": "sys_val", "ed_ip": "ed_or_ip_i"}
PRESPEC = ["hba1c_cat_control", "rural_site", "male", "ip_i_0mo", "ed_i_0mo",
           "drug_insulin_0mo", "drug_metformin_0mo"]


def load_raw():
    raw = pd.read_csv(ROOT / "data/raw/dataverse_files/analysis_sample.csv")
    raw.loc[raw["trigly_val_0mo"] >= 9_999_999, "trigly_val_0mo"] = np.nan
    for m in (0, 6, 12):
        raw[f"weight_val_{m}mo"] = raw[f"weight_lb_val_{m}mo"] / 2.20462
    raw["days_a1c_to_consent"] = raw["consent_date"] - raw["hba1c_date_0mo"]
    raw["stratum"] = 2 * raw["rural_site"] + raw["hba1c_cat_control"]
    return raw.copy()


def features(outcome, ldl=False):
    return (CONT + (["weight_val_0mo"] if outcome == "weight" else []) + BIN
            + (["ldl_val_0mo"] if ldl else []))


def analytic(raw, outcome, ldl=False):
    v = VAL[outcome]
    r = raw.copy()
    if outcome == "ed_ip":
        r["Y"] = 100.0 * r[f"{v}_6mo"]
    else:
        r["Y"] = r[f"{v}_6mo"] - r[f"{v}_0mo"]
    F = features(outcome, ldl)
    keep = r["Y"].notna() & r[F].notna().all(axis=1)
    df = r[keep].copy().reset_index(drop=True)
    df["E_HAT"] = df.groupby("stratum")["treat_now"].transform("mean")
    for k in (1, 2, 3):
        df[f"st{k}"] = (df.stratum == k).astype(int)
    return df, F
