"""Agent C: why does the AIPW average effect differ (HbA1c -0.053 original vs -0.085 B)?

Reproduces both implementations' cross-fitted AIPW scores exactly (random-forest outcome
models only, no causal forests, n_jobs=1), then swaps one ingredient at a time:
  - fold partition + RF seeds: original scheme vs B scheme
  - RF max_features: 1.0 (original, sklearn default) vs 1/3 (B)
  - aggregation: median of per-repeat means (original) vs mean of repeat-averaged scores (B)
Usage: python c_reconcile_aipw.py hba1c [sbp ...]   -> c_aipw_<outcome>.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import StratifiedKFold

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SEED_O, SEED_B = 20260621, 20261029
CONT = ["age", "hba1c_val_0mo", "days_a1c_to_consent", "bmi_val_0mo", "sys_val_0mo",
        "dias_val_0mo", "op_n_0mo", "ed_or_ip_n_0mo"]
BIN = ["male", "hispanic", "race_nonhisp_black", "rural_site", "hba1c_cat_control",
       "drug_insulin_0mo", "drug_metformin_0mo", "drug_glp1_0mo", "drug_sulf_0mo",
       "drug_other_0mo", "drug_any_db_0mo", "ip_i_0mo", "ed_i_0mo"]
VAL = {"hba1c": "hba1c_val", "weight": "weight_val", "sbp": "sys_val", "ed_ip": "ed_or_ip_i"}


def orig_sample(outcome):
    """Analytic sample and X columns exactly as code/hte_food_as_medicine.py builds them."""
    raw = pd.read_csv(ROOT / "data" / "raw" / "dataverse_files" / "analysis_sample.csv")
    for m in (0, 6, 12):
        raw[f"weight_val_{m}mo"] = raw[f"weight_lb_val_{m}mo"] / 2.20462
    raw["days_a1c_to_consent"] = raw["consent_date"] - raw["hba1c_date_0mo"]
    raw["stratum"] = 2 * raw["rural_site"] + raw["hba1c_cat_control"]
    v = VAL[outcome]
    raw["Y"] = (100.0 * raw[f"{v}_6mo"]) if outcome == "ed_ip" else raw[f"{v}_6mo"] - raw[f"{v}_0mo"]
    cols = CONT + (["weight_val_0mo"] if outcome == "weight" else []) + BIN   # original column order
    keep = raw["Y"].notna() & raw[cols].notna().all(axis=1)
    return raw[keep].reset_index(drop=True), cols


def scores(Y, W, X, e, scheme, mf):
    reps = []
    strat = df.stratum.to_numpy() * 2 + W
    for r in range(5):
        if scheme == "O":
            skf = StratifiedKFold(5, shuffle=True, random_state=SEED_O + r)
            splits = list(skf.split(X, strat))
        else:
            skf = StratifiedKFold(5, shuffle=True, random_state=SEED_B + 1000 * r)
            splits = list(skf.split(np.zeros(len(Y)), strat))
        mu0, mu1 = np.empty(len(Y)), np.empty(len(Y))
        for k, (tr, te) in enumerate(splits):
            for arm, mu in ((1, mu1), (0, mu0)):
                if scheme == "O":
                    seed = SEED_O + 100 * r + k
                else:
                    seed = SEED_B + 10 * r + k + (0 if arm == 1 else 50)
                idx = tr[W[tr] == arm]
                m = RandomForestRegressor(n_estimators=500, min_samples_leaf=5, max_features=mf,
                                          random_state=seed, n_jobs=1).fit(X[idx], Y[idx])
                mu[te] = m.predict(X[te])
        reps.append(mu1 - mu0 + W * (Y - mu1) / e - (1 - W) * (Y - mu0) / (1 - e))
    return np.array(reps)


rows = []
for o in sys.argv[1:]:
    df, cols = orig_sample(o)
    Y, W = df.Y.to_numpy(float), df.treat_now.to_numpy(int)
    X = df[cols].to_numpy(float)
    e = df.groupby("stratum")["treat_now"].transform("mean").to_numpy()
    pl = pd.read_csv(ROOT / "outputs" / "v3" / o / "patient_level_cate.csv")
    for scheme in ("O", "B"):
        for mf in (1.0, 1 / 3):
            G = scores(Y, W, X, e, scheme, mf)
            means = G.mean(axis=1)
            rows.append(dict(outcome=o, folds_and_seeds=scheme, max_features=round(mf, 3),
                             ate_median_of_repeat_means=float(np.median(means)),
                             ate_mean_of_avg_scores=float(G.mean()),
                             se_sd_over_sqrt_n=float(G.mean(0).std(ddof=1) / np.sqrt(len(Y))),
                             repeat_means=";".join(f"{m:.3f}" for m in means),
                             matches_original_gamma=bool(np.allclose(G.mean(0), pl.gamma_aipw))
                             if (scheme, mf) == ("O", 1.0) else None))
            print(rows[-1], flush=True)
    pd.DataFrame([r for r in rows if r["outcome"] == o]).to_csv(HERE / f"c_aipw_{o}.csv", index=False)
