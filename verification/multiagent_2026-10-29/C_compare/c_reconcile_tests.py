"""Agent C: reconcile the parametric tests between the original pipeline and B's rebuild.

Light computation only (OLS on <=456 rows; no forests). For each main outcome:
  1. Rebuild the analytic sample the original way and check it equals B's (subject ids, Y).
  2. Interaction test (Y ~ W + X + W*X, HC1 joint Wald on the interactions) under
     covariate set O (original: 19 terms, drug_any_db_0mo, hba1c_cat_control and baseline
     weight dropped) and B (all 21/22 modifiers), asymptotic p by F (original) and chi2 (B),
     and randomization p (2000 within-stratum re-randomizations, constant-effect null), with
     the original's seed/algorithm for set O (should reproduce outputs/v3) and B's for set B
     (should reproduce B). Also O + each of the two dropped binaries separately.
  3. Joint BLP of the ORIGINAL repeat-averaged AIPW scores on covariates, set O vs set B,
     F vs chi2 (isolates the spec effect from the AIPW-score effect).
  4. Collinearity diagnostics for the two dropped binaries.
Writes c_tests.csv and c_collinearity.csv next to this file.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / "B_rebuild"))
from b_data import load_raw, analytic  # noqa: E402  (B's data prep, for the equality check)

CONT = ["age", "hba1c_val_0mo", "days_a1c_to_consent", "bmi_val_0mo", "sys_val_0mo",
        "dias_val_0mo", "op_n_0mo", "ed_or_ip_n_0mo"]
BIN = ["male", "hispanic", "race_nonhisp_black", "rural_site", "hba1c_cat_control",
       "drug_insulin_0mo", "drug_metformin_0mo", "drug_glp1_0mo", "drug_sulf_0mo",
       "drug_other_0mo", "drug_any_db_0mo", "ip_i_0mo", "ed_i_0mo"]
DROP_O = ["drug_any_db_0mo", "hba1c_cat_control", "weight_val_0mo"]
SEED_O, SEED_B = 20260621, 20261029
N_PERM = 2000
VAL = {"hba1c": "hba1c_val", "weight": "weight_val", "sbp": "sys_val", "ed_ip": "ed_or_ip_i"}


def orig_sample(outcome):
    raw = pd.read_csv(ROOT / "data" / "raw" / "dataverse_files" / "analysis_sample.csv")
    for m in (0, 6, 12):
        raw[f"weight_val_{m}mo"] = raw[f"weight_lb_val_{m}mo"] / 2.20462
    raw["days_a1c_to_consent"] = raw["consent_date"] - raw["hba1c_date_0mo"]
    raw["stratum"] = 2 * raw["rural_site"] + raw["hba1c_cat_control"]
    v = VAL[outcome]
    raw["Y"] = (100.0 * raw[f"{v}_6mo"]) if outcome == "ed_ip" else raw[f"{v}_6mo"] - raw[f"{v}_0mo"]
    cols = CONT + (["weight_val_0mo"] if outcome == "weight" else []) + BIN
    keep = raw["Y"].notna() & raw[cols].notna().all(axis=1)
    return raw[keep].reset_index(drop=True), cols


def hc1(y, Z):
    n, k = Z.shape
    ZtZi = np.linalg.pinv(Z.T @ Z)
    b = ZtZi @ (Z.T @ y)
    r = y - Z @ b
    cov = ZtZi @ ((Z * (r ** 2)[:, None]).T @ Z) @ ZtZi * n / (n - k)
    return b, cov


def inter_stat(y, X, w):
    D = np.column_stack([np.ones(len(y)), X, w, X * w[:, None]])
    b, cov = hc1(y, D)
    k = X.shape[1] + 2
    bb, V = b[k:], cov[k:, k:]
    return float(bb @ np.linalg.pinv(V) @ bb), D.shape[1]


def interaction_test(y, X, W, strata, seed, algo):
    stat, ncol = inter_stat(y, X, W.astype(float))
    q, n = X.shape[1], len(y)
    p_chi2 = stats.chi2.sf(stat, q)
    p_F = stats.f.sf(stat / q, q, n - ncol)
    Z0 = np.column_stack([np.ones(n), W, X])
    tau0 = hc1(y, Z0)[0][1]
    y0 = y - tau0 * W
    rng = np.random.default_rng(seed)
    cnt = 0
    for _ in range(N_PERM):
        wp = W.copy()
        for s in np.unique(strata):
            idx = np.where(strata == s)[0]
            if algo == "O":
                wp[idx] = rng.permutation(W[idx])
            else:
                wp[idx] = W[rng.permutation(idx)]
        cnt += inter_stat(y0 + tau0 * wp, X, wp.astype(float))[0] >= stat
    return dict(stat=stat, q=q, p_chi2=p_chi2, p_F=p_F, tau0=tau0, p_rand=(1 + cnt) / (1 + N_PERM))


def blp_cov(g, X):
    Xs = (X - X.mean(0)) / X.std(0)
    Z = np.column_stack([np.ones(len(g)), Xs])
    b, cov = hc1(g, Z)
    stat = float(b[1:] @ np.linalg.pinv(cov[1:, 1:]) @ b[1:])
    q = X.shape[1]
    return stats.chi2.sf(stat, q), stats.f.sf(stat / q, q, len(g) - q - 1)


rows, coll = [], []
dB, _, _ = load_raw()
for o in ["hba1c", "weight", "sbp", "ed_ip"]:
    df, cols = orig_sample(o)
    bdf, bcols = analytic(dB, o)
    same_ids = set(df.subject_id) == set(bdf.subject_id)
    m = df[["subject_id", "Y"]].merge(bdf[["subject_id", "Y"]], on="subject_id")
    same_y = bool(np.allclose(m.Y_x, m.Y_y))
    pl = pd.read_csv(ROOT / "outputs" / "v3" / o / "patient_level_cate.csv")
    assert (pl.subject_id.to_numpy() == df.subject_id.to_numpy()).all()
    assert np.allclose(pl.Y.to_numpy(), df.Y.to_numpy())
    Y, W, S = df.Y.to_numpy(float), df.treat_now.to_numpy(int), df.stratum.to_numpy(int)
    colsO = [c for c in cols if c not in DROP_O]
    sets = {"O (19, original)": (colsO, SEED_O + 99, "O"),
            "B (all, rebuild)": (cols, SEED_B + 77, "B"),
            "O + hba1c_cat_control": (colsO + ["hba1c_cat_control"], SEED_O + 99, "O"),
            "O + drug_any_db_0mo": (colsO + ["drug_any_db_0mo"], SEED_O + 99, "O")}
    if o == "weight":
        sets["O + weight_val_0mo"] = (colsO + ["weight_val_0mo"], SEED_O + 99, "O")
    g = pl.gamma_aipw.to_numpy()
    for name, (cs, seed, algo) in sets.items():
        X = df[cs].to_numpy(float)
        it = interaction_test(Y, X, W, S, seed, algo)
        pc2, pF = blp_cov(g, X)
        rows.append(dict(outcome=o, covariate_set=name, n=len(Y), n_terms=it["q"], wald=it["stat"],
                         inter_p_F=it["p_F"], inter_p_chi2=it["p_chi2"], inter_p_rand=it["p_rand"],
                         tau0=it["tau0"], blpcov_origscores_p_F=pF, blpcov_origscores_p_chi2=pc2,
                         same_sample_as_B=same_ids, same_Y_as_B=same_y))
        print(rows[-1], flush=True)
    # collinearity of the two dropped binaries
    Xall = df[cols].to_numpy(float)
    for c in ["drug_any_db_0mo", "hba1c_cat_control"]:
        j = cols.index(c)
        others = np.column_stack([np.ones(len(df)), np.delete(Xall, j, axis=1)])
        bb = np.linalg.lstsq(others, Xall[:, j], rcond=None)[0]
        r2 = 1 - np.var(Xall[:, j] - others @ bb) / np.var(Xall[:, j])
        coll.append(dict(outcome=o, variable=c, r2_on_other_modifiers=r2, vif=1 / (1 - r2)))
    anydrug = (df[["drug_insulin_0mo", "drug_metformin_0mo", "drug_glp1_0mo", "drug_sulf_0mo",
                   "drug_other_0mo"]].max(axis=1) != df.drug_any_db_0mo).sum()
    coll.append(dict(outcome=o, variable="drug_any_db != OR(5 drug flags), n", r2_on_other_modifiers=np.nan,
                     vif=float(anydrug)))
    coll.append(dict(outcome=o, variable="rank(full interaction design) / ncol",
                     r2_on_other_modifiers=np.nan,
                     vif=float(np.linalg.matrix_rank(np.column_stack([np.ones(len(Y)), Xall, W, Xall * W[:, None]])))
                     / (2 * Xall.shape[1] + 2)))
pd.DataFrame(rows).to_csv(HERE / "c_tests.csv", index=False)
pd.DataFrame(coll).to_csv(HERE / "c_collinearity.csv", index=False)
print(pd.DataFrame(coll).to_string())
