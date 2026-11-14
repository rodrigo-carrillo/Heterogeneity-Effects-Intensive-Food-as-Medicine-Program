"""Agent C: interaction-test randomization P in the *_ldl_cc runs (quoted in Supplement 1),
original covariate set (drops drug_any_db_0mo, hba1c_cat_control, weight) vs all modifiers.
OLS only. -> c_ldl_interaction.csv"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
from c_reconcile_aipw import orig_sample  # noqa: E402  (pure function; no forests at import)

SEED_O, N_PERM = 20260621, 2000
DROP_O = ["drug_any_db_0mo", "hba1c_cat_control", "weight_val_0mo"]


def hc1(y, Z):
    n, k = Z.shape
    A = np.linalg.pinv(Z.T @ Z)
    b = A @ (Z.T @ y)
    r = y - Z @ b
    return b, A @ ((Z * (r ** 2)[:, None]).T @ Z) @ A * n / (n - k)


def stat(y, X, w):
    D = np.column_stack([np.ones(len(y)), X, w, X * w[:, None]])
    b, V = hc1(y, D)
    k = X.shape[1] + 2
    return float(b[k:] @ np.linalg.pinv(V[k:, k:]) @ b[k:])


rows = []
for o in ["hba1c", "weight", "sbp", "ed_ip"]:
    df, cols = orig_sample(o)
    df = df[df.ldl_val_0mo.notna()].reset_index(drop=True)
    cols = cols + ["ldl_val_0mo"]
    Y, W, S = df.Y.to_numpy(float), df.treat_now.to_numpy(int), df.stratum.to_numpy(int)
    for name, cs in [("O (original)", [c for c in cols if c not in DROP_O]), ("all modifiers", cols)]:
        X = df[cs].to_numpy(float)
        s0 = stat(Y, X, W.astype(float))
        tau0 = hc1(Y, np.column_stack([np.ones(len(Y)), W, X]))[0][1]
        rng = np.random.default_rng(SEED_O + 99)
        cnt = 0
        for _ in range(N_PERM):
            wp = W.copy()
            for s in np.unique(S):
                idx = np.where(S == s)[0]
                wp[idx] = rng.permutation(W[idx])
            cnt += stat(Y - tau0 * W + tau0 * wp, X, wp.astype(float)) >= s0
        q = X.shape[1]
        rows.append(dict(run=o + "_ldl_cc", covariate_set=name, n=len(Y), n_terms=q,
                         p_F=stats.f.sf(s0 / q, q, len(Y) - 2 * q - 2), p_rand=(1 + cnt) / (1 + N_PERM)))
        print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(HERE / "c_ldl_interaction.csv", index=False)
