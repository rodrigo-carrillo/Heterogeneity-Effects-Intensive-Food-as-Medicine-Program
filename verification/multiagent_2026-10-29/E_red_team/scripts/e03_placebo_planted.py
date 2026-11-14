"""E03: placebo (within-stratum permutation of treatment) and planted-effect simulations,
running the project's own cross-fitting core (code/hte_lib.crossfit_repeat) and the same
BLP / RATE / 7-subgroup Holm formulas as code/hte_food_as_medicine.py §7 and §9.

Fast settings (500 trees, 200 nuisance trees, 5 folds, 1 repeat, 200 bootstrap draws).
Each task runs in its own process (max_tasks_per_child=1) because EconML can segfault
in long-lived processes (see code/hte_lib.py).

Scenarios (all start from W permuted within the 4 randomization strata on the 465 randomized
patients, so the true effect of the permuted W is exactly zero):
  null          : nothing added
  pc05 / pc10   : pipeline positive control, +/- delta/2 by age >= median, delta = 0.5 / 1 SD of Y
  plant_ed      : weight only, treated with prior-year ED get -3.8 kg (the observed interaction)
  plant_lin     : HbA1c only, tau(x) = -0.5 * z(baseline HbA1c) pp (linear, SD 0.5 pp)

usage: python e03_placebo_planted.py <n_perm_null> <n_perm_sim> <workers>
"""
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def run_task(spec):
    import warnings
    warnings.filterwarnings("ignore")
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    from scipy.stats import norm
    from statsmodels.stats.multitest import multipletests
    from common import load_raw, analytic, ROOT, PRESPEC
    sys.path.insert(0, str(ROOT / "code"))
    import hte_lib

    outcome, perm, scen = spec["outcome"], spec["perm"], spec["scen"]
    raw = load_raw()
    rng = np.random.default_rng(1000 + perm)
    w = raw["treat_now"].to_numpy().copy()
    if perm >= 0:
        for s in np.unique(raw["stratum"]):
            idx = np.where(raw["stratum"].to_numpy() == s)[0]
            w[idx] = rng.permutation(w[idx])
    raw["treat_now"] = w
    df, F = analytic(raw, outcome)
    X = df[F].to_numpy(float)
    W = df["treat_now"].to_numpy(int)
    Y = df["Y"].to_numpy(float).copy()
    e = df["E_HAT"].to_numpy()
    strata = df["stratum"].to_numpy(int)
    sdy = Y.std(ddof=1)
    true_contrast = np.nan
    if scen in ("pc05", "pc10"):
        k = 0.5 if scen == "pc05" else 1.0
        older = (df["age"] >= df["age"].median()).to_numpy()
        delta = k * sdy
        Y = Y + W * np.where(older, -delta / 2, delta / 2)
        true_contrast = delta
    elif scen == "plant_ed":
        Y = Y + W * np.where(df["ed_i_0mo"].to_numpy() == 1, -3.8, 0.0)
        true_contrast = -3.8
    elif scen == "plant_lin":
        z = (df["hba1c_val_0mo"] - df["hba1c_val_0mo"].mean()) / df["hba1c_val_0mo"].std()
        Y = Y + W * (-0.5 * z.to_numpy())
        true_contrast = -0.5
    cfg = {"n_trees": 500, "n_trees_nuis": 200, "k_folds": 5, "min_leaf_cf": 10,
           "stratum_cols": [F.index("rural_site"), F.index("hba1c_cat_control")]}
    t0 = time.time()
    tau, gam = hte_lib.crossfit_repeat(Y, W, X, e, strata, 0, 20260621 + 7919 * (perm + 1), cfg)
    # BLP (same as §7 blp_test)
    Xb = sm.add_constant(tau - tau.mean())
    m = sm.OLS(gam, Xb).fit(cov_type="HC1")
    b2, se2 = m.params[1], m.bse[1]
    # RATE / AUTOC (same as §7)
    def toc(b, p):
        o = np.argsort(-p, kind="stable"); bb = b[o]
        return np.cumsum(bb) / np.arange(1, len(bb) + 1) - bb.mean()
    ben, pri = -gam, -tau
    autoc = toc(ben, pri).mean()
    rb = np.random.default_rng(perm + 5)
    boots = [toc(ben[i], pri[i]).mean() for i in (rb.integers(0, len(Y), len(Y)) for _ in range(200))]
    ase = np.std(boots, ddof=1)
    # 7 Doyle subgroups, AIPW interaction (same as §9)
    ps, ests = [], {}
    for v in PRESPEC:
        vv = df[v].to_numpy()
        r = {}
        for lev in (0, 1):
            g = gam[vv == lev]; r[lev] = (g.mean(), g.std(ddof=1) / np.sqrt(len(g)))
        d = r[1][0] - r[0][0]; dse = np.hypot(r[1][1], r[0][1])
        ps.append(2 * (1 - norm.cdf(abs(d / dse)))); ests[v] = d
    holm = multipletests(ps, method="holm")[1]
    out = dict(outcome=outcome, perm=perm, scen=scen, n=len(Y), true_contrast=true_contrast,
               blp_b2=b2, blp_se=se2, blp_p1=1 - norm.cdf(b2 / se2), autoc=autoc, autoc_se=ase,
               autoc_p1=1 - norm.cdf(autoc / ase), ate=gam.mean(), tau_sd=tau.std(),
               min_p_sub=min(ps), min_holm=holm.min(), argmin_sub=PRESPEC[int(np.argmin(ps))],
               p_ed=ps[PRESPEC.index("ed_i_0mo")], holm_ed=holm[PRESPEC.index("ed_i_0mo")],
               est_ed=ests["ed_i_0mo"], secs=time.time() - t0)
    # pc truth: AIPW contrast older vs younger
    if scen in ("pc05", "pc10"):
        older = (df["age"] >= df["age"].median()).to_numpy()
        out["est_contrast"] = gam[older].mean() - gam[~older].mean()
    return out


def main():
    n_null, n_sim, workers = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    specs = []
    for p in range(n_null):
        for oc in ("hba1c", "weight", "sbp", "ed_ip"):
            specs.append(dict(outcome=oc, perm=p, scen="null"))
    for p in range(n_sim):
        for oc in ("hba1c", "weight"):
            for sc in ("pc05", "pc10"):
                specs.append(dict(outcome=oc, perm=p, scen=sc))
        specs.append(dict(outcome="weight", perm=p, scen="plant_ed"))
        specs.append(dict(outcome="hba1c", perm=p, scen="plant_lin"))
    outf = HERE.parent / "e03_results.jsonl"
    done = set()
    if outf.exists():
        for line in outf.read_text().splitlines():
            r = json.loads(line); done.add((r["outcome"], r["perm"], r["scen"]))
    specs = [s for s in specs if (s["outcome"], s["perm"], s["scen"]) not in done]
    # interleave so partial results cover every scenario
    specs.sort(key=lambda s: s["perm"])
    print(f"{len(specs)} tasks to run with {workers} workers", flush=True)
    attempts = {}
    with ProcessPoolExecutor(max_workers=workers, max_tasks_per_child=1) as ex:
        futs = {ex.submit(run_task, s): s for s in specs}
        while futs:
            for f in as_completed(list(futs)):
                s = futs.pop(f)
                try:
                    r = f.result()
                    with open(outf, "a") as fh:
                        fh.write(json.dumps({k: (float(v) if isinstance(v, (np.floating,)) else v)
                                             for k, v in r.items()}) + "\n")
                except Exception:
                    key = (s["outcome"], s["perm"], s["scen"])
                    attempts[key] = attempts.get(key, 0) + 1
                    print("FAILED", key, traceback.format_exc()[-300:], flush=True)
                    if attempts[key] < 3:
                        futs[ex.submit(run_task, s)] = s
                break
    print("done", flush=True)


if __name__ == "__main__":
    main()
