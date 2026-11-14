"""Figure for the manuscript: program effect within Doyle's pre-specified subgroups,
HbA1c (primary) and weight, from outputs/v3/<run>/table_subgroups.csv (main runs).

    python code/manuscript_figure.py                 ->  manuscript/figure1.png (+ .pdf), P = .01
    python code/manuscript_figure.py --leading-zero  ->  manuscript/figure1_leading_zero.png (+ .pdf), P = 0.01
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript"
INK, MUTED, GRID, ACCENT = "#1a1a1a", "#6b6b6b", "#d9d9d9", "#1f5fa8"

PANELS = [("hba1c", "A. HbA1c, percentage points"), ("weight", "B. Body weight, kg")]
ORDER = ["Baseline HbA1c", "Site", "Sex", "Prior-year inpatient admission",
         "Prior-year ED visit", "Insulin prescribed", "Metformin prescribed"]
LABEL = {"Baseline HbA1c": "Baseline HbA1c", "Prior-year inpatient admission": "Inpatient stay, prior year",
         "Prior-year ED visit": "ED visit, prior year", "Insulin prescribed": "Insulin prescribed",
         "Metformin prescribed": "Metformin prescribed", "Site": "Site", "Sex": "Sex"}


def load(run):
    d = pd.read_csv(ROOT / "outputs" / "v3" / run / "table_subgroups.csv")
    return d[d.family == "Doyle pre-specified"]


LEADING_ZERO = "--leading-zero" in sys.argv[1:]   # zero before the decimal point (P = 0.01)
NO_X = 0.82 if LEADING_ZERO else 0.86     # x of the No. column (axes fraction); P labels are wider


def fmt_p(p):  # 2 decimals, 3 if P < .01; leading zero dropped unless --leading-zero
    lead = (lambda s: s) if LEADING_ZERO else (lambda s: s.lstrip("0"))
    if p < 0.001:
        return "< " + lead("0.001")
    if p > 0.99:
        return "> " + lead("0.99")
    return "= " + lead(f"{p:.3f}" if p < 0.01 else f"{p:.2f}")


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "font.size": 8})
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 4.4), sharey=True)
    for ax, (run, title) in zip(axes, PANELS):
        d = load(run)
        y, ticks, labels = 0, [], []
        for sg in ORDER:
            rows = d[d.subgroup == sg]
            holm = rows.p_holm.iloc[0]
            hi_light = holm < 0.05
            for _, r in rows.iterrows():
                col = ACCENT if hi_light else INK
                ax.plot([r.lo, r.hi], [y, y], color=col, lw=1.2, solid_capstyle="butt")
                ax.plot(r.est, y, "s", ms=4.5, color=col)
                ax.text(NO_X, y, f"{r.n}", transform=ax.get_yaxis_transform(), ha="right",
                        va="center", fontsize=7, color=MUTED)
                ticks.append(y); labels.append(f"{LABEL[sg]}: {r.level.replace('>=', '≥')}")
                y -= 1
            ax.text(1.0, y + 1.5, f"P {fmt_p(holm)}", transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=7, color=ACCENT if hi_light else MUTED)
            y -= 0.6
        ax.text(NO_X, ticks[0] + 0.9, "No.", transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=7, color=MUTED)
        ax.text(1.0, ticks[0] + 0.9, "Holm P", transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=7, color=MUTED)
        ax.axvline(0, color=MUTED, lw=0.8)
        ax.set_title(title, loc="left", fontsize=9, fontweight="bold", color=INK)
        ax.set_xlabel("Intervention effect at 6 months (95% CI)", color=INK)
        ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(ticks, labels, color=INK)
    for ax in axes:  # room on the right for the No. and Holm P columns
        lo, hi = ax.get_xlim()
        ax.set_xticks([t for t in ax.get_xticks() if lo <= t <= hi])
        ax.set_xlim(lo, hi + 0.55 * (hi - lo))
        ax.spines["bottom"].set_bounds(lo, hi)
        ax.xaxis.set_label_coords(0.5 / 1.55, -0.08)
    fig.tight_layout()
    stem = "figure1_leading_zero" if LEADING_ZERO else "figure1"
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=300)
    print("written", OUT / f"{stem}.png")


if __name__ == "__main__":
    main()
