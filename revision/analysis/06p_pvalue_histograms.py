#!/usr/bin/env python
"""
Module 06p - p-value histograms behind the Figure 4C volcanoes, and the
permuted-label calibration check (module 02i).

  Rows 1-2  real data: ASV level, 5% prevalence filter, full cohort adjusted,
            each group vs all 233 IgE-negative (the four volcano contrasts),
            LinDA (row 1) and ANCOM-BC2 (row 2). A calibrated test on data with
            little signal gives a near-flat histogram; a spike at 0 means
            either many true effects or anti-conservative p-values.
  Row 3     permuted IgE labels (20 permutations, primary contrast, ASV, 02i):
            pooled p-value histograms for LinDA, ANCOM-BC2 and ANCOM-BC (v1) -
            no true effect exists, so these MUST be flat if the method is
            calibrated - and, right, taxa called at q < 0.05 per permutation.

Per histogram: fraction p < 0.05 and Storey's pi0 (lambda = 0.5).

Runs in scont_plots, Roboto font.
Output: results/figures/review/figS_pvalue_histograms.{png,svg,pdf}
        results/tables/pvalue_calibration_summary.csv
Usage:  python revision/analysis/06p_pvalue_histograms.py
"""

from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
OUT = RES / "figures" / "review"
NAVY, RED = "#2a3f5f", "#d62728"
CONTRASTS = [("any_ige", "Any IgE"), ("high_score", "High score"),
             ("medium_score", "Medium score"), ("low_score", "Low score")]
MLAB = {"LinDA": "LinDA", "ANCOMBC2": "ANCOM-BC2", "ANCOMBC": "ANCOM-BC (v1)"}
BINS = np.linspace(0, 1, 21)


def pi0(p, lam=0.5):
    p = p[~np.isnan(p)]
    return min(1.0, (p > lam).mean() / (1 - lam))


def hist(ax, p, title, color):
    p = np.asarray(p, float); p = p[~np.isnan(p)]
    ax.hist(p, bins=BINS, color=color, ec="white", lw=0.4,
            weights=np.full(len(p), len(BINS) - 1) / len(p))   # density: flat = 1
    ax.axhline(1, color="#888888", lw=0.8, ls="--")
    ax.set_title(title, fontsize=9, color=NAVY)
    ax.text(0.97, 0.95, f"p<0.05: {100 * (p < 0.05).mean():.0f}%\n$\\pi_0$ = {pi0(p):.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "none", "alpha": 0.9})
    ax.set_xlim(0, 1)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.3)   # headroom for the text box
    ax.spines[["top", "right"]].set_visible(False)
    return {"n": len(p), "frac_p_lt_0.05": (p < 0.05).mean(), "pi0": pi0(p)}


cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "savefig.transparent": False,
                     "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "font.size": 8.5, "svg.fonttype": "none"})

fig, axes = plt.subplots(3, 4, figsize=(11, 7.6), gridspec_kw={"hspace": 0.6, "wspace": 0.3})
rows = []
for r, (m, col) in enumerate([("LinDA", NAVY), ("ANCOMBC2", RED)]):
    for c, (part, lab) in enumerate(CONTRASTS):
        d = pd.read_csv(RES / "da_full" / f"da_long_asv_prev05_{part}.csv")
        p = d.loc[d.method == m, "p_raw"].values
        s = hist(axes[r, c], p, f"{MLAB[m]}: {lab}", col)
        rows.append({"data": "real", "method": m, "contrast": part, **s})
    axes[r, 0].set_ylabel("density (flat = 1)")

cal = sorted((RES / "calibration").glob("calib_asv_perm*.csv"))
perm = pd.concat([pd.read_csv(f) for f in cal], ignore_index=True)
perm = perm[perm.perm > 0]
n_perm = perm.perm.nunique()
for c, (m, col) in enumerate([("LinDA", NAVY), ("ANCOMBC2", RED), ("ANCOMBC", "#7f7f7f")]):
    s = hist(axes[2, c], perm.loc[perm.method == m, "p_raw"].values,
             f"{MLAB[m]}: permuted labels", col)
    rows.append({"data": f"permuted ({n_perm})", "method": m, "contrast": "any_ige", **s})
axes[2, 0].set_ylabel("density (flat = 1)")
for ax in axes[-1, :3]:
    ax.set_xlabel("p (unadjusted)")

ax = axes[2, 3]
calls = (perm.assign(sig=perm.q_within < 0.05).groupby(["method", "perm"]).sig.sum()
             .unstack("method").reindex(columns=list(MLAB)))
rng = np.random.default_rng(1)
for i, m in enumerate(MLAB):
    v = calls[m].values
    ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=12,
               c={"LinDA": NAVY, "ANCOMBC2": RED, "ANCOMBC": "#7f7f7f"}[m], lw=0)
    ax.text(i, v.max() + 12, f"median {np.median(v):.0f}", ha="center", va="bottom", fontsize=7)
ax.set_xticks(range(3)); ax.set_xticklabels([MLAB[m] for m in MLAB], rotation=20)
ax.set_ylabel("ASVs with q < 0.05")
ax.set_ylim(-10, calls.values.max() * 1.25)
ax.set_title(f"calls per permutation (n = {n_perm})", fontsize=9, color=NAVY)
ax.spines[["top", "right"]].set_visible(False)
for ax, t in zip(axes[:, 0], ["A", "B", "C"]):
    ax.text(-0.3, 1.18, t, transform=ax.transAxes, fontsize=14, fontweight="bold")

OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"figS_pvalue_histograms.{ext}", dpi=300, bbox_inches="tight",
                transparent=False, facecolor="white")
summ = pd.DataFrame(rows)
summ.to_csv(RES / "tables" / "pvalue_calibration_summary.csv", index=False)
print(summ.round(3).to_string(index=False))
print("\ncalls per permutation (ASV, q<0.05):\n", calls.describe().round(1).to_string())
