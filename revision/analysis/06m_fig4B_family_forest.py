#!/usr/bin/env python
"""
Module 06m - Figure 4B: effect sizes for the highlighted families, per IgE group
(fills the space under panel A; panel C = LinDA / ANCOM-BC2 volcanoes).

One row per family (the four families of the submitted paper + Prevotellaceae,
row labels in the volcano colours so they double as the panel C key). Within a
row one point per contrast vs non-IgE: any, high, medium, low score, in the
panel A group colours (06d_alpha_boxplot.py).
  x = LinDA log2 fold change, bar = 95% CI (estimate +/- 1.96 SE),
  right margin = BH q within that analysis (bold if q < 0.05).

Design: full cohort, each group vs all 233 IgE-negative, ~ IgE status + sex +
age + BMI (02f --design full), family level, 5% prevalence filter - the same
analyses as the volcano panels.

Runs in scont_plots (cnsplots), Roboto font.
Output: results/figures/review/fig4B_family_forest.{png,svg,pdf}
        results/tables/fig4B_family_forest.csv (source data)
Usage:  python revision/analysis/06m_fig4B_family_forest.py
"""

from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
OUT = RES / "figures" / "review"
LEVEL, PREV, METHOD, ALPHA, Z = "family", "prev05", "LinDA", 0.05, 1.96
FAM = {"Lachnospiraceae": "#fcc853", "Bacteroidaceae": "#078a07",
       "Oscillospiraceae": "#1f77b4", "Veillonellaceae": "#ff7f0e",
       "Prevotellaceae": "#c2185b"}
GROUPS = [("any_ige", "any", "#f74e6b"), ("high_score", "high", "#fb9a99"),
          ("medium_score", "medium", "#56b4e9"), ("low_score", "low", "#d55e00")]
DODGE = 0.17
NAVY = "#2a3f5f"

# ── data ──────────────────────────────────────────────────────────────────
rows = []
for part, lab, _ in GROUPS:
    d = pd.read_csv(RES / "da_full" / f"da_long_{LEVEL}_{PREV}_{part}.csv")
    d = d[(d.method == METHOD) & d.taxon.isin(FAM)]
    assert len(d) == len(FAM) and (d.n_ctrl == 233).all(), part
    rows.append(d.assign(group=lab))
d = pd.concat(rows, ignore_index=True)
d["lo"], d["hi"] = d.effect_log2 - Z * d.se, d.effect_log2 + Z * d.se
d[["taxon", "group", "n_ige", "n_ctrl", "effect_log2", "se", "lo", "hi", "p_raw", "q_within",
   "formula"]].to_csv(RES / "tables" / "fig4B_family_forest.csv", index=False)

# ── style ─────────────────────────────────────────────────────────────────
cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "savefig.transparent": False,
                     "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "font.size": 9, "axes.labelsize": 9.5,
                     "xtick.labelsize": 8.5, "ytick.labelsize": 9, "svg.fonttype": "none"})

# ── figure ────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7.0, 3.6))
for i, fam in enumerate(FAM):
    if i % 2 == 0:   # light banding separates the family rows
        ax.axhspan(i - 0.5, i + 0.5, color="#f4f5f7", lw=0, zorder=0)
    for k, (_, lab, col) in enumerate(GROUPS):
        r = d[(d.taxon == fam) & (d.group == lab)].iloc[0]
        y = i + (k - 1.5) * DODGE
        ax.plot([r.lo, r.hi], [y, y], color=col, lw=1.6, solid_capstyle="round", zorder=2)
        ax.scatter(r.effect_log2, y, s=30, color=col, ec="white", lw=0.6, zorder=3)
        ax.text(1.015, y, f"q = {r.q_within:.2f}", transform=ax.get_yaxis_transform(),
                va="center", fontsize=7, color="#555555",
                fontweight="bold" if r.q_within < ALPHA else "normal")
ax.axvline(0, color="#999999", lw=0.8, zorder=1)
ax.set_yticks(range(len(FAM))); ax.set_yticklabels(list(FAM))
for t, c in zip(ax.get_yticklabels(), FAM.values()):
    t.set_color(c); t.set_fontweight("bold")
ax.set_ylim(len(FAM) - 0.5, -0.5)
lim = np.ceil(max(d.lo.abs().max(), d.hi.abs().max()) * 10) / 10
ax.set_xlim(-lim, lim)
ax.set_xlabel("log$_2$ fold change vs non-IgE (LinDA, 95% CI)")
ax.tick_params(axis="y", length=0)
ax.spines[["top", "right", "left"]].set_visible(False)
ax.text(-lim * 0.97, len(FAM) - 0.5, "lower in sensitized", fontsize=7.5, color="#777777",
        va="bottom", ha="left")
ax.text(lim * 0.97, len(FAM) - 0.5, "higher in sensitized", fontsize=7.5, color="#777777",
        va="bottom", ha="right")

n = d.groupby("group").n_ige.first()
h = [Line2D([], [], color=col, marker="o", ms=5, lw=1.6, mec="white",
            label=f"{lab} (n = {n[lab]})") for _, lab, col in GROUPS]
ax.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, frameon=False,
          fontsize=8.5, handletextpad=0.4, columnspacing=1.4,
          title="IgE group vs non-IgE (n = 233)", title_fontsize=8.5)

OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"fig4B_family_forest.{ext}", dpi=300, bbox_inches="tight",
                transparent=False, facecolor="white")
print("wrote", OUT / "fig4B_family_forest.{png,svg,pdf}")
print(d.pivot_table(index="taxon", columns="group", values="q_within").round(3))
