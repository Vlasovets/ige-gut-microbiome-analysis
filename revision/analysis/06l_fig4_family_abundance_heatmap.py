#!/usr/bin/env python
"""
Module 06l - Figure 4B: per-sample abundance of the highlighted families by
IgE group (fills the space under panel A; panel C = LinDA / ANCOM-BC2 volcanoes).

Samples: full cohort, the same 275 sensitized + 233 IgE-negative used by the
full-cohort DA (02f inputs, family level, 5% prevalence filter). The burden
groups high / medium / low are disjoint and together equal any_ige (asserted).

Values: CLR over all families retained by the filter (pseudocount 0.5), then
z-scored per family across all 508 samples and clipped to +/-2.5 for colour.
Left: one column per sample, grouped high | medium | low | non-IgE; within a
group samples are ordered by Prevotellaceae CLR (a display choice: it is the
family with a consistent direction across methods, and its bimodality is the
visible structure). Right: group mean z, value printed.
Unadjusted display of the data; the inference is the covariate-adjusted DA.

Runs in scont_plots (cnsplots), Roboto font.
Output: results/figures/review/fig4B_family_abundance.{png,svg,pdf}
        results/tables/fig4B_family_abundance_group_means.csv (source data)
Usage:  python revision/analysis/06l_fig4_family_abundance_heatmap.py
"""

from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib.colors import TwoSlopeNorm

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
OUT = RES / "figures" / "review"
SRC = RES / "full_cohort" / "inputs" / "family" / "prev05"
FAM = {"Lachnospiraceae": "#fcc853", "Bacteroidaceae": "#078a07",
       "Oscillospiraceae": "#1f77b4", "Veillonellaceae": "#ff7f0e",
       "Prevotellaceae": "#c2185b"}
# group colours of panel A (06d_alpha_boxplot.py)
GROUPS = [("high_score", "high", "#fb9a99"), ("medium_score", "medium", "#56b4e9"),
          ("low_score", "low", "#d55e00"), ("non_ige", "non-IgE", "#854e88")]
SORT_BY, PSEUDO, ZCLIP, GAP = "Prevotellaceae", 0.5, 2.5, 6
NAVY = "#2a3f5f"

# ── data ──────────────────────────────────────────────────────────────────
read = lambda p: pd.read_csv(p, index_col=0)
otu = pd.concat([read(SRC / "any_ige" / "source" / "otu_ige.csv"),
                 read(SRC / "any_ige" / "source" / "otu_non_ige.csv")])
group = pd.Series("non_ige", index=read(SRC / "any_ige" / "source" / "otu_non_ige.csv").index)
for part, _, _ in GROUPS[:3]:
    ids = read(SRC / part / "source" / "otu_ige.csv").index
    assert not group.index.intersection(ids).size, f"{part} overlaps"
    group = pd.concat([group, pd.Series(part, index=ids)])
assert otu.index.is_unique and set(group.index) == set(otu.index), "groups != any_ige + non-IgE"
assert len(otu) == 508 and all(f in otu.columns for f in FAM)

log = np.log(otu + PSEUDO)
clr = log.sub(log.mean(axis=1), axis=0)[list(FAM)]
z = (clr - clr.mean()) / clr.std()

order = []
for part, _, _ in GROUPS:
    ids = group.index[group == part]
    order.append(clr.loc[ids, SORT_BY].sort_values(ascending=False).index)
means = pd.DataFrame({lab: z.loc[ids].mean() for (_, lab, _), ids in zip(GROUPS, order)})
means.to_csv(RES / "tables" / "fig4B_family_abundance_group_means.csv")

# matrix with white gaps between groups
blocks, xticks = [], []
pos = 0
for k, ids in enumerate(order):
    blocks.append(z.loc[ids].T.values)
    xticks.append(pos + len(ids) / 2)
    pos += len(ids)
    if k < len(order) - 1:
        blocks.append(np.full((len(FAM), GAP), np.nan))
        pos += GAP
mat = np.hstack(blocks)

# ── style ─────────────────────────────────────────────────────────────────
cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "savefig.transparent": False,
                     "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "font.size": 9, "axes.labelsize": 9.5,
                     "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "svg.fonttype": "none"})

# ── figure ────────────────────────────────────────────────────────────────
cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("white")
norm = TwoSlopeNorm(vmin=-ZCLIP, vcenter=0, vmax=ZCLIP)
fig = plt.figure(figsize=(8.0, 3.0))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 0.26, 0.025], height_ratios=[0.07, 1],
                      wspace=0.06, hspace=0.04)
ax_bar = fig.add_subplot(gs[0, 0]); ax = fig.add_subplot(gs[1, 0], sharex=ax_bar)
ax_mbar = fig.add_subplot(gs[0, 1]); ax_m = fig.add_subplot(gs[1, 1], sharey=ax)
cax = fig.add_subplot(gs[1, 2])

ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
ax.set_yticks(range(len(FAM))); ax.set_yticklabels(list(FAM))
for t, c in zip(ax.get_yticklabels(), FAM.values()):
    t.set_color(c); t.set_fontweight("bold")
ax.set_xticks(xticks)
ax.set_xticklabels([f"{lab}\n(n = {len(ids)})" for (_, lab, _), ids in zip(GROUPS, order)])
ax.tick_params(axis="x", length=0)
ax.spines[:].set_visible(False)

pos = 0
for k, ((_, _, col), ids) in enumerate(zip(GROUPS, order)):
    ax_bar.add_patch(plt.Rectangle((pos - 0.5, 0), len(ids), 1, color=col, lw=0))
    pos += len(ids) + GAP
ax_bar.set_ylim(0, 1); ax_bar.axis("off")

# same colour scale as the per-sample matrix: group means are near 0, so the
# printed value carries the information
ax_m.imshow(means.values, aspect="auto", cmap=cmap, norm=norm)
for i in range(means.shape[0]):
    for j in range(means.shape[1]):
        ax_m.text(j, i, f"{means.iat[i, j]:+.2f}", ha="center", va="center", fontsize=7.5)
ax_m.tick_params(axis="both", length=0, labelleft=False, labelbottom=False)
ax_m.spines[:].set_visible(False)
for j, (_, _, col) in enumerate(GROUPS):
    ax_mbar.add_patch(plt.Rectangle((j - 0.5, 0), 1, 1, color=col, lw=0, ec="white"))
ax_mbar.set_xlim(-0.5, len(GROUPS) - 0.5); ax_mbar.set_ylim(0, 1); ax_mbar.axis("off")
ax_mbar.set_title("group mean", fontsize=8.5, color=NAVY, pad=3)

cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
cb.set_label("CLR abundance (z)", fontsize=8.5); cb.ax.tick_params(labelsize=8)

OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"fig4B_family_abundance.{ext}", dpi=300, bbox_inches="tight",
                transparent=False, facecolor="white")
print("wrote", OUT / "fig4B_family_abundance.{png,svg,pdf}")
print("group sizes:", {lab: len(ids) for (_, lab, _), ids in zip(GROUPS, order)})
print("group mean z:\n", means.round(3))
