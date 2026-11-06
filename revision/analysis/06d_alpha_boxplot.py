#!/usr/bin/env python
"""
Module 06d - revised Figure 3A: alpha diversity boxplots, Shannon and Faith's PD.

Layout follows the submitted Figure 3A (same group order, same per-group
colours, sampled from figures/png/alpha_diversity_clusters.png), but plots the
analysis actually used in the revision: each IgE group next to ITS OWN matched
controls (1:1 on sex, age +-5 y, BMI +-4 kg/m2). The submitted panel compared
every subgroup with all 233 controls without adjustment, which is what
Reviewer 3 asked to change.

Annotations are read from the result files, never recomputed here:
  Shannon   paired Wilcoxon p, raw counts  results/alpha/alpha_tests_long.csv
  Faith PD  paired Wilcoxon p              results/phylo/faithpd_<group>.csv
q = Benjamini-Hochberg across the 15 groups within each panel.

Plotting uses cnsplots (author's standard). Runs in the scont_plots env.
Output: results/figures/review/fig3A_alpha_matched.{png,svg,pdf}

Usage:  python revision/analysis/06d_alpha_boxplot.py
"""

from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
OUT = RES / "figures" / "review"
OUT.mkdir(parents=True, exist_ok=True)

GROUPS = [("any_ige", "any"), ("single", "single"), ("multiple", "multiple"),
          ("high_score", "high"), ("medium_score", "medium"), ("low_score", "low")] + \
         [(f"cluster_{c}", c) for c in "ABCDEFGHI"]
# Colours sampled from the submitted figure, so readers can map old to new.
COL = {"any": "#f74e6b", "single": "#b3de69", "multiple": "#fdb462", "high": "#fb9a99",
       "medium": "#56b4e9", "low": "#d55e00", "A": "#66a61e", "B": "#8dd3c7",
       "C": "#f0e442", "D": "#d55e00", "E": "#cc79a7", "F": "#e69f00",
       "G": "#0072b2", "H": "#10cdcd", "I": "#84e200"}
CTRL = "#854e88"            # the submitted "non-IgE" colour, now = matched controls
HUE = ["IgE-positive", "matched non-IgE"]

# ── data: membership + pairs from the Faith per-sample files ──────────────
shan = pd.read_csv(RES / "alpha" / "alpha_persample.csv", dtype={"sample": str})
shan = shan[shan.variant == "raw"].set_index("sample")["shannon"]
rows = []
for part, lab in GROUPS:
    d = pd.read_csv(RES / "phylo" / f"faithpd_persample_{part}.csv", dtype={"sample": str})
    d["shannon"] = d["sample"].map(shan)
    d["group"] = d["group"].map({"ige": HUE[0], "non_ige": HUE[1]})
    d["part"] = lab
    rows.append(d)
long = pd.concat(rows, ignore_index=True)
assert long.shannon.notna().all(), "Shannon missing for some matched samples"

# ── annotations from the result files ─────────────────────────────────────
at = pd.read_csv(RES / "alpha" / "alpha_tests_long.csv")
at = at[(at.analysis == "paired_wilcoxon") & (at.variant == "raw") & (at["index"] == "shannon")]
p_sh = [float(at.loc[at.partition == part, "p_value"].iloc[0]) for part, _ in GROUPS]
p_pd = [float(pd.read_csv(RES / "phylo" / f"faithpd_{part}.csv").p_value.iloc[0])
        for part, _ in GROUPS]
npair = [int(at.loc[at.partition == part, "n_pairs"].iloc[0]) for part, _ in GROUPS]
q_sh = multipletests(p_sh, method="fdr_bh")[1]
q_pd = multipletests(p_pd, method="fdr_bh")[1]

# ── plot ──────────────────────────────────────────────────────────────────
cns.setup_matplotlib()
# cnsplots sets savefig.transparent=True: transparent pixels render BLACK in
# dark-mode viewers. Force an opaque white background for figure and axes.
plt.rcParams.update({"savefig.transparent": False, "savefig.facecolor": "white",
                     "figure.facecolor": "white", "axes.facecolor": "white"})
plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "xtick.labelsize": 10,
                     "ytick.labelsize": 10, "svg.fonttype": "none"})
order = [lab for _, lab in GROUPS]
fig, axes = plt.subplots(2, 1, figsize=(11, 7.2), sharex=True)
for ax, col, ylab, qs, tag in [(axes[0], "shannon", "Shannon", q_sh, "A"),
                               (axes[1], "pd", "Faith PD", q_pd, "B")]:
    cns.boxplot(data=long, x="part", y=col, hue="group", order=order, hue_order=HUE,
                showoutliers=True, ax=ax)
    # recolour by box position: IgE box = group colour, control box = CTRL
    for patch in ax.patches:
        verts = patch.get_path().vertices
        xc = (verts[:, 0].min() + verts[:, 0].max()) / 2
        i = int(round(xc))
        if not 0 <= i < len(order):
            continue
        patch.set_facecolor(COL[order[i]] if xc < i else CTRL)
        patch.set_edgecolor("#2a3f5f"); patch.set_alpha(1)
    # lines: whiskers, caps and medians in the submitted figure's navy; outlier
    # dots take the colour of the box they belong to
    for ln in ax.lines:
        xs = np.asarray(ln.get_xdata(), dtype=float)
        if not len(xs):
            continue
        xc = float(np.nanmean(xs)); i = int(round(xc))
        if ln.get_marker() not in ("None", None, "", " ") and ln.get_linestyle() in ("None", ""):
            if 0 <= i < len(order):
                c = COL[order[i]] if xc < i else CTRL
                ln.set_markerfacecolor(c); ln.set_markeredgecolor(c)
        else:
            ln.set_color("#2a3f5f"); ln.set_linewidth(1.0)
    ax.set_xlabel(""); ax.set_ylabel(ylab)
    if ax.get_legend():
        ax.get_legend().remove()
    top = ax.get_ylim()[1]
    span = top - ax.get_ylim()[0]
    for i, q in enumerate(qs):
        ax.text(i, top + 0.02 * span, f"q={q:.2f}", ha="center", va="bottom", fontsize=9,
                color="#2a3f5f")
    ax.set_ylim(ax.get_ylim()[0], top + 0.10 * span)
    ax.yaxis.grid(True, color="#e5e5e5", lw=0.8); ax.set_axisbelow(True)
    ax.text(-0.06, 1.04, tag, transform=ax.transAxes, fontsize=16, fontweight="bold")
axes[1].set_xticks(range(len(order)))
axes[1].set_xticklabels([f"{lab}\n({n})" for lab, n in zip(order, npair)])
from matplotlib.patches import Patch
fig.legend(handles=[Patch(facecolor="#bbbbbb", edgecolor="#2a3f5f", label="IgE-positive (group colour)"),
                    Patch(facecolor=CTRL, edgecolor="#2a3f5f", label="matched non-IgE")],
           loc="upper right", ncol=2, frameon=False, fontsize=10)
fig.text(0.5, 0.005, "IgE group (number of matched pairs)", ha="center", fontsize=11)
fig.tight_layout(rect=(0, 0.02, 1, 0.96))
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"fig3A_alpha_matched.{ext}", dpi=300, bbox_inches="tight",
                transparent=False, facecolor="white")
print("wrote", *(OUT / f"fig3A_alpha_matched.{e}" for e in ("png", "svg", "pdf")), sep="\n  ")
print("min q Shannon = %.3f | min q Faith PD = %.3f" % (q_sh.min(), q_pd.min()))
