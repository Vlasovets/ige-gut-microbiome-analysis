#!/usr/bin/env python
"""
Module 06o - Figure 4A: alpha diversity (Shannon, Faith's PD) per IgE group,
full cohort.

Each IgE group is shown with all its members; non-IgE = all 233 IgE-negative
participants. Colour = stratum (fig4_groups.py), with a gap between strata.
q above each box: BH across the 15 groups, from lm(index ~ group + sex + age +
bmi) vs all 233 controls (module 05d) - the same adjustment as the DA.

Runs in scont_plots (cnsplots), Roboto font.
Output: results/figures/review/fig4A_alpha.{png,svg,pdf}
Usage:  python revision/analysis/06o_fig4A_alpha.py
"""

from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

from fig4_groups import STRATA, layout

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
OUT = RES / "figures" / "review"
NAVY = "#2a3f5f"

ps = pd.read_csv(RES / "alpha" / "alpha_full_persample.csv", dtype={"sample": str})
base = pd.read_csv(RES / "nested_models" / "designs" / "base.csv", dtype=str)
tests = pd.read_csv(RES / "alpha" / "alpha_full_tests.csv")
members = {p: set(g.sample_id) for p, g in base[base.role == "case"].groupby("partition")}
members["non_ige"] = set(base.sample_id[base.partition == "__controls__"])
assert len(members["non_ige"]) == 233 and len(members["any_ige"]) == 275

cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "savefig.transparent": False,
                     "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "font.size": 9, "axes.labelsize": 10,
                     "xtick.labelsize": 9, "ytick.labelsize": 9, "svg.fonttype": "none"})

L = layout()
fig, axes = plt.subplots(2, 1, figsize=(8.0, 5.6), sharex=True, gridspec_kw={"hspace": 0.12})
for ax, (col, ylab) in zip(axes, [("shannon", "Shannon"), ("pd", "Faith PD")]):
    data = [ps.loc[ps["sample"].isin(members[p]), col].values for p, _, _, _ in L]
    bp = ax.boxplot(data, positions=[x for *_, x in L], widths=0.68, patch_artist=True,
                    flierprops={"marker": "o", "markersize": 2.5, "markeredgewidth": 0},
                    medianprops={"color": "white", "lw": 1.4},
                    whiskerprops={"color": NAVY, "lw": 0.9}, capprops={"color": NAVY, "lw": 0.9})
    for patch, fl, (_, _, c, _) in zip(bp["boxes"], bp["fliers"], L):
        patch.set_facecolor(c); patch.set_edgecolor(NAVY); patch.set_linewidth(0.8)
        fl.set_markerfacecolor(c)
    q = tests[tests["index"] == col].set_index("partition").q_value
    top = max(max(v) for v in data)
    for p, _, _, x in L:
        if p in q:
            ax.text(x, top * 1.03, f"{q[p]:.2f}", ha="center", va="bottom", fontsize=7,
                    color="#555555")
    ax.set_ylim(top=top * 1.12)
    ax.set_ylabel(ylab)
    ax.yaxis.grid(True, color="#ececec", lw=0.8); ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
axes[1].set_xticks([x for *_, x in L])
axes[1].set_xticklabels([lab for _, lab, _, _ in L], rotation=45, ha="right",
                        rotation_mode="anchor")
for t, (_, _, c, _) in zip(axes[1].get_xticklabels(), L):
    t.set_color(c); t.set_fontweight("bold")
axes[1].tick_params(axis="x", length=0)
fig.legend(handles=[Patch(color=c, label=s) for s, c, _ in STRATA], loc="upper center",
           bbox_to_anchor=(0.5, 0.955), ncol=5, frameon=False, fontsize=8.5,
           handlelength=1.0, columnspacing=1.2)

OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"fig4A_alpha.{ext}", dpi=300, bbox_inches="tight",
                transparent=False, facecolor="white")
print("wrote", OUT / "fig4A_alpha.{png,svg,pdf}")
print(tests.round(3).to_string(index=False))
