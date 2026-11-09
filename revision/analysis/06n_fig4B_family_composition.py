#!/usr/bin/env python
"""
Module 06n - Figure 4B (alternative to 06m): mean gut community composition by
IgE group, with the highlighted families coloured.

One bar per IgE group, ordered and coloured by stratum with gaps (fig4_groups.py): the mean
over participants of per-sample relative abundance. The four families of the
submitted paper + Prevotellaceae keep the volcano colours; all other families
are pooled in grey. Percentages printed inside segments >= MIN_LABEL %.

Samples: full cohort, the same 275 sensitized + 233 IgE-negative as the
full-cohort DA (02f inputs, family level, 5% prevalence filter). Relative
abundance is computed over the families retained by that filter. Descriptive
and unadjusted; the inference is the covariate-adjusted DA (panel C, 06m).

Runs in scont_plots (cnsplots), Roboto font.
Output: results/figures/review/fig4B_family_composition.{png,svg,pdf}
        results/tables/fig4B_family_composition.csv (source data)
Usage:  python revision/analysis/06n_fig4B_family_composition.py
"""

from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

from fig4_groups import layout

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
OUT = RES / "figures" / "review"
SRC = RES / "full_cohort" / "inputs" / "family" / "prev05"
FAM = {"Lachnospiraceae": "#fcc853", "Bacteroidaceae": "#078a07",
       "Oscillospiraceae": "#1f77b4", "Veillonellaceae": "#ff7f0e",
       "Prevotellaceae": "#c2185b"}
OTHER = ("other families", "#c9ccd3")
# group order, stratum colours and gaps shared with panel A (fig4_groups.py)
L = layout()
GROUPS = [(p, lab, col) for p, lab, col, _ in L]
XPOS = [x for *_, x in L]
MIN_LABEL = 3.0

# ── data ──────────────────────────────────────────────────────────────────
read = lambda p: pd.read_csv(p, index_col=0)
ctrl = read(SRC / "any_ige" / "source" / "otu_non_ige.csv")
members = {"non_ige": ctrl}
for part, _, _ in GROUPS[:-1]:
    members[part] = read(SRC / part / "source" / "otu_ige.csv")
    assert members[part].index.isin(members.get("any_ige", members[part]).index).all(), part
n_of = lambda parts: sum(len(members[p]) for p in parts)
assert n_of(["high_score", "medium_score", "low_score"]) == len(members["any_ige"]) == 275
assert n_of(["single", "multiple"]) == 275 and len(ctrl) == 233

comp = {}
for part, lab, _ in GROUPS:
    otu = members[part]
    ra = otu.div(otu.sum(axis=1), axis=0) * 100
    mean = ra.mean()
    comp[f"{lab} (n = {len(otu)})"] = pd.concat([mean[list(FAM)],
                                                 pd.Series({OTHER[0]: mean.drop(list(FAM)).sum()})])
comp = pd.DataFrame(comp).T
assert (comp.sum(axis=1).sub(100).abs() < 1e-6).all()
comp.to_csv(RES / "tables" / "fig4B_family_composition.csv")

# ── style ─────────────────────────────────────────────────────────────────
cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "savefig.transparent": False,
                     "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "font.size": 9, "axes.labelsize": 9.5,
                     "xtick.labelsize": 8.5, "ytick.labelsize": 9, "svg.fonttype": "none"})

# ── figure ────────────────────────────────────────────────────────────────
colours = list(FAM.values()) + [OTHER[1]]
fig, ax = plt.subplots(figsize=(8.0, 3.6))
for i, (grp, row) in zip(XPOS, comp.iterrows()):
    bottom = 0.0
    for (seg, val), col in zip(row.items(), colours):
        ax.bar(i, val, bottom=bottom, color=col, ec="white", lw=0.6, width=0.78)
        if val >= MIN_LABEL:
            ax.text(i, bottom + val / 2, f"{val:.0f}", ha="center", va="center", fontsize=6.5,
                    color="white" if seg in ("Bacteroidaceae", "Oscillospiraceae", "Prevotellaceae")
                    else "#222222")
        bottom += val
ax.set_xticks(XPOS)
# names only, as in panel A; group sizes are in the source-data table
ax.set_xticklabels([g[1] for g in GROUPS], rotation=45, ha="right", rotation_mode="anchor")
for t, (_, _, col) in zip(ax.get_xticklabels(), GROUPS):
    t.set_color(col); t.set_fontweight("bold")
ax.set_xlim(XPOS[0] - 0.5, XPOS[-1] + 0.5)
ax.set_ylim(0, 100)
ax.set_ylabel("mean relative abundance (%)")
ax.tick_params(axis="x", length=0)
ax.spines[["top", "right", "bottom"]].set_visible(False)
h = [Patch(color=c, label=f) for f, c in FAM.items()] + [Patch(color=OTHER[1], label=OTHER[0])]
ax.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=6, frameon=False,
          fontsize=8, handlelength=1.0, columnspacing=1.0, handletextpad=0.4)

OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"fig4B_family_composition.{ext}", dpi=300, bbox_inches="tight",
                transparent=False, facecolor="white")
print("wrote", OUT / "fig4B_family_composition.{png,svg,pdf}")
print(comp.round(2).to_string())
