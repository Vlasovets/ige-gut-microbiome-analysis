#!/usr/bin/env python
"""
Module 06g - Figure 1A rebuilt in the PUBLISHED ordering, with every allergen
source category annotated (R: "good to show all allergen categories, not only
pollen/food/mite/animal").

Supersedes the fig1A panel of 06_figures.py, which re-sorted IgEs by category
and so no longer matched the published clustering.

Ordering (recomputed, then asserted):
  IgEs         complete-Jaccard dendrogram of notebook cell 18, leaf 0 at bottom
  individuals  complete-Euclidean LAC dendrogram of cell 21, with the branch
               flips in lib/fig1A_lac_flips.csv (from 06f, which proves the
               result matches the published heatmap cell for cell)

Annotation (rows):
  Source          12 categories, lib/allergen_categories.csv (from 06_figures.py)
  Protein family  lib/allergen_components.csv; cross-reactive families coloured,
                  species-specific/other in grey
  Brackets        contiguous blocks of the published row order, each checked
                  against the two annotations before it is drawn

** Both annotation files need sign-off from an allergy domain expert. **

Output: results/figures/review/fig1A_all_categories.{png,svg,pdf}
Env:    scont_plots (cnsplots)
"""

import sys
from pathlib import Path

import cnsplots as cns
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig1a_order import apply_flips, leaf_ids, recompute  # noqa: E402

REV = Path(__file__).resolve().parent.parent
LIB = REV / "lib"
OUT = REV / "results" / "figures" / "review"
OUT.mkdir(parents=True, exist_ok=True)

cns.setup_matplotlib()
# cnsplots sets savefig.transparent=True (renders black in dark-mode viewers)
BASE_PT = 9
plt.rcParams.update({"savefig.transparent": False, "savefig.facecolor": "white",
                     "figure.facecolor": "white", "axes.facecolor": "white",
                     "font.size": BASE_PT, "legend.fontsize": BASE_PT,
                     "svg.fonttype": "none"})

PURPLE, YELLOW = "#440154", "#FDE725"        # viridis ends, as published
SEP = "#E8505B"                               # dotted cluster separators, as published
NAVY = "#2a3f5f"
# Cluster colours of published Fig 1D, read off its legend in the published PNG.
# Not notebook cell 44, which swaps G (brown there) and I (grey there).
CLUSTER_COLORS = {
    "A": "#0072B2", "B": "#E69F00", "C": "#CC79A7", "D": "#D55E00", "E": "#F0E442",
    "F": "#009E73", "G": "#999999", "H": "#56B4E9", "I": "#A6761D",
}
CATEGORY_COLORS = dict(pd.read_csv(LIB / "category_palette.csv").values)
FAMILY_COLORS = {
    "PR-10": "#1B9E77", "Profilin": "#D95F02", "Polcalcin": "#7570B3",
    "nsLTP": "#E7298A", "Seed storage protein": "#66A61E", "Lipocalin": "#E6AB02",
    "Serum albumin": "#A6761D", "Tropomyosin": "#1F78B4", "Parvalbumin": "#A6CEE3",
    "Carbohydrate determinant": "#000000", "Other / species-specific": "#DDDDDD",
}

# Left-hand brackets: (label, first row, last row, check). Rows are named, not
# indexed, and each block is verified to be contiguous and to hold what the
# label says. Wording may change after expert review.
BLOCKS = [
    ("Seed storage\nproteins", "peanut_2s_albumin_saver", "brazil_nut_2s_albumin",
     ("family", {"Seed storage protein"})),
    ("Animal dander", "dog_f5", "cat_d2", ("category", {"Animal dander"})),
    ("nsLTPs", "walnut_lipid", "hazel_nut_lipid_", ("family", {"nsLTP"})),
    ("Cypress, cedar\n+ CCD", "jpn_red_cedar", "bromelain",
     ("category", {"Tree pollen", "CCD marker (not an allergen)"})),
    ("Grass pollen", "bermuda_grass", "timothy_grass_p2", ("category", {"Grass pollen"})),
    ("PR-10: birch, alder\n+ plant foods", "soy_bean_pr10", "kiwi_pr10",
     ("family", {"PR-10"})),
    ("Mites", "dust_mite_pter_p1", "storage_mite",
     ("category", {"House dust mite", "Storage mite"})),
    ("Profilins", "birch_v2", "timothy_grass_p12", ("family", {"Profilin"})),
]


def enforce_min_font(fig):
    """Same rule as 06_figures.py: no text below the BASE_PT floor."""
    for t in fig.findobj(match=matplotlib.text.Text):
        if t.get_text() and t.get_fontsize() < BASE_PT:
            t.set_fontsize(BASE_PT)


def strip(ax, values, colors):
    cats = list(colors)
    codes = np.array([cats.index(v) for v in values]).reshape(-1, 1)
    ax.imshow(codes, aspect="auto", interpolation="nearest", extent=(0, 1, 0, len(values)),
              cmap=ListedColormap(list(colors.values())), vmin=0, vmax=len(cats) - 1)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def main():
    r = recompute()
    vis, lac, cluster = r["vis"], r["lac"], r["cluster"]
    flips = pd.read_csv(LIB / "fig1A_lac_flips.csv").linkage_row.tolist()
    lac_link = apply_flips(r["lac_link"], flips)
    cols = leaf_ids(lac_link, lac.index)                 # individuals, left -> right
    rows = list(vis.columns)[::-1]                       # IgEs, top -> bottom
    n_ind, n_ige = len(cols), len(rows)

    runs = "".join(pd.Series(cluster.reindex(cols).values).loc[lambda s: s != s.shift()])
    assert runs == "ABCDEFGHI", runs
    assert leaf_ids(r["ige_link"], r["ige_names"]).tolist() == rows[::-1]

    ann = (pd.read_csv(LIB / "allergen_components.csv").set_index("component")
           .join(pd.read_csv(LIB / "allergen_categories.csv").set_index("component")))
    assert set(ann.index) == set(rows) and ann.category.notna().all()
    assert set(ann.family_group) <= set(FAMILY_COLORS)
    ann = ann.reindex(rows)
    for lab, a, b, (field, allowed) in BLOCKS:
        blk = ann.iloc[rows.index(a):rows.index(b) + 1]
        assert len(blk) >= 3 and set(blk[field if field == "category" else "family_group"]) \
            <= allowed, (lab, blk[[field if field == "category" else "family_group"]])

    M = vis.loc[cols, rows].values.T                     # n_ige x n_ind, row 0 = top

    # ── layout, in inches ────────────────────────────────────────────────
    W, H = 7.2, 7.6
    fig = plt.figure(figsize=(W, H))
    def ax_in(x, y, w, h):
        return fig.add_axes([x / W, y / H, w / W, h / H])
    hx, hw = 1.62, 4.55                  # heatmap left, width
    hy, hh = 2.55, 3.9                   # heatmap bottom, height
    s1x, s2x, sw = hx - 0.33, hx - 0.17, 0.13

    ax = ax_in(hx, hy, hw, hh)
    ax.imshow(M, aspect="auto", interpolation="nearest", extent=(0, n_ind, 0, n_ige),
              cmap=ListedColormap([PURPLE, YELLOW]), vmin=0, vmax=1)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(f"Individuals (n = {n_ind})", color=NAVY, fontsize=BASE_PT + 2, pad=44,
                 fontweight="normal")

    # cluster separators + letters, as published
    lab = cluster.reindex(cols).values
    edges = [0] + [i for i in range(1, n_ind) if lab[i] != lab[i - 1]] + [n_ind]
    for e in edges[1:-1]:
        ax.axvline(e, color=SEP, ls=(0, (1, 1.5)), lw=0.8)
    axl = ax_in(hx, hy - 0.30, hw, 0.28)
    axl.set_xlim(0, n_ind); axl.set_ylim(0, 1); axl.axis("off")
    for a, b in zip(edges[:-1], edges[1:]):
        c = lab[a]
        axl.add_patch(plt.Rectangle((a, 0.62), b - a, 0.38, color=CLUSTER_COLORS[c], lw=0))
        axl.text((a + b) / 2, 0.22, c, ha="center", va="center", color=NAVY)

    # dendrograms; branch colour = cluster when the whole branch is one cluster
    leaf_cl = cluster.reindex(lac.index).values
    n = len(leaf_cl)
    members = {}
    for i, (a, b, *_ ) in enumerate(lac_link):
        members[n + i] = members.get(int(a), {int(a)}) | members.get(int(b), {int(b)})
    def link_color(k):
        cs = {leaf_cl[j] for j in members[k]}
        return CLUSTER_COLORS[cs.pop()] if len(cs) == 1 else "#9AA5B4"
    axd = ax_in(hx, hy + hh + 0.02, hw, 0.55)
    sch.dendrogram(lac_link, ax=axd, no_labels=True, link_color_func=link_color)
    axd.set_xlim(0, 10 * n_ind); axd.axis("off")
    for c in axd.collections:
        c.set_linewidth(0.5)

    axr = ax_in(hx + hw + 0.02, hy, 0.42, hh)
    sch.dendrogram(r["ige_link"], ax=axr, orientation="right", no_labels=True,
                   link_color_func=lambda k: "#9AA5B4")
    axr.set_ylim(0, 10 * n_ige); axr.axis("off")
    for c in axr.collections:
        c.set_linewidth(0.5)
    fig.text((hx + hw + 0.55) / W, (hy + hh / 2) / H, f"IgE (n = {n_ige})", rotation=90,
             ha="center", va="center", color=NAVY, fontsize=BASE_PT + 2)

    # annotation strips
    a1 = ax_in(s1x, hy, sw, hh); strip(a1, ann.category, CATEGORY_COLORS)
    a2 = ax_in(s2x, hy, sw, hh); strip(a2, ann.family_group, FAMILY_COLORS)
    for x, t in [(s1x, "Source"), (s2x, "Protein family")]:
        fig.text((x + sw / 2) / W, (hy + hh + 0.06) / H, t, rotation=90,
                 ha="center", va="bottom", color=NAVY)

    # brackets
    axb = ax_in(0.05, hy, s1x - 0.10, hh)
    axb.set_xlim(0, 1); axb.set_ylim(0, n_ige); axb.axis("off")
    for lab_txt, a, b, _ in BLOCKS:
        top = n_ige - rows.index(a); bot = n_ige - rows.index(b) - 1
        axb.plot([0.97, 0.97], [bot + 0.15, top - 0.15], color=NAVY, lw=1.2,
                 solid_capstyle="butt")
        axb.text(0.92, (top + bot) / 2, lab_txt, ha="right", va="center",
                 color=NAVY, linespacing=0.95)

    # LAC strip
    lac_y, lac_h = hy - 0.80, 0.42
    axc = ax_in(hx, lac_y, hw, lac_h)
    im = axc.imshow(lac.loc[cols].values.T, aspect="auto", interpolation="nearest",
                    cmap="viridis", vmin=0, vmax=3, extent=(0, n_ind, 3, 0))
    for e in edges[1:-1]:
        axc.axvline(e, color=SEP, ls=(0, (1, 1.5)), lw=0.8)
    axc.set_xticks([]); axc.set_yticks([0.5, 1.5, 2.5])
    axc.set_yticklabels(["LAC$_1$", "LAC$_2$", "LAC$_3$"], color=NAVY)
    axc.tick_params(length=0)
    for s in axc.spines.values():
        s.set_visible(False)
    axc.set_xlabel("Latent allergy components", color=NAVY, fontsize=BASE_PT + 2,
                   labelpad=4)
    cb = fig.colorbar(im, cax=ax_in(hx + hw + 0.10, lac_y, 0.08, lac_h), ticks=[0, 1, 2, 3])
    cb.outline.set_visible(False)
    cb.set_label("score", color=NAVY)

    # legends
    cnt_c = ann.category.value_counts(); cnt_f = ann.family_group.value_counts()
    fig.legend(handles=[Patch(color=YELLOW, label="IgE detected"),
                        Patch(color=PURPLE, label="not detected")],
               loc="upper left", bbox_to_anchor=((hx + hw + 0.05) / W, (hy + hh + 0.60) / H),
               frameon=False, handlelength=1.0, handleheight=1.0)
    short = {"CCD marker (not an allergen)": "CCD marker, MUXF3"}   # caption explains
    fig.legend(handles=[Patch(color=CATEGORY_COLORS[c], label=f"{short.get(c, c)} ({cnt_c[c]})")
                        for c in CATEGORY_COLORS if c in cnt_c],
               title="Source category (all 12)", loc="upper left", ncol=2,
               bbox_to_anchor=(0.03 / W, (lac_y - 0.32) / H), frameon=False,
               handlelength=1.0, handleheight=1.0, columnspacing=0.8, alignment="left")
    fig.legend(handles=[Patch(color=FAMILY_COLORS[c], label=f"{c} ({cnt_f[c]})")
                        for c in FAMILY_COLORS if c in cnt_f],
               title="Protein family", loc="upper left", ncol=2,
               bbox_to_anchor=(3.85 / W, (lac_y - 0.32) / H), frameon=False,
               handlelength=1.0, handleheight=1.0, columnspacing=0.8, alignment="left")

    enforce_min_font(fig)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"fig1A_all_categories.{ext}", dpi=300, bbox_inches="tight",
                    transparent=False, facecolor="white")
    print("wrote", OUT / "fig1A_all_categories.{png,svg,pdf}")


if __name__ == "__main__":
    main()
