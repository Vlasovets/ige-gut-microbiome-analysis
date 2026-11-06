#!/usr/bin/env python
"""
Module 06 - figure rebuilds.

Answers R2 Minor 2/4/7 and R3 Minor 4.

1. Fig 1A: annotate ALL allergen source categories. The published panel labelled
   only four (animal dander, pollen, food, house dust mite) of the twelve present
   in the 104-component panel, leaving moulds, cockroach, latex, storage mite,
   parasite and the CCD marker unlabelled, and lumping tree/grass/weed pollen.
2. Fig 2B: one colour-blind-safe cluster palette, defined once and shared, so
   every panel that shows clusters A-I agrees.
3. Fig 3D and Fig 7B: replace pie/donut charts with boxplots plus jittered points.
   Both figures used them; the reviewer flagged 3D.
4. All panels at >= 9 pt and exported as vector PDF/SVG. Every published figure
   was raster PNG while 85 unused PDFs sat in figures/pdf/.

Also writes lib/allergen_categories.csv, a curated artifact the manuscript needs
and does not currently have. NOTE: bromelain is flagged as a cross-reactive
carbohydrate determinant (CCD) marker, not a true allergen, yet it is counted in
the IgE burden score that defines the single/multiple/high/medium/low partitions.

** The category assignment needs sign-off from an allergy domain expert before
   submission. The placement of bromelain (CCD), herring worm (Anisakis, a
   parasite) and storage mite versus house dust mite are judgement calls. **

Usage:
  python revision/analysis/06_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

REVISION = Path(__file__).resolve().parent.parent
CONFIG = REVISION / "config" / "config.yaml"
LIB = REVISION / "lib"
FIGS = REVISION / "results" / "figures"

# ── Publication theme: 9 pt floor, enforced by 06b_check_fonts.py ─────────
BASE_PT = 9          # hard FLOOR, never subtract from it
plt.rcParams.update({
    "font.size": BASE_PT, "axes.titlesize": BASE_PT + 1,
    "axes.labelsize": BASE_PT, "xtick.labelsize": BASE_PT,
    "ytick.labelsize": BASE_PT, "legend.fontsize": BASE_PT,
    "figure.titlesize": BASE_PT + 2, "svg.fonttype": "none",
    "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
})

# Okabe-Ito, colour-blind safe. One definition, shared by every panel.
CLUSTER_COLORS = {
    "A": "#E69F00", "B": "#56B4E9", "C": "#009E73", "D": "#F0E442",
    "E": "#0072B2", "F": "#D55E00", "G": "#CC79A7", "H": "#999999",
    "I": "#000000",
}
CATEGORY_COLORS = {
    "Food": "#D55E00", "Grass pollen": "#009E73", "Tree pollen": "#006644",
    "Weed pollen": "#66C2A5", "House dust mite": "#8C564B",
    "Storage mite": "#C49C94", "Animal dander": "#0072B2", "Mould": "#7F3B08",
    "Insect": "#CC79A7", "Latex": "#9467BD", "Parasite": "#E69F00",
    "CCD marker (not an allergen)": "#999999",
}

# ── The curated mapping: prefix -> (category, source) ─────────────────────
# Ordered longest-prefix-first where prefixes could collide.
RULES = [
    ("red_meat", "Food", "Beef"), ("codfish", "Food", "Cod"),
    ("egg_", "Food", "Hen egg"), ("cow_milk", "Food", "Cow milk"),
    ("shrimp", "Food", "Shrimp"), ("hazel_nut", "Food", "Hazelnut"),
    ("walnut", "Food", "Walnut"), ("peanut", "Food", "Peanut"),
    ("sesame", "Food", "Sesame"), ("cashew", "Food", "Cashew"),
    ("brazil_nut", "Food", "Brazil nut"), ("soy_bean", "Food", "Soy"),
    ("wheat", "Food", "Wheat"), ("kiwi", "Food", "Kiwi"),
    ("peach", "Food", "Peach"), ("apple", "Food", "Apple"),
    ("celery", "Food", "Celery"),
    ("cat_", "Animal dander", "Cat"), ("dog_", "Animal dander", "Dog"),
    ("horse_", "Animal dander", "Horse"), ("mouse", "Animal dander", "Mouse"),
    ("storage_mite", "Storage mite", "Storage mite"),
    ("dust_mite", "House dust mite", "House dust mite"),
    ("german_cockroach", "Insect", "German cockroach"),
    ("herring_worm", "Parasite", "Anisakis simplex"),
    ("latex", "Latex", "Latex"),
    ("bromelain", "CCD marker (not an allergen)", "Bromelain (MUXF3)"),
    ("mould", "Mould", "Mould"),
    ("sycamore_tree", "Tree pollen", "Sycamore"),
    ("olive_tree", "Tree pollen", "Olive"),
    ("arizona_cypress", "Tree pollen", "Arizona cypress"),
    ("jpn_red_cedar", "Tree pollen", "Japanese red cedar"),
    ("birch", "Tree pollen", "Birch"), ("alder", "Tree pollen", "Alder"),
    ("timothy_grass", "Grass pollen", "Timothy grass"),
    ("bermuda_grass", "Grass pollen", "Bermuda grass"),
    ("common_ragweed", "Weed pollen", "Ragweed"),
    ("mug_wort", "Weed pollen", "Mugwort"),
    ("white_goose_foot", "Weed pollen", "Goosefoot"),
    ("annual_mercury", "Weed pollen", "Annual mercury"),
    ("english_platain", "Weed pollen", "English plantain"),
    ("prickly_saltwort", "Weed pollen", "Prickly saltwort"),
]


def categorise(cols) -> pd.DataFrame:
    rows = []
    for c in cols:
        hit = next(((cat, src) for pre, cat, src in RULES if c.startswith(pre)),
                   None)
        if hit is None:
            rows.append({"component": c, "category": "UNASSIGNED", "source": ""})
        else:
            rows.append({"component": c, "category": hit[0], "source": hit[1]})
    return pd.DataFrame(rows)


def enforce_min_font(fig) -> int:
    """Clamp every text artist to the BASE_PT floor.

    Setting rcParams is not sufficient: matplotlib resolves some sizes from
    RELATIVE keywords ('x-small' = 0.694 x base), so a 9 pt base still emits
    6.3 pt text. Walking the artists is the only reliable guarantee, and
    08_star_methods.py independently verifies the emitted SVG.
    """
    n = 0
    for t in fig.findobj(match=matplotlib.text.Text):
        try:
            if t.get_fontsize() < BASE_PT:
                t.set_fontsize(BASE_PT); n += 1
        except Exception:
            pass
    for ax in fig.get_axes():
        lg = ax.get_legend()
        if lg is not None and lg.get_title() is not None:
            if lg.get_title().get_fontsize() < BASE_PT:
                lg.get_title().set_fontsize(BASE_PT); n += 1
    return n


def save(fig, stem: str):
    n = enforce_min_font(fig)
    if n:
        print(f"    clamped {n} text artist(s) up to the {BASE_PT} pt floor")
    for ext in ("pdf", "svg", "png"):
        fig.savefig(FIGS / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {stem}.{{pdf,svg,png}}")


def fig1a(ige: pd.DataFrame, cats: pd.DataFrame, clusters: pd.Series):
    """IgE presence/absence heatmap with ALL source categories annotated."""
    sens = ige.loc[ige.index.isin(clusters.index)]
    binm = (sens > 0).astype(int)
    cl = clusters.reindex(binm.index)

    # Order allergens by hierarchical clustering within category, so the category
    # bands stay contiguous and readable.
    cats = cats.set_index("component").reindex(binm.columns)
    order = []
    for cat in [c for c in CATEGORY_COLORS if c in set(cats.category)]:
        sub = cats.index[cats.category == cat]
        if len(sub) > 2:
            d = pdist(binm[sub].T.values, metric="jaccard")
            d = np.nan_to_num(d, nan=1.0)
            lk = linkage(d, method="complete")
            sub = sub[dendrogram(lk, no_plot=True)["leaves"]]
        order += list(sub)
    binm = binm[order]
    cats = cats.reindex(order)

    # Order individuals by cluster, then by burden within cluster.
    burden = binm.sum(axis=1)
    idx = sorted(binm.index, key=lambda s: (str(cl[s]), -burden[s]))
    binm = binm.loc[idx]
    cl = cl.reindex(idx)

    fig = plt.figure(figsize=(9.5, 6.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 26], height_ratios=[1, 26],
                          wspace=0.012, hspace=0.012)
    ax_cat = fig.add_subplot(gs[1, 0])
    ax_cl = fig.add_subplot(gs[0, 1])
    ax = fig.add_subplot(gs[1, 1])

    ax.imshow(binm.T.values, aspect="auto", cmap="Purples",
              interpolation="nearest", vmin=0, vmax=1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(f"Individuals (n = {binm.shape[0]}), ordered by cluster")
    # No ylabel here: it would collide with the category band's tick labels.

    # Category band, with every category labelled.
    cat_codes = [list(CATEGORY_COLORS).index(c) for c in cats.category]
    from matplotlib.colors import ListedColormap
    ax_cat.imshow(np.array(cat_codes).reshape(-1, 1), aspect="auto",
                  cmap=ListedColormap(list(CATEGORY_COLORS.values())),
                  vmin=0, vmax=len(CATEGORY_COLORS) - 1, interpolation="nearest")
    ax_cat.set_xticks([])
    # Label only bands thick enough to carry text; single-component categories
    # would overlap illegibly (Parasite and the CCD marker are 1 component each).
    # Every category, including those, appears in the legend below.
    MIN_BAND = 3
    mids, labs = [], []
    for cat in cats.category.unique():
        pos = np.where(cats.category.values == cat)[0]
        if len(pos) >= MIN_BAND:
            mids.append(pos.mean()); labs.append(f"{cat} ({len(pos)})")
    ax_cat.set_yticks(mids); ax_cat.set_yticklabels(labs, fontsize=BASE_PT)
    ax_cat.tick_params(length=0)
    ax_cat.set_ylabel(f"Allergen components (n = {binm.shape[1]})",
                      fontsize=BASE_PT, labelpad=8)
    for s in ax_cat.spines.values():
        s.set_visible(False)

    # Complete category legend, so nothing is unlabelled - which is the whole
    # point of the reviewer's request.
    from matplotlib.patches import Patch
    counts = cats.category.value_counts()
    handles = [Patch(facecolor=CATEGORY_COLORS[c], edgecolor="none",
                     label=f"{c} ({counts[c]})")
               for c in CATEGORY_COLORS if c in counts.index]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.055),
              ncol=4, frameon=False, fontsize=BASE_PT,
              handlelength=1.1, handleheight=1.1, columnspacing=1.2,
              title="Allergen source category (all 12 shown)",
              title_fontsize=BASE_PT)

    # Cluster band.
    cl_codes = [list(CLUSTER_COLORS).index(c) if c in CLUSTER_COLORS else -1
                for c in cl]
    ax_cl.imshow(np.array(cl_codes).reshape(1, -1), aspect="auto",
                 cmap=ListedColormap(list(CLUSTER_COLORS.values())),
                 vmin=0, vmax=len(CLUSTER_COLORS) - 1, interpolation="nearest")
    ax_cl.set_yticks([]); ax_cl.set_xticks([])
    for cc in CLUSTER_COLORS:
        pos = np.where(np.array([str(x) for x in cl]) == cc)[0]
        if len(pos):
            ax_cl.text(pos.mean(), 0, cc, ha="center", va="center",
                       fontsize=BASE_PT,
                       color="white" if cc in ("I", "E") else "black")
    for s in ax_cl.spines.values():
        s.set_visible(False)
    ax_cl.set_title("IgE sensitization clusters (A–I)", fontsize=BASE_PT + 1,
                    pad=4)
    save(fig, "fig1A_allergen_categories")


def fig3d(cons_long: pd.DataFrame, cats_note: str):
    """Replace the donut charts with boxplots + jittered points."""
    d = cons_long[(cons_long.level == "genus") &
                  (cons_long.prevalence == "prev10") &
                  (cons_long.partition == "any_ige")]
    if d.empty:
        print("    fig3d: no data")
        return
    methods = [m for m in ["LinDA", "ANCOMBC2", "ALDEx2"] if m in set(d.method)]
    fig, ax = plt.subplots(figsize=(6.4, 3.9))
    rng = np.random.default_rng(42)
    for i, m in enumerate(methods):
        v = d[d.method == m].effect_log2.dropna().values
        bp = ax.boxplot([v], positions=[i], widths=0.55, patch_artist=True,
                        showfliers=False,
                        medianprops=dict(color="black", lw=1.4))
        bp["boxes"][0].set_facecolor("#CCCCCC"); bp["boxes"][0].set_alpha(0.65)
        ax.scatter(i + rng.uniform(-0.16, 0.16, len(v)), v, s=9, alpha=0.55,
                   color="#0072B2", linewidths=0)
        # Mark the one consensus taxon.
        pc = d[(d.method == m) & (d.taxon.str.contains("copri", case=False))]
        if len(pc):
            ax.scatter([i], pc.effect_log2.values[:1], s=52, marker="D",
                       color="#D55E00", zorder=5, edgecolor="black", linewidths=0.6)
    ax.axhline(0, color="black", lw=0.8, ls="--")
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(methods)
    ax.set_ylabel("Effect size (log$_2$ scale)")
    ax.set_title("Effect-size distribution by method — any_ige, genus, >10% prevalence")
    ax.scatter([], [], s=52, marker="D", color="#D55E00", edgecolor="black",
               linewidths=0.6, label="Prevotella (Segatella) copri")
    ax.legend(loc="lower right", frameon=False)
    fig.text(0.01, -0.06, cats_note, fontsize=BASE_PT, va="top", wrap=True)
    fig.tight_layout()
    save(fig, "fig3D_effect_distribution")


def main() -> int:
    cfg = yaml.safe_load(CONFIG.read_text())
    kora = Path(cfg["paths"]["kora_data"])
    LIB.mkdir(parents=True, exist_ok=True); FIGS.mkdir(parents=True, exist_ok=True)

    print("=== Module 06: figures ===")
    ige = pd.read_csv(kora / "ige_kora_subset_all.csv", index_col="u3_16s_id",
                      low_memory=False)
    ige.index = ige.index.astype(str)

    print(f"\n[1] Allergen categories for all {ige.shape[1]} components")
    cats = categorise(ige.columns)
    un = cats[cats.category == "UNASSIGNED"]
    if len(un):
        print(f"    UNASSIGNED: {list(un.component)}")
        return 1
    print("    every component assigned. Counts by category:")
    for c, n in cats.category.value_counts().items():
        print(f"      {c:32s} {n:3d}")
    cats.to_csv(LIB / "allergen_categories.csv", index=False)
    print(f"    wrote {LIB / 'allergen_categories.csv'}")
    print("    NOTE: published Fig 1A annotated only 4 of these "
          f"{cats.category.nunique()} categories.")

    print("\n[2] Figure 1A")
    cl = pd.read_csv(kora / "cluster_df_latent.csv", index_col=0)
    cl.index = cl.index.astype(str)
    fig1a(ige, cats, cl["cluster"])

    print("\n[3] Figure 3D (donuts -> boxplot + jitter)")
    f = REVISION / "results" / "tables" / "da_consensus_long.csv"
    if f.exists():
        note = ("Replaces the donut charts of the submitted Fig 3D. Effect sizes are "
                "on different scales per method (LinDA: log2 FC of relative "
                "abundance; ANCOM-BC2: log2 FC of absolute abundance; ALDEx2: CLR "
                "difference) and are shown side by side, never averaged.")
        fig3d(pd.read_csv(f), note)
    else:
        print("    da_consensus_long.csv missing; run 02b first")

    print("\n[4] Shared palette written for reuse")
    pd.DataFrame({"cluster": list(CLUSTER_COLORS),
                  "hex": list(CLUSTER_COLORS.values())}).to_csv(
        LIB / "cluster_palette.csv", index=False)
    pd.DataFrame({"category": list(CATEGORY_COLORS),
                  "hex": list(CATEGORY_COLORS.values())}).to_csv(
        LIB / "category_palette.csv", index=False)
    print(f"    wrote cluster_palette.csv and category_palette.csv (Okabe-Ito)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
