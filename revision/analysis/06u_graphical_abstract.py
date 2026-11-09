#!/usr/bin/env python
"""
Module 06u - graphical abstract of the revised manuscript.

The submitted graphical abstract showed a volcano plot with many hits, TaxSEA
enrichment of folic acid / vitamin A producers and four family-level changes,
none of which hold in the revised analysis. This version shows:
  top    cohort -> IgE profiles (detected IgE, individuals ordered by cluster, IgEs by
         the Figure 1 dendrogram; fig1a_order.recompute) -> 3 latent allergy components
  middle differential abundance (LinDA, genus, 10% prevalence; no q < 0.05),
         prediction (nested CV AUC, genus), taxon set enrichment (no set),
         the most consistent single signal (Prevotella, direction in 5 DA methods)
  bottom take-home sentence (= last sentence of the Summary)

All data panels read saved result tables (no participant-level data):
  results/da_full/da_long_genus_prev10_any_ige.csv
  results/da_extended_full/daext_long_genus_prev10_any_ige.csv
  results/tables/ml_benchmark_all_analyses.csv
  results/tables/taxsea_recomputed_full.csv
Output : results/figures/review/graphical_abstract_revised.{png,svg,pdf}
         plus each data panel separately (ga_panel_*.{png,svg}) for re-assembly
Env    : scont_plots (cnsplots), Roboto
Usage  : python revision/analysis/06u_graphical_abstract.py
"""

from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle

from fig1a_order import recompute

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
OUT = RES / "figures" / "review"

NAVY, GREY, LIGHT = "#2a3f5f", "#BBBBBB", "#EEF1F6"
IGE, NONIGE = "#E64B6E", "#5B1A5E"
PREV = "#0072B2"
CLUSTER = ["#0072B2", "#E69F00", "#CC79A7", "#D55E00", "#F0E442", "#009E73", "#999999", "#56B4E9", "#A6761D"]

# ── data ────────────────────────────────────────────────────────────────
da = pd.concat([pd.read_csv(RES / "da_full" / "da_long_genus_prev10_any_ige.csv"),
                pd.read_csv(RES / "da_extended_full" / "daext_long_genus_prev10_any_ige.csv")])
linda = da[da.method == "LinDA"]
prev = da[da.taxon == "Prevotella"].set_index("method")
METHODS = ["LinDA", "ALDEx2", "ZicoSeq", "limma_voom", "ANCOMBC2"]
assert set(METHODS) <= set(prev.index), prev.index.tolist()
assert (prev.loc[METHODS, "effect_log2"] < 0).all()

ml = pd.read_csv(RES / "tables" / "ml_benchmark_all_analyses.csv")
ml = ml[(ml.partition == "any_ige") & (ml.level == "genus")].set_index("model")
MODELS = {"L1-logistic": "Logistic (L1)", "RandomForest": "Random forest", "XGBoost": "XGBoost",
          "classo-StabSel": "c-lasso"}

ts = pd.read_csv(RES / "tables" / "taxsea_recomputed_full.csv")
ts = ts[(ts.sub_database == "Metabolite_producers") & (ts.prevalence == "prev05") & (ts.rank_by == "stat")]

# ── style ───────────────────────────────────────────────────────────────
cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "font.size": 8,
                     "savefig.transparent": False, "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "svg.fonttype": "none"})


def panel_volcano(ax):
    ax.scatter(linda.effect_log2, -np.log10(linda.q_within), s=9, c=GREY, edgecolors="none")
    p = linda[linda.taxon == "Prevotella"].iloc[0]
    ax.scatter(p.effect_log2, -np.log10(p.q_within), s=26, c=PREV, zorder=3)
    ax.annotate("Prevotella", (p.effect_log2, -np.log10(p.q_within)), xytext=(5, 5),
                textcoords="offset points", ha="left", fontsize=7, fontstyle="italic", color=PREV)
    ax.axhline(-np.log10(0.05), color=IGE, ls="--", lw=0.9)
    ax.text(linda.effect_log2.max(), -np.log10(0.05) + 0.06, "q = 0.05",
            ha="right", fontsize=6.5, color=IGE)
    ax.set_ylim(0, 1.6)
    ax.set_xlabel(r"$\log_2$ fold change", fontsize=7)
    ax.set_ylabel(r"$-\log_{10}$ q", fontsize=7)
    ax.tick_params(labelsize=6.5)


def panel_auc(ax):
    y = np.arange(len(MODELS))[::-1]
    for yi, (k, lab) in zip(y, MODELS.items()):
        r = ml.loc[k]
        ax.plot([r.auc_ci_low, r.auc_ci_high], [yi, yi], color=NAVY, lw=1.2)
        ax.scatter(r.auc_pooled, yi, s=24, color=NAVY, zorder=3)
    ax.axvline(0.5, color=IGE, ls="--", lw=0.9)
    ax.text(0.495, -0.45, "chance", ha="right", va="bottom", fontsize=6.5, color=IGE)
    ax.set_yticks(y)
    ax.set_yticklabels(list(MODELS.values()), fontsize=6.5)
    ax.set_xlim(0.35, 0.75)
    ax.set_ylim(-0.6, len(MODELS) - 0.2)
    ax.set_xlabel("AUC (nested CV)", fontsize=7)
    ax.tick_params(axis="x", labelsize=6.5)


def panel_taxsea(ax):
    s = ts.sort_values("PValue")
    x = np.arange(len(s))
    ax.bar(x, -np.log10(s.PValue), color=GREY, width=0.75)
    for xi, name in zip(x, s.taxonSetName):
        if name in ("MiMeDB_producers_of_Folic_acid", "MiMeDB_producers_of_Vitamin_A"):
            ax.bar(xi, -np.log10(s.set_index("taxonSetName").loc[name, "PValue"]) + 0.02, color=PREV, width=0.75)
    ax.axhline(-np.log10(0.05), color=IGE, ls="--", lw=0.9)
    ax.text(len(s) - 0.5, -np.log10(0.05) + 0.06, "p = 0.05", ha="right", fontsize=6.5, color=IGE)
    pv = s.set_index("taxonSetName").PValue
    xs = {n: i for i, n in enumerate(s.taxonSetName)}
    for name, lab, dy in (("MiMeDB_producers_of_Vitamin_A", "vitamin A", 0.45), ("MiMeDB_producers_of_Folic_acid", "folic acid", 0.45)):
        x0, y0 = xs[name], -np.log10(pv[name]) + 0.03
        ax.annotate(lab, (x0, y0), xytext=(x0, y0 + dy), ha="center", fontsize=6.5, color=PREV,
                    arrowprops={"arrowstyle": "-", "color": PREV, "lw": 0.8})
    ax.set_ylim(0, 1.6)
    ax.set_xticks([])
    ax.set_xlabel(f"{len(s)} metabolite-producer sets", fontsize=7)
    ax.set_ylabel(r"$-\log_{10}$ p", fontsize=7)
    ax.tick_params(labelsize=6.5)


def panel_prevotella(ax):
    lab = {"LinDA": "LinDA", "ALDEx2": "ALDEx2", "ZicoSeq": "ZicoSeq", "limma_voom": "limma-voom",
           "ANCOMBC2": "ANCOM-BC2"}
    y = np.arange(len(METHODS))[::-1]
    ax.axvline(0, color="#777777", lw=0.8)
    ax.scatter(prev.loc[METHODS, "effect_log2"], y, s=26, color=PREV, zorder=3)
    for yi, m in zip(y, METHODS):
        ax.plot([0, prev.loc[m, "effect_log2"]], [yi, yi], color=PREV, lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels([lab[m] for m in METHODS], fontsize=6.5)
    ax.set_xlim(-2.1, 0.6)
    ax.set_ylim(-0.6, len(METHODS) - 0.4)
    ax.set_xlabel(r"$\log_2$ fold change (IgE vs. non-IgE)", fontsize=7)
    ax.tick_params(axis="x", labelsize=6.5)


def box(fig, x, y, w, h, fc=LIGHT):
    fig.patches.append(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.012",
                                      transform=fig.transFigure, fc=fc, ec="none", zorder=-1))


def arrow(fig, x0, y0, x1, y1):
    fig.patches.append(FancyArrowPatch((x0, y0), (x1, y1), transform=fig.transFigure, arrowstyle="-|>",
                                       mutation_scale=12, color=NAVY, lw=1.4))


def person(fig, cx, cy, h, color):
    fig.patches.append(Circle((cx, cy + 0.38 * h), 0.14 * h, transform=fig.transFigure, fc=color, ec="none"))
    fig.patches.append(FancyBboxPatch((cx - 0.15 * h, cy - 0.45 * h), 0.30 * h, 0.62 * h,
                                      boxstyle="round,pad=0,rounding_size=0.02", transform=fig.transFigure,
                                      fc=color, ec="none"))


# ── composite ───────────────────────────────────────────────────────────
fig = plt.figure(figsize=(6, 6))
T = dict(color=NAVY, ha="center", va="center")

# top row: cohort -> IgE profiles -> LACs
box(fig, 0.02, 0.70, 0.96, 0.28)
fig.text(0.17, 0.955, "KORA FF4 adults", fontsize=9, fontweight="bold", **T)
person(fig, 0.11, 0.835, 0.13, IGE)
person(fig, 0.23, 0.835, 0.13, NONIGE)
fig.text(0.11, 0.745, "IgE-sensitized\nn = 275", fontsize=7, **T)
fig.text(0.23, 0.745, "non-sensitized\nn = 233", fontsize=7, **T)
arrow(fig, 0.315, 0.84, 0.365, 0.84)

fig.text(0.51, 0.955, "Allergen-specific IgE profiles", fontsize=9, fontweight="bold", **T)
hm = fig.add_axes([0.38, 0.735, 0.26, 0.185])
ige = recompute()
vis = ige["vis"].astype(bool)
cl = ige["cluster"].reindex(vis.index)
cols = sorted(cl.unique())
assert vis.shape == (275, 104) and [str(c) for c in cols] == list("ABCDEFGHI"), (vis.shape, cols)
order = [i for c in cols for i in vis.index[cl == c][np.argsort(-vis.loc[cl == c].sum(1).values, kind="stable")]]
# vis columns are already in dendrogram leaf order; Figure 1 shows them bottom-up
mat = vis.loc[order].values.T[::-1]
grid = np.where(mat[..., None], [0.99, 0.91, 0.15], [0.27, 0.0, 0.33])
hm.imshow(grid, aspect="auto", interpolation="nearest")
nr = grid.shape[0]
sizes = [int((cl == c).sum()) for c in cols]
c0 = 0
for n, c in zip(sizes, CLUSTER):
    hm.add_patch(plt.Rectangle((c0 - 0.5, nr - 0.5), n, nr * 0.09, color=c, clip_on=False))
    c0 += n
hm.set_ylim(nr * 1.09 - 0.5, -0.5)
hm.axis("off")
fig.text(0.51, 0.712, "9 sensitization clusters", fontsize=7, **T)
arrow(fig, 0.655, 0.84, 0.705, 0.84)

fig.text(0.84, 0.955, "3 latent allergy\ncomponents", fontsize=9, fontweight="bold", **T)
for i, (lab, c) in enumerate([("pollen", "#009E73"), ("food", "#D55E00"), ("house dust mite", "#8C564B")]):
    yy = 0.865 - i * 0.055
    fig.patches.append(FancyBboxPatch((0.735, yy - 0.019), 0.21, 0.038, boxstyle="round,pad=0,rounding_size=0.018",
                                      transform=fig.transFigure, fc=c, ec="none"))
    fig.text(0.84, yy, lab, fontsize=8, color="white", fontweight="bold", ha="center", va="center")

# middle: microbiome comparisons
fig.text(0.5, 0.665, "Gut microbiome (16S rRNA) of sensitized vs. non-sensitized individuals",
         fontsize=9, fontweight="bold", **T)
PANELS = [
    ("Differential abundance", "no taxon with q < 0.05", panel_volcano, (0.02, 0.345)),
    ("Prediction", "close to chance", panel_auc, (0.51, 0.345)),
    ("Taxon set enrichment", "no metabolite-producer set enriched", panel_taxsea, (0.02, 0.055)),
    ("Most consistent signal", "Prevotella lower in all 5 methods,\nnot significant after FDR correction",
     panel_prevotella, (0.51, 0.055)),
]
for title, sub, fn, (x, y) in PANELS:
    box(fig, x, y, 0.47, 0.295)
    fig.text(x + 0.235, y + 0.275, title, fontsize=8.5, fontweight="bold", **T)
    fig.text(x + 0.235, y + 0.243 - (0.008 if "\n" in sub else 0), sub, fontsize=7,
             fontstyle="italic" if title.startswith("Most") else "normal", **T)
    left = 0.12 if fn in (panel_auc, panel_prevotella) else 0.075
    ax = fig.add_axes([x + left, y + 0.06, 0.47 - left - 0.03, 0.15 - (0.012 if "\n" in sub else 0)])
    fn(ax)

fig.text(0.5, 0.022, "IgE sensitization in adults is at most weakly linked to gut microbiome composition",
         fontsize=8.5, fontweight="bold", **T)

OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"graphical_abstract_revised.{ext}", dpi=300, transparent=False, facecolor="white")

# separate data panels for re-assembly in a figure editor
for name, fn, size in (("volcano", panel_volcano, (2.4, 1.7)), ("auc", panel_auc, (2.4, 1.7)),
                       ("taxsea", panel_taxsea, (2.4, 1.7)), ("prevotella", panel_prevotella, (2.4, 1.7))):
    f2, a2 = plt.subplots(figsize=size, layout="constrained")
    fn(a2)
    for ext in ("png", "svg"):
        f2.savefig(OUT / f"ga_panel_{name}.{ext}", dpi=300, transparent=False, facecolor="white")
    plt.close(f2)
print("wrote", OUT / "graphical_abstract_revised.{png,svg,pdf}", "and ga_panel_*.{png,svg}")
print("LinDA genus prev10: min q", round(linda.q_within.min(), 3), "| Prevotella effects",
      prev.loc[METHODS, "effect_log2"].round(2).to_dict())
