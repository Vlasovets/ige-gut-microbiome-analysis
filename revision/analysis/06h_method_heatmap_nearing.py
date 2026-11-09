#!/usr/bin/env python
"""
Module 06h - DA method heatmap in the layout of Nearing et al. 2022, Figure 1
(Nat Commun 13:342).

Nearing's design, reproduced here:
  * one row per "dataset" - here one analysis = IgE group x taxonomic rank
    (15 groups x 4 ranks = 60 rows)
  * one column per DA method
  * cell text = NUMBER of taxa called significant (q_within < 0.05, BH within
    that analysis)
  * cell colour = the same NUMBER in fixed, discrete buckets
    (0, 1, 2, 3-5, 6-10, 11-50, 51-100, >100), so equal counts always get equal
    colours. DEVIATION from Nearing (row-standardized %): a co-author found the
    standardized colours misleading (zeros coloured differently by row)
  * left-hand annotation columns = analysis characteristics that could drive
    the variation (darker = higher): sensitized group size, features tested
  * DEVIATION from Nearing: rows are ordered by taxonomic rank (family, genus,
    species, ASV) and, within rank, by IgE group in a fixed order, not
    hierarchically clustered - so the same analysis sits on the same row in
    every panel and rank effects are readable
  * three panels, one per prevalence filter: a = 10%, b = 5%, c = 1%
    (Nearing: unfiltered vs 10%; no unfiltered table was run here)

Data: the FULL-COHORT, covariate-adjusted design (all sensitized vs all 233
controls, ~ IgE status + sex + age + BMI; 02f + --design full), five methods.

Runs in scont_plots (cnsplots). Output:
  results/figures/review/figS_da_methods_nearing.{png,svg,pdf}
  results/tables/da_full_method_counts_by_analysis.csv   (source data, counts only)

Usage:  python revision/analysis/06h_method_heatmap_nearing.py
"""

from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib.colors import BoundaryNorm, ListedColormap

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
OUT = RES / "figures" / "review"
OUT.mkdir(parents=True, exist_ok=True)
ALPHA = 0.05
METHODS = ["ALDEx2", "ZicoSeq", "LinDA", "limma_voom", "ANCOMBC2"]   # fewest -> most calls
MLAB = {"LinDA": "LinDA", "ANCOMBC2": "ANCOM-BC2", "ALDEx2": "ALDEx2",
        "ZicoSeq": "ZicoSeq", "limma_voom": "limma-voom"}
GLAB = {"any_ige": "any", "single": "single", "multiple": "multiple",
        "high_score": "high", "medium_score": "medium", "low_score": "low",
        **{f"cluster_{c}": f"cluster {c}" for c in "ABCDEFGHI"}}
RLAB = {"family": "family", "genus": "genus", "species": "species", "asv": "ASV"}
PANELS = [("prev10", "A", "10% prevalence filter"), ("prev05", "B", "5% prevalence filter"),
          ("prev01", "C", "1% prevalence filter")]

# ── counts per analysis x method ──────────────────────────────────────────
files = sorted((RES / "da_full").glob("da_long_*.csv")) + \
        sorted((RES / "da_extended_full").glob("daext_long_*.csv"))
cols = ["method", "q_within", "level", "prevalence", "partition", "n_ige", "n_taxa"]
d = pd.concat([pd.read_csv(f, usecols=cols) for f in files], ignore_index=True)
d = d[d.method.isin(METHODS) & d.prevalence.isin([p for p, _, _ in PANELS])]   # excludes prev00 (unfiltered, primary contrast only)
c = (d.assign(sig=d.q_within < ALPHA)
      .groupby(["level", "prevalence", "partition", "method"])
      .agg(n_sig=("sig", "sum"), n_taxa=("n_taxa", "first"), n_ige=("n_ige", "first"))
      .reset_index())
c["pct_sig"] = 100 * c.n_sig / c.n_taxa
assert len(c) == 4 * 3 * 15 * len(METHODS), f"incomplete grid: {len(c)} rows"
c.to_csv(RES / "tables" / "da_full_method_counts_by_analysis.csv", index=False)

# ── plot: one shared row order, rank blocks, one panel per prevalence ─────
cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "savefig.transparent": False,
                     "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "font.size": 9, "svg.fonttype": "none"})
# Discrete buckets of the significant-taxon count; 0 is near-white
EDGES = [0, 1, 2, 3, 6, 11, 51, 101, 10**6]          # [lo, hi) bins
BLABS = ["0", "1", "2", "3\u20135", "6\u201310", "11\u201350", "51\u2013100", ">100"]
BCOLS = ["#f4f4f4"] + [plt.get_cmap("YlOrRd")(v) for v in np.linspace(0.12, 0.95, 7)]
CMAP = ListedColormap(BCOLS)
norm = BoundaryNorm(EDGES, CMAP.N)
GAP_R, GAP_P, W_PANEL = 0.8, 0.9, 1 + 0.15 + len(METHODS)


def cell(ax, x, y, face, text, dark):
    ax.add_patch(plt.Rectangle((x, y), 1, 1, color=face, ec="white", lw=0.6))
    ax.text(x + 0.5, y + 0.5, text, ha="center", va="center", fontsize=7.5,
            color="white" if dark else "#222222")


def grey(v, lo, hi):
    s = (np.log10(v) - np.log10(lo)) / (np.log10(hi) - np.log10(lo))
    return plt.cm.Greys(0.12 + 0.75 * s), s > 0.55


fig, ax = plt.subplots(figsize=(15, 14))
n_lo, n_hi = c.n_ige.min(), c.n_ige.max()
f_lo, f_hi = c.n_taxa.min(), c.n_taxa.max()
y = 0.0
for level in ["family", "genus", "species", "asv"]:
    y_top = y
    for part, glab in GLAB.items():
        r = c[(c.level == level) & (c.partition == part)]
        ax.text(-0.15, y + 0.5, glab, ha="right", va="center", fontsize=8.5)
        face, dark = grey(r.n_ige.iloc[0], n_lo, n_hi)
        cell(ax, 0, y, face, str(int(r.n_ige.iloc[0])), dark)
        for k, (prev, _, _) in enumerate(PANELS):
            x0 = 1 + GAP_P + k * (W_PANEL + GAP_P)
            q = r[r.prevalence == prev].set_index("method").loc[METHODS]
            face, dark = grey(q.n_taxa.iloc[0], f_lo, f_hi)
            cell(ax, x0, y, face, str(int(q.n_taxa.iloc[0])), dark)
            for j, n in enumerate(q.n_sig.values):
                b = int(np.searchsorted(EDGES, n, side="right") - 1)
                cell(ax, x0 + 1.15 + j, y, BCOLS[b], str(int(n)), b >= 5)
        y += 1
    ax.text(-2.6, (y_top + y) / 2, RLAB[level], rotation=90, ha="center", va="center",
            fontsize=10, fontweight="bold")
    ax.plot([-2.3, -2.3], [y_top + 0.1, y - 0.1], color="#444444", lw=1.2)
    y += GAP_R

# column headers and panel titles
def header(x, s):
    ax.text(x + 0.5, -0.25, s, rotation=45, ha="left", va="bottom",
            rotation_mode="anchor", fontsize=9)

header(0, "n sensitized")
for k, (prev, tag, title) in enumerate(PANELS):
    x0 = 1 + GAP_P + k * (W_PANEL + GAP_P)
    header(x0, "features tested")
    for j, m in enumerate(METHODS):
        header(x0 + 1.15 + j, MLAB[m])
    ax.text(x0 + W_PANEL / 2, -6.4, title, ha="center", va="bottom", fontsize=10.5,
            fontweight="bold")
    ax.text(x0 - 0.2, -6.4, tag, ha="right", va="bottom", fontsize=15, fontweight="bold")

ax.set_xlim(-3.2, 1 + GAP_P + len(PANELS) * (W_PANEL + GAP_P))
ax.set_ylim(y, -7.6)
ax.axis("off")
cax = fig.add_axes([0.36, 0.075, 0.32, 0.013])
for i, (col, lab) in enumerate(zip(BCOLS, BLABS)):
    cax.add_patch(plt.Rectangle((i, 0), 1, 1, color=col, ec="white", lw=1))
    cax.text(i + 0.5, -0.4, lab, ha="center", va="top", fontsize=8.5)
cax.set_xlim(0, len(BCOLS)); cax.set_ylim(0, 1); cax.axis("off")
cax.text(len(BCOLS) / 2, 1.5, "significant taxa (q < 0.05)", ha="center", va="bottom", fontsize=9)
fig.text(0.52, 0.035,
         "Cell numbers: taxa with q < 0.05 (BH within analysis). Full cohort (275 sensitized, "
         "233 controls), adjusted for sex, age and BMI.\nGrey columns: darker = higher "
         "(log scale).",
         ha="center", fontsize=8.5, color="#444444")
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"figS_da_methods_nearing.{ext}", dpi=300, bbox_inches="tight",
                transparent=False, facecolor="white")
print("wrote", OUT / "figS_da_methods_nearing.{png,svg,pdf}")

tot = c.groupby(["prevalence", "method"])[["n_sig", "n_taxa"]].sum()
tot["pct"] = (100 * tot.n_sig / tot.n_taxa).round(2)
print(tot.pct.unstack("method")[METHODS])
