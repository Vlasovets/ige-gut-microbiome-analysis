#!/usr/bin/env python
"""
Module 06i - revised Figure 4B-C (panel A, alpha diversity, is kept as in the
manuscript). Replaces the 61-ASV heatmap, the old volcanoes and the pie charts.

  B  volcano per taxonomic rank (family, genus, species, ASV) for the primary
     contrast: all sensitized (n=275) vs all IgE-negative (n=233), full cohort,
     ~ IgE status + sex + age + BMI, 5% prevalence filter.
       x = LinDA log2 fold change, y = -log10 RAW p (no BH plateau);
       dashed line = largest raw p with LinDA BH q < 0.05 (drawn only if any);
       colour = family of the taxon (the four families of the submitted paper
       + Prevotellaceae; grey = other); black ring = consensus (q < 0.05 in
       >= 2 of 5 methods, concordant sign). Family names are labelled in the
       family panel; consensus taxa are labelled in every panel.
  C  method agreement for the highlighted families: signed -log10 p
     (sign = direction in sensitized) per method, full cohort and the matched
     sensitivity design side by side; cell text = BH q within that analysis
     (bold = q < 0.05). No taxon reaches LinDA q < 0.05 at any rank, so the
     volcanoes carry no BH line.

Runs in scont_plots (cnsplots), Roboto font.
Output: results/figures/review/fig4BC_da.{png,svg,pdf}
        results/tables/fig4C_family_agreement.csv (source data)
Usage:  python revision/analysis/06i_fig4_da.py
"""

from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
OUT = RES / "figures" / "review"
KORA = Path("/home/itg/oleg.vlasovets/slr_example/KORA_data")
PREV, PART, ALPHA = "prev05", "any_ige", 0.05
METHODS = ["LinDA", "ANCOMBC2", "ALDEx2", "ZicoSeq", "limma_voom"]
MLAB = {"LinDA": "LinDA", "ANCOMBC2": "ANCOM-BC2", "ALDEx2": "ALDEx2",
        "ZicoSeq": "ZicoSeq", "limma_voom": "limma-voom"}
# Colours of the submitted figure, so readers can map old to new.
FAM = {"Lachnospiraceae": "#fcc853", "Bacteroidaceae": "#078a07",
       "Oscillospiraceae": "#1f77b4", "Veillonellaceae": "#ff7f0e",
       "Prevotellaceae": "#c2185b"}
OTHER, NAVY = "#c9ccd3", "#2a3f5f"
LEVELS = [("family", "Family"), ("genus", "Genus"), ("species", "Species"), ("asv", "ASV")]

# ── data ──────────────────────────────────────────────────────────────────
def long(design, level):
    da, ext = ("da_full", "da_extended_full") if design == "full" else ("da", "da_extended")
    f1 = RES / da / f"da_long_{level}_{PREV}_{PART}.csv"
    f2 = RES / ext / f"daext_long_{level}_{PREV}_{PART}.csv"
    d = pd.concat([pd.read_csv(f1), pd.read_csv(f2)], ignore_index=True)
    return d[d.method.isin(METHODS)]


def consensus(d):
    g = d.assign(sig=d.q_within < ALPHA).groupby("taxon")
    n_sig = g.sig.sum()
    conc = g.effect_log2.apply(lambda e: np.sign(e.dropna()).nunique() == 1)
    return (n_sig >= 2) & conc


tax = pd.read_csv(KORA / "taxonomy_clean.csv", index_col=0)
clean = lambda s, p: str(s).replace(p, "").strip().rstrip(";").strip()
tax["fam"] = tax.family.map(lambda s: clean(s, "f__"))
tax["gen"] = tax.genus.map(lambda s: clean(s, "g__"))
tax["spe"] = tax.species.map(lambda s: clean(s, "s__"))
fam_of = {"family": lambda t: t,
          "genus": tax.groupby("gen").fam.agg(lambda x: x.mode().iat[0]).to_dict().get,
          "species": tax.groupby("spe").fam.agg(lambda x: x.mode().iat[0]).to_dict().get,
          "asv": tax.fam.to_dict().get}

# ── figure ────────────────────────────────────────────────────────────────
cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
# Roboto (latin subset) lacks some symbols (e.g. ≥); DejaVu Sans fills those glyphs.
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "savefig.transparent": False,
                     "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "font.size": 9, "axes.labelsize": 9.5,
                     "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "svg.fonttype": "none"})
fig = plt.figure(figsize=(14, 8.4))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 0.78], hspace=0.62, wspace=0.28)

for k, (level, title) in enumerate(LEVELS):
    ax = fig.add_subplot(gs[0, k])
    d = long("full", level)
    cons = consensus(d)
    lin = d[d.method == "LinDA"].set_index("taxon")
    lin["fam"] = [fam_of[level](t) for t in lin.index]
    lin["y"] = -np.log10(lin.p_raw)
    lin["cons"] = cons.reindex(lin.index).fillna(False)
    other = lin[~lin.fam.isin(FAM)]
    ax.scatter(other.effect_log2, other.y, s=9, c=OTHER, lw=0, zorder=1)
    for fname, col in FAM.items():
        h = lin[lin.fam == fname]
        ax.scatter(h.effect_log2, h.y, s=36 if level == "family" else 16, c=col,
                   ec="white", lw=0.4, zorder=3)
    c = lin[lin.cons]
    ax.scatter(c.effect_log2, c.y, s=70, facecolors="none", ec="black", lw=1.0, zorder=4)
    sig = lin[lin.q_within < ALPHA]
    if len(sig):
        ax.axhline(-np.log10(sig.p_raw.max()), ls="--", lw=0.8, color="#d62728")
    ax.axvline(0, lw=0.6, color="#999999", zorder=0)
    labels = list(lin.index[lin.fam.isin(FAM)]) if level == "family" else []
    labels += [t for t in c.index if t not in labels and level != "asv"]
    for t in labels:
        r = lin.loc[t]
        ax.annotate(t.replace("_", " "), (r.effect_log2, r.y), xytext=(4, 3),
                    textcoords="offset points", fontsize=7.5,
                    style="italic" if level == "species" else "normal",
                    color=FAM.get(r.fam, "black") if level == "family" else "black")
    ax.set_title(f"{title} ({len(lin)} tested)", fontsize=10, color=NAVY)
    ax.set_xlabel("log$_2$ fold change (LinDA)")
    if k == 0:
        ax.set_ylabel("$-$log$_{10}$ p")
    ax.set_ylim(0, max(lin.y.max() * 1.12, 2.0))
    lo, hi = ax.get_xlim(); ax.set_xlim(lo, hi + 0.25 * (hi - lo))   # room for labels
    ax.spines[["top", "right"]].set_visible(False)
    if k == 0:
        ax.text(-0.28, 1.12, "B", transform=ax.transAxes, fontsize=16, fontweight="bold")

h = [Line2D([], [], marker="o", ls="", mfc=c, mec="white", ms=7, label=f) for f, c in FAM.items()]
h += [Line2D([], [], marker="o", ls="", mfc=OTHER, mec="white", ms=6, label="other families"),
      Line2D([], [], marker="o", ls="", mfc="none", mec="black", ms=8,
             label="consensus (≥2 of 5 methods)")]
fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, 0.505), ncol=8,
           frameon=False, fontsize=8.5, handletextpad=0.3, columnspacing=1.2)

# ── C: agreement matrix ───────────────────────────────────────────────────
rows = []
for design in ["full", "matched"]:
    d = long(design, "family")
    d = d[d.taxon.isin(FAM)]
    for _, r in d.iterrows():
        rows.append({"design": design, "family": r.taxon, "method": r.method,
                     "effect_log2": r.effect_log2, "p_raw": r.p_raw, "q_within": r.q_within,
                     "signed_logp": np.sign(r.effect_log2) * -np.log10(r.p_raw)})
m = pd.DataFrame(rows)
m.to_csv(RES / "tables" / "fig4C_family_agreement.csv", index=False)

ax = fig.add_subplot(gs[1, 1:3])
fams = list(FAM)
norm = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
cmap = plt.get_cmap("RdBu_r")
for gi, (design, glab) in enumerate([("full", "Full cohort (275 vs 233)"),
                                     ("matched", "Matched pairs (221 vs 221)")]):
    x0 = gi * (len(METHODS) + 0.6)
    for i, f in enumerate(fams):
        for j, meth in enumerate(METHODS):
            r = m[(m.design == design) & (m.family == f) & (m.method == meth)]
            if r.empty:
                ax.add_patch(plt.Rectangle((x0 + j, i), 1, 1, color="#f2f2f2", ec="white"))
                continue
            r = r.iloc[0]
            v = float(np.clip(r.signed_logp, -3, 3))
            ax.add_patch(plt.Rectangle((x0 + j, i), 1, 1, color=cmap(norm(v)), ec="white", lw=1))
            q = r.q_within
            txt = "<0.01" if q < 0.01 else f"{q:.2f}"
            ax.text(x0 + j + 0.5, i + 0.5, txt, ha="center", va="center", fontsize=8,
                    fontweight="bold" if q < ALPHA else "normal",
                    color="white" if abs(v) > 2 else "#333333")
    ax.text(x0 + len(METHODS) / 2, -2.9, glab, ha="center", fontsize=9.5, color=NAVY)
    for j, meth in enumerate(METHODS):
        ax.text(x0 + j + 0.5, -0.15, MLAB[meth], rotation=40, ha="left", va="bottom",
                rotation_mode="anchor", fontsize=8.5)
for i, f in enumerate(fams):
    ax.text(-0.15, i + 0.5, f, ha="right", va="center", fontsize=9, color=FAM[f],
            fontweight="bold")
ax.set_xlim(-0.05, 2 * len(METHODS) + 0.65); ax.set_ylim(len(fams), -3.4)
ax.axis("off")
ax.text(-0.30, 1.25, "C", transform=ax.transAxes, fontsize=16, fontweight="bold")
cax = ax.inset_axes([1.04, 0.05, 0.025, 0.9])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
cb.set_label("signed $-$log$_{10}$ p\n(red = higher in sensitized)", fontsize=8.5)
cb.ax.tick_params(labelsize=8)

OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"fig4BC_da.{ext}", dpi=300, bbox_inches="tight",
                transparent=False, facecolor="white")
print("wrote", OUT / "fig4BC_da.{png,svg,pdf}")
print(m.pivot_table(index=["design", "family"], columns="method", values="q_within").round(3))
