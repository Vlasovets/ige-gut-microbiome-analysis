#!/usr/bin/env python
"""
Module 06t - supplementary figure of the re-computed taxon set enrichment
(TaxSEA), replacing the submitted TaxSEA figures (fig:taxsea, enrichment_big_figure).

Primary analysis (module 07 --design full): species ranked by the signed LinDA
statistic of all 275 sensitized vs 233 non-sensitized individuals adjusted for
sex, age and BMI, 5% prevalence filter, TaxSEA 1.2.0 bundled database.

  A  all metabolite-producer sets: -log10 p of the two-sided KS test, coloured by
     the median statistic of the set members; FDR within the sub-database.
  B  cumulative distribution of the LinDA statistic for the folic acid and the
     vitamin A producer sets against the TaxSEA background (what the KS test compares).
  C  species-level DA (LinDA log2 fold change vs -log10 p); set members marked.

Inputs : results/tables/taxsea_recomputed_full.csv           (module 07)
         results/tables/taxsea_species_membership_full.csv   (module 07b)
Output : results/figures/review/figS_taxsea_null.{png,svg,pdf}
Env    : scont_plots (cnsplots), Roboto
Usage  : python revision/analysis/06t_taxsea_null.py
"""

import re
from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

REV = Path(__file__).resolve().parent.parent
TAB = REV / "results" / "tables"
OUT = REV / "results" / "figures" / "review"

UP, DOWN, GREY = "#C44E52", "#4C72B0", "#BBBBBB"
FOLIC, VITA = "#0072B2", "#E69F00"
HEADLINE = {"MiMeDB_producers_of_Folic_acid": "Folic acid", "MiMeDB_producers_of_Vitamin_A": "Vitamin A"}

t = pd.read_csv(TAB / "taxsea_recomputed_full.csv")
sets = t[(t.sub_database == "Metabolite_producers") & (t.prevalence == "prev05") & (t.rank_by == "stat")].copy()
sp = pd.read_csv(TAB / "taxsea_species_membership_full.csv")
bg = sp[sp.in_background]
assert len(sets) == 18 and len(bg) == 53, (len(sets), len(bg))


def set_label(s):
    src = "MiMeDB" if s.startswith("MiMeDB") else "gutMGene"
    v = re.sub(r"^(MiMeDB|GutMGene)_producers_of_", "", s).strip("_").replace("_", " ")
    v = {"3IAA": "Indole-3-acetic acid", "TMAO": "TMAO", "GABA": "GABA"}.get(v, v)
    return f"{v[0].upper() + v[1:]} ({src})"


def ecdf(x):
    x = np.sort(np.asarray(x))
    return np.r_[x[0], x], np.r_[0, np.arange(1, len(x) + 1) / len(x)]


cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "font.size": 9, "axes.titlesize": 10,
                     "savefig.transparent": False, "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "svg.fonttype": "none"})

fig = plt.figure(figsize=(11, 4.1))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1], wspace=0.42, left=0.2, right=0.98, top=0.88, bottom=0.14)

# A - all metabolite-producer sets
ax = fig.add_subplot(gs[0])
s = sets.sort_values("PValue", ascending=False).reset_index(drop=True)
y = np.arange(len(s))
col = np.where(s.median_rank_of_set_members > 0, UP, DOWN)
ax.barh(y, -np.log10(s.PValue), color=col, height=0.7)
ax.set_yticks(y)
labels = ax.set_yticklabels([set_label(v) for v in s.taxonSetName], fontsize=8)
for lab, name in zip(labels, s.taxonSetName):
    if name in HEADLINE:
        lab.set_fontweight("bold")
ax.axvline(-np.log10(0.05), color="#555555", ls="--", lw=0.8)
ax.text(-np.log10(0.05) + 0.02, len(s) - 0.6, "p = 0.05", fontsize=7.5, color="#555555", va="top")
ax.set_xlabel(r"$-\log_{10}$ p (two-sided KS test)")
ax.set_xlim(0, 1.45)
ax.set_title(f"Metabolite producer sets (n = {len(s)})\nsmallest FDR = {s.FDR.min():.2f}", loc="left")
ax.legend(handles=[Patch(color=UP, label="higher in sensitized"), Patch(color=DOWN, label="lower in sensitized")],
          title="median of set members", fontsize=7.5, title_fontsize=7.5, frameon=False, loc="lower right")
ax.text(0.005, 0.975, "A", transform=fig.transFigure, fontsize=14, fontweight="bold", ha="left", va="top")

# B - ECDF: set members vs background
ax = fig.add_subplot(gs[1])
x, yv = ecdf(bg.stat)
ax.step(x, yv, where="post", color=GREY, lw=2.2, label=f"all species (n = {len(bg)})")
for key, name, c in (("in_folic_acid", "MiMeDB_producers_of_Folic_acid", FOLIC),
                     ("in_vitamin_a", "MiMeDB_producers_of_Vitamin_A", VITA)):
    m = bg[bg[key]]
    p = sets.loc[sets.taxonSetName == name, "PValue"].iloc[0]
    x, yv = ecdf(m.stat)
    ax.step(x, yv, where="post", color=c, lw=1.6,
            label=f"{HEADLINE[name]} producers\n(n = {len(m)}), KS p = {p:.2f}")
ax.axvline(0, color="#999999", lw=0.6)
ax.set_xlabel("LinDA statistic (IgE vs. non-IgE)")
ax.set_ylabel("Cumulative fraction of species")
ax.set_title("Set members vs. background", loc="left")
ax.legend(fontsize=7.5, frameon=False, loc="upper left", handlelength=1.2)
ax.text(-0.2, 1.1, "B", transform=ax.transAxes, fontsize=14, fontweight="bold", ha="left", va="top")

# C - species volcano with set members
ax = fig.add_subplot(gs[2])
ly = -np.log10(sp.p_raw)
ax.scatter(sp.effect_log2, ly, s=16, c=GREY, edgecolors="none", label="other species")
for key, c, lab in (("in_folic_acid", FOLIC, "folic acid producer"), ("in_vitamin_a", VITA, "vitamin A producer")):
    m = sp[sp[key]]
    ax.scatter(m.effect_log2, -np.log10(m.p_raw), s=24 if key == "in_folic_acid" else 55,
               facecolors=c if key == "in_folic_acid" else "none", edgecolors=c,
               linewidths=0 if key == "in_folic_acid" else 1.3, label=lab, zorder=3)
ax.axhline(-np.log10(0.05), color="#555555", ls="--", lw=0.8)
for _, r in sp[sp.p_raw < 0.05].iterrows():
    ax.annotate(r.taxon.replace("_", " ").replace("Prevotella ", "P. ").replace("Bacteroides ", "B. "),
                (r.effect_log2, -np.log10(r.p_raw)), xytext=(6, -2), textcoords="offset points",
                fontsize=7.5, fontstyle="italic")
ax.set_xlabel(r"LinDA $\log_2$ fold change")
ax.set_ylabel(r"$-\log_{10}$ p")
ax.set_title(f"Species (n = {len(sp)}); smallest q = {sp.q_within.min():.2f}", loc="left")
ax.legend(fontsize=7.5, frameon=False, loc="upper right")
ax.text(-0.2, 1.1, "C", transform=ax.transAxes, fontsize=14, fontweight="bold", ha="left", va="top")

OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"figS_taxsea_null.{ext}", dpi=300, transparent=False, facecolor="white")
print("wrote", OUT / "figS_taxsea_null.{png,svg,pdf}")
print(f"folic acid members {int(bg.in_folic_acid.sum())}/{len(bg)} background species; "
      f"vitamin A {int(bg.in_vitamin_a.sum())}; min set FDR {sets.FDR.min():.2f}; min species q {sp.q_within.min():.3f}")
