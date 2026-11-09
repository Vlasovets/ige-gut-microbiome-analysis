#!/usr/bin/env python
"""
Module 06k - corrected versions of the submitted Figure 4 volcano plots, for the
same four contrasts: any IgE, medium score, low score and cluster A, each vs
non-IgE. ASV level, LinDA, as in the submitted panels (notebook 03, cell 26).

Reviewer: "What's the reason for significance plateau around 0.6?" The submitted
panels plotted -log10(BH q); BH q-values are a step function (running minimum
of p*m/rank), so hundreds of ASVs share one q and stack on a horizontal line.

Corrections:
  y = -log10 UNADJUSTED p (no plateau);
  dashed line = BH q = 0.05 within the analysis, drawn in EVERY panel: at the
      largest raw p that passes, or, when nothing passes, at 0.05/m - the raw p
      the top-ranked taxon would need (BH rank 1);
  x = the method's log2 fold change;
  one threshold (0.05; the submitted panels drew 0.1);
  colour = family of the ASV, ONLY for q < 0.05 (four families of the submitted
      paper + Prevotellaceae in their original colours, other families navy);
      every non-significant ASV is grey (co-author review, 2026-10-22).

Design: full cohort, each group vs all 233 IgE-negative, ~ IgE status + sex +
age + BMI (02f --design full), 5% prevalence filter. Cluster A is 12 vs 233:
small-group hits were found to be group-size imbalance artifacts, so its panel
carries a caveat.

Options:
  --method  LinDA (default, as submitted) or any other revision method;
  --set     submitted (default: any, medium, low, cluster A) or
            score (any, high, medium, low IgE score).
Hits are labelled only when a panel has at most MAX_LABELS of them.

Runs in scont_plots (cnsplots), Roboto font.
Output: results/figures/review/fig4_volcanoes_contrasts[_score][_<method>].{png,svg,pdf}
        results/tables/fig4_volcanoes_contrasts[_score][_<method>].csv (source data)
        (no suffix for the submitted set / LinDA)
        --compact adds suffix _compact (layout for the Figure 4 composite, 06r)
Usage:  python revision/analysis/06k_fig4_volcanoes_contrasts.py [--method ANCOMBC2] [--set score] [--compact]
"""

import argparse
from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
OUT = RES / "figures" / "review"
KORA = Path("/home/itg/oleg.vlasovets/slr_example/KORA_data")
LEVEL, PREV, ALPHA = "asv", "prev05", 0.05
SETS = {"submitted": [("any_ige", "Any IgE"), ("medium_score", "Medium score"),
                      ("low_score", "Low score"), ("cluster_A", "Cluster A")],
        "score": [("any_ige", "Any IgE"), ("high_score", "High score"),
                  ("medium_score", "Medium score"), ("low_score", "Low score")]}
FAM = {"Lachnospiraceae": "#fcc853", "Bacteroidaceae": "#078a07",
       "Oscillospiraceae": "#1f77b4", "Veillonellaceae": "#ff7f0e",
       "Prevotellaceae": "#c2185b"}
OTHER, NAVY = "#c9ccd3", "#2a3f5f"
MLAB = {"LinDA": "LinDA", "ANCOMBC2": "ANCOM-BC2", "ALDEx2": "ALDEx2",
        "ZicoSeq": "ZicoSeq", "limma_voom": "limma-voom"}
ap = argparse.ArgumentParser()
ap.add_argument("--method", default="LinDA", choices=list(MLAB))
ap.add_argument("--set", default="submitted", choices=list(SETS))
# wider-than-tall layout for the Figure 4 composite (06r): two methods stacked
# beside panels A+B; legend only when a panel has hits (the ANCOM-BC2 legend
# then serves both methods)
ap.add_argument("--compact", action="store_true")
args = ap.parse_args()
METHOD, CONTRASTS = args.method, SETS[args.set]
STEM = ("fig4_volcanoes_contrasts" + ("" if args.set == "submitted" else f"_{args.set}")
        + ("" if METHOD == "LinDA" else f"_{METHOD.lower()}")
        + ("_compact" if args.compact else ""))
MAX_LABELS = 10

# ── data ──────────────────────────────────────────────────────────────────
tax = pd.read_csv(KORA / "taxonomy_clean.csv", index_col=0)
clean = lambda s, p: str(s).replace(p, "").strip().rstrip(";").strip()
tax["fam"] = tax.family.map(lambda s: clean(s, "f__"))
tax["gen"] = tax.genus.map(lambda s: clean(s, "g__"))

frames = []
for part, _ in CONTRASTS:
    d = pd.concat([pd.read_csv(RES / "da_full" / f"da_long_{LEVEL}_{PREV}_{part}.csv"),
                   pd.read_csv(RES / "da_extended_full" / f"daext_long_{LEVEL}_{PREV}_{part}.csv")])
    d = d[d.method == METHOD].copy()
    assert len(d) and d.taxon.is_unique and (d.n_ctrl == 233).all(), part
    frames.append(d)
d = pd.concat(frames, ignore_index=True)
d["family"] = d.taxon.map(tax.fam)
d["genus"] = d.taxon.map(tax.gen)
assert d.family.notna().all(), "ASV missing from taxonomy"
d["y"] = -np.log10(d.p_raw)
d["sig"] = d.q_within < ALPHA
d[["partition", "taxon", "family", "genus", "effect_log2", "p_raw", "q_within",
   "n_ige", "n_ctrl", "formula", "method", "effect_type", "passed_ss"]] \
  .to_csv(RES / "tables" / f"{STEM}.csv", index=False)

# ── style ─────────────────────────────────────────────────────────────────
cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "savefig.transparent": False,
                     "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "font.size": 9, "axes.labelsize": 9.5,
                     "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "svg.fonttype": "none"})

# ── figure ────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(6.8, 4.8) if args.compact else (8.4, 7.4), sharey=True,
                         gridspec_kw={"wspace": 0.15, "hspace": 0.75 if args.compact else 0.5})
ymax = max(d.y.max(), -np.log10(ALPHA / d.groupby('partition').size().max())) * 1.1
ringed = False
for k, (ax, (part, title)) in enumerate(zip(axes.flat, CONTRASTS)):
    t = d[d.partition == part]
    ns = t[~t.sig]
    ax.scatter(ns.effect_log2, ns.y, s=10, c=OTHER, lw=0, zorder=1)
    s = t[t.sig]
    so = s[~s.family.isin(FAM)]
    ax.scatter(so.effect_log2, so.y, s=14, c=NAVY, lw=0, zorder=2)
    for fname, col in FAM.items():
        h = s[s.family == fname]
        ax.scatter(h.effect_log2, h.y, s=18, c=col, ec="white", lw=0.4, zorder=3)
    if len(s) and len(s) <= MAX_LABELS:   # with many hits the BH line alone marks them
        ax.scatter(s.effect_log2, s.y, s=34, facecolors="none", ec="black", lw=1.0, zorder=4)
        ringed = True
    bh_p = s.p_raw.max() if len(s) else ALPHA / len(t)   # rank-1 BH requirement
    ax.axhline(-np.log10(bh_p), ls="--", lw=0.8, color="#d62728", zorder=0)
    # labels left of the points, staggered so close hits do not overlap
    label_y = np.inf
    for _, r in (s if len(s) <= MAX_LABELS else s.iloc[:0]).sort_values("y", ascending=False).iterrows():
        label_y = min(r.y, label_y - 0.42)
        ax.annotate(r.genus.replace("_", " "), (r.effect_log2, r.y),
                    xytext=(r.effect_log2 - 0.12, label_y), textcoords="data",
                    ha="right", va="center", fontsize=7, style="italic", zorder=5,
                    bbox={"boxstyle": "square,pad=0.1", "fc": "white", "ec": "none"},
                    arrowprops={"arrowstyle": "-", "lw": 0.5, "color": "#666666",
                                "shrinkA": 1, "shrinkB": 3})
    ax.axvline(0, lw=0.6, color="#999999", zorder=0)
    lo, hi = ax.get_xlim(); r = max(abs(lo), abs(hi)); ax.set_xlim(-r, r)
    ax.set_ylim(0, ymax)
    ax.set_title(f"{title} vs non-IgE\n(n = {t.n_ige.iat[0]} vs {t.n_ctrl.iat[0]};"
                 f" {int(t.sig.sum())} of {len(t)} ASVs q < 0.05)",
                 fontsize=8.5 if args.compact else 9.5, color=NAVY)
    ax.set_xlabel(f"log$_2$ fold change ({MLAB[METHOD]})")
    if k % 2 == 0:
        ax.set_ylabel("$-$log$_{10}$ p (unadjusted)")
    ax.spines[["top", "right"]].set_visible(False)
    if part == "cluster_A":
        ax.text(0.03, 0.97, "small group: interpret with caution", transform=ax.transAxes,
                fontsize=7.5, color="#666666", va="top", style="italic")

h = []
if d.sig.any():
    h += [Line2D([], [], marker="o", ls="", mfc=c, mec="white", ms=7, label=f"{f}, q < 0.05")
          for f, c in FAM.items() if (d.sig & (d.family == f)).any()]
    h += [Line2D([], [], marker="o", ls="", mfc=NAVY, mec=NAVY, ms=5, label="other family, q < 0.05")]
h += [Line2D([], [], marker="o", ls="", mfc=OTHER, mec=OTHER, ms=5, label="not significant"),
      Line2D([], [], ls="--", color="#d62728", label="BH q = 0.05")]
if ringed:
    h += [Line2D([], [], marker="o", ls="", mfc="none", mec="black", ms=7, label="labelled hit")]
if not args.compact or d.sig.any():
    fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.02 if args.compact else 0.03),
               ncol=4, frameon=False, fontsize=8 if args.compact else 8.5, handletextpad=0.3,
               columnspacing=1.2)

OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"{STEM}.{ext}", dpi=300, bbox_inches="tight",
                transparent=False, facecolor="white")
print("wrote", OUT / f"{STEM}.{{png,svg,pdf}}")
print(d.groupby("partition").agg(n_ige=("n_ige", "first"), tested=("taxon", "size"),
                                 n_sig=("sig", "sum"), min_p=("p_raw", "min"),
                                 min_q=("q_within", "min")))
print(d[d.sig][["partition", "family", "genus", "effect_log2", "p_raw", "q_within"]]
      .sort_values(["partition", "q_within"]).to_string(index=False))
