#!/usr/bin/env python
"""
Module 06r - revised Figure 5 (log-contrast classification), same layout as the
submitted figure, rebuilt from the module 04 --classo results.

  A  family-level stability: for each IgE group (columns) the fraction of the 50
     training folds (5 folds x 10 repeats, matched pairs) in which c-lasso with
     stability selection retained the family (rows). Families retained in at
     least 5% of folds in any group are shown, plus the five families coloured
     in Figure 4 (same colours), which are always shown.
  B  out-of-sample AUC of the 46 analyses (15 IgE groups x family, genus,
     species; ASV for all sensitized) for each lambda selection scheme;
     dashed line = 0.5 (chance), red line = median.

Submitted version for comparison: per-group selection probabilities from one
stability-selection fit and MCR histograms over 315 models (7 ranks, relative
abundances, sex-only matching).

Runs in scont_plots (cnsplots), Roboto font.
Output: results/figures/review/fig5_classification.{png,svg,pdf}
        results/tables/fig5A_family_selection.csv (source data)
Usage:  python revision/analysis/06r_fig5_classification.py
"""

from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
OUT = RES / "figures" / "review"
NAVY = "#2a3f5f"
GROUPS = [("any_ige", "Any"), ("single", "Single"), ("multiple", "Multiple"),
          ("high_score", "High"), ("medium_score", "Medium"), ("low_score", "Low")] + \
         [(f"cluster_{c}", c) for c in "ABCDEFGHI"]
FAM = {"Lachnospiraceae": "#fcc853", "Bacteroidaceae": "#078a07",
       "Oscillospiraceae": "#1f77b4", "Veillonellaceae": "#ff7f0e",
       "Prevotellaceae": "#c2185b"}
SCHEMES = [("classo-theoretical", r"theoretical $\lambda$"),
           ("classo-CV", "cross-validation"), ("classo-StabSel", "stability selection")]
MIN_FREQ = 0.05

# ── data ──────────────────────────────────────────────────────────────────
cols = {}
for part, lab in GROUPS:
    f = RES / "ml" / f"ml_family_prev10_{part}_classo_selection.csv"
    cols[lab] = pd.read_csv(f, index_col=0)["classo-StabSel"]
A = pd.DataFrame(cols).fillna(0.0)
# highlighted families are always shown (also when never selected), so the
# figure can be read against Figure 4
A = A.reindex(A.index.union(list(FAM))).fillna(0.0)
A = A[(A >= MIN_FREQ).any(axis=1) | A.index.isin(list(FAM))].sort_index(ascending=False)
A.to_csv(RES / "tables" / "fig5A_family_selection.csv")

t = pd.read_csv(RES / "tables" / "ml_benchmark_all_analyses.csv")
t = t[t.model.str.startswith("classo")]
assert t.groupby("model").size().eq(46).all(), "c-lasso grid incomplete"

# ── style ─────────────────────────────────────────────────────────────────
cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "savefig.transparent": False,
                     "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "font.size": 9, "svg.fonttype": "none"})

fig = plt.figure(figsize=(12, 0.32 * len(A) + 2.6))
gs = fig.add_gridspec(3, 2, width_ratios=[2.3, 1], wspace=0.25, hspace=0.55)

# A: heatmap
ax = fig.add_subplot(gs[:, 0])
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list("w_red", ["#ffffff", "#ffcccc", "#ff6666", "#990000"])   # as the submitted panel
im = ax.imshow(A.values, cmap=cmap, vmin=0, vmax=1, aspect="auto")
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        v = A.values[i, j]
        ax.text(j, i, f"{v:.2f}".rstrip("0").rstrip(".") if v else "0", ha="center",
                va="center", fontsize=7, color="white" if v > 0.6 else "#333333")
ax.set_xticks(range(A.shape[1])); ax.set_xticklabels(A.columns, rotation=45, ha="left",
                                                       rotation_mode="anchor")
ax.xaxis.tick_top(); ax.tick_params(length=0)
ax.set_yticks(range(A.shape[0])); ax.set_yticklabels([n.replace("_", " ") for n in A.index])
for tl in ax.get_yticklabels():
    if tl.get_text() in FAM:
        tl.set_bbox({"facecolor": FAM[tl.get_text()], "edgecolor": "none", "pad": 1.5})
        tl.set_color("white" if tl.get_text() != "Lachnospiraceae" else "#222222")
for s in ax.spines.values():
    s.set_visible(False)
cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.03, pad=0.02, aspect=30)
cb.set_label("fraction of training folds selecting the family (stability selection)", fontsize=8.5)
ax.text(-0.42, 1.04, "A", transform=ax.transAxes, fontsize=16, fontweight="bold")

# B: AUC histograms per scheme
bins = np.arange(0.35, 0.775, 0.025)
for k, (m, lab) in enumerate(SCHEMES):
    a = fig.add_subplot(gs[k, 1])
    v = t.loc[t.model == m, "auc_pooled"].values
    a.hist(v, bins=bins, color=NAVY, ec="white", lw=0.5)
    a.axvline(0.5, color="#888888", ls="--", lw=0.9)
    a.axvline(np.median(v), color="#d62728", lw=1.2)
    a.set_title(f"Log-contrast classification ({lab})", fontsize=9, color=NAVY)
    a.set_ylabel("analyses")
    a.set_xlim(bins[0], bins[-1])
    a.spines[["top", "right"]].set_visible(False)
    a.text(0.98, 0.92, f"median {np.median(v):.2f}", transform=a.transAxes, ha="right",
           va="top", fontsize=8, color="#d62728")
    if k == 0:
        a.text(-0.25, 1.15, "B", transform=a.transAxes, fontsize=16, fontweight="bold")
a.set_xlabel("out-of-sample AUC")

OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"fig5_classification.{ext}", dpi=300, bbox_inches="tight",
                transparent=False, facecolor="white")
print("wrote", OUT / "fig5_classification.{png,svg,pdf}", "| families shown:", len(A))
print(t.groupby("model").auc_pooled.median().round(3).to_string())
