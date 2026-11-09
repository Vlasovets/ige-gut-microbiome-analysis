#!/usr/bin/env python
"""
Module 06s - Figure 2B redrawn: the 275 sensitized individuals in the space of
the three latent allergy components (LACs), coloured by IgE cluster.

Reviewer 2: cluster colours were hard to tell apart in the 3D panel. Each
cluster is now also labelled with its letter at its centroid, connected by a
leader line to a label placed outside the point cloud, so the colour is not the
only cue. Colours = published cluster colours (Figure 1; hardcoded in 06g).

LAC scores and clusters are recomputed with fig1a_order.recompute() (notebook
01_IgE_analysis.ipynb, verified against KORA_data/cluster_df_latent.csv); no
participant-level output is written.

Runs in scont_plots (cnsplots), Roboto font.
Output: results/figures/review/fig2B_lac_clusters.{png,svg,pdf}
Usage:  python revision/analysis/06s_fig2B_lac_clusters.py
"""

from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from mpl_toolkits.mplot3d import proj3d

from fig1a_order import recompute

OUT = Path(__file__).resolve().parent.parent / "results" / "figures" / "review"
COL = {"A": "#0072B2", "B": "#E69F00", "C": "#CC79A7", "D": "#D55E00", "E": "#F0E442",
       "F": "#009E73", "G": "#999999", "H": "#56B4E9", "I": "#A6761D"}
NAVY = "#2a3f5f"

r = recompute()
lac, cluster = r["lac"], r["cluster"].reindex(r["lac"].index)
assert cluster.notna().all() and len(lac) == 275

cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "savefig.transparent": False,
                     "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "font.size": 12, "svg.fonttype": "none"})

fig = plt.figure(figsize=(8.2, 7.2))
ax = fig.add_subplot(projection="3d")
ax.view_init(elev=18, azim=45)
ax.set_box_aspect(None, zoom=0.88)
for c in sorted(COL):
    v = lac[cluster == c]
    ax.scatter(v.LAC1, v.LAC2, v.LAC3, s=34, c=COL[c], edgecolors="#333333", linewidths=0.4,
               depthshade=False, label=f"{c} (n = {len(v)})")

# Leader lines: centroid (open circle) -> letter label. Offsets (points) are set per
# cluster in screen space so labels sit outside the point cloud for this view angle.
OFFSET = {"A": (-35, 30), "B": (-25, -25), "C": (-25, -35), "D": (55, -55), "E": (45, -15),
          "F": (45, 20), "G": (0, 45), "H": (30, 75), "I": (-80, 25)}
cent = lac.groupby(cluster).mean()
for c, xyz in cent.iterrows():
    x2, y2, _ = proj3d.proj_transform(*xyz.values, ax.get_proj())
    ax.scatter(*xyz.values, s=60, c="white", edgecolors=NAVY, linewidths=1.4, zorder=11, depthshade=False)
    ax.annotate(c, xy=(x2, y2), xytext=OFFSET[c], textcoords="offset points", fontsize=15,
                fontweight="bold", ha="center", va="center", zorder=12,
                bbox={"boxstyle": "round,pad=0.25", "fc": COL[c], "ec": NAVY, "lw": 0.8},
                arrowprops={"arrowstyle": "-", "color": NAVY, "lw": 1.0, "shrinkA": 0, "shrinkB": 4})

ax.set_xlabel("LAC$_1$", labelpad=8); ax.set_ylabel("LAC$_2$", labelpad=8); ax.set_zlabel("LAC$_3$", labelpad=6)
for a in (ax.xaxis, ax.yaxis, ax.zaxis):
    a.pane.set_facecolor("white"); a.pane.set_edgecolor("#dddddd")
ax.tick_params(labelsize=10)
ax.legend(title="Cluster", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False,
          fontsize=11, title_fontsize=12, markerscale=1.3)

OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"fig2B_lac_clusters.{ext}", dpi=300, bbox_inches="tight",
                transparent=False, facecolor="white")
print("wrote", OUT / "fig2B_lac_clusters.{png,svg,pdf}")
print(cent.round(2).to_string())
