#!/usr/bin/env python
"""
Module 06e - revised Figure 3B-D (Figure 3A is kept as the author's panel).

  B  method heatmap in the style of Nearing et al. 2022: % of tested taxa
     called significant (q < 0.05), per IgE group x DA method, pooled over all
     ranks and prevalence filters; top bar = overall rate per method
  C  corrected volcano, primary contrast (all sensitized vs matched controls,
     species rank, 5% prevalence): y = -log10(RAW p) so the BH step-function
     plateau disappears; BH q < 0.05 drawn as the equivalent p cut-off; points
     coloured by how many of the 6 methods call the taxon
  D  stability of the one consensus taxon (P. copri aggregate) across the 6
     methods and 3 prevalence filters: effect +- 95% CI (method-specific
     scales); filled = q < 0.05 within that analysis

All values are read from result files. Answers R1 (more DA methods) and R2
(Fig 3C volcano plateau; Fig 3D pie charts). Runs in scont_plots (cnsplots).
Output: results/figures/review/fig3{B,C,D}_*.{png,svg,pdf} and fig3BCD.*
"""

from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

REV = Path(__file__).resolve().parent.parent
T = REV / "results" / "tables"
DAEXT = REV / "results" / "da_extended"
OUT = REV / "results" / "figures" / "review"
OUT.mkdir(parents=True, exist_ok=True)

cns.setup_matplotlib()
# cnsplots sets savefig.transparent=True (renders black in dark-mode viewers)
plt.rcParams.update({"savefig.transparent": False, "savefig.facecolor": "white",
                     "figure.facecolor": "white", "axes.facecolor": "white",
                     "font.size": 10, "axes.labelsize": 10, "xtick.labelsize": 9.5,
                     "ytick.labelsize": 9.5, "legend.fontsize": 9, "svg.fonttype": "none"})
NAVY = "#2a3f5f"
METHODS = ["LinDA", "ANCOMBC2", "ALDEx2", "ZicoSeq", "Wilcoxon_CLR_paired", "limma_voom"]
MLAB = {"LinDA": "LinDA", "ANCOMBC2": "ANCOM-BC2", "ALDEx2": "ALDEx2", "ZicoSeq": "ZicoSeq",
        "Wilcoxon_CLR_paired": "Wilcoxon (CLR)", "limma_voom": "limma-voom"}
PARTS = [("any_ige", "any"), ("single", "single"), ("multiple", "multiple"),
         ("high_score", "high"), ("medium_score", "medium"), ("low_score", "low")] + \
        [(f"cluster_{c}", c) for c in "ABCDEFGHI"]
ALPHA = 0.05

t = pd.read_csv(T / "da_all_methods_long.csv", low_memory=False)


def save(fig, stem):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=300, bbox_inches="tight",
                    transparent=False, facecolor="white")
    print("  wrote", stem)


# ── B: method heatmap ─────────────────────────────────────────────────────
def panel_B(ax_bar, ax):
    ev = t[t.q_within.notna()]
    rate = (ev.assign(sig=ev.q_within < ALPHA)
              .groupby(["partition", "method"]).sig.mean().mul(100)
              .unstack("method").reindex(index=[p for p, _ in PARTS], columns=METHODS))
    cmap = LinearSegmentedColormap.from_list("w2r", ["#ffffff", "#fcbba1", "#ef3b2c", "#67000d"])
    im = ax.imshow(rate.values, aspect="auto", cmap=cmap, vmin=0, vmax=max(40, rate.values.max()))
    for i in range(rate.shape[0]):
        for j in range(rate.shape[1]):
            v = rate.values[i, j]
            ax.text(j, i, f"{v:.1f}" if v >= 0.05 else "0", ha="center", va="center",
                    fontsize=8.5, color="white" if v > 22 else NAVY)
    ax.set_xticks(range(len(METHODS)))
    ax.set_xticklabels([MLAB[m] for m in METHODS], rotation=35, ha="right")
    ax.set_yticks(range(len(PARTS))); ax.set_yticklabels([lab for _, lab in PARTS])
    ax.set_ylabel("IgE group")
    for s in ax.spines.values():
        s.set_visible(False)
    overall = ev.assign(sig=ev.q_within < ALPHA).groupby("method").sig.mean().mul(100).reindex(METHODS)
    ax_bar.bar(range(len(METHODS)), overall.values, color=NAVY, width=0.7)
    for j, v in enumerate(overall.values):
        ax_bar.text(j, v + 0.6, f"{v:.1f}%" if v >= 0.05 else "0%", ha="center", va="bottom", fontsize=8.5)
    ax_bar.set_xlim(-0.5, len(METHODS) - 0.5); ax_bar.set_xticks([])
    ax_bar.set_ylim(0, overall.max() * 1.3); ax_bar.set_ylabel("% significant\n(all groups)")
    ax_bar.spines[["top", "right"]].set_visible(False)
    return im


# ── C: corrected volcano ──────────────────────────────────────────────────
def panel_C(ax):
    sp = t[(t.level == "species") & (t.prevalence == "prev05") & (t.partition == "any_ige")]
    nsig = sp[sp.q_within < ALPHA].groupby("taxon").method.nunique()
    li = sp[sp.method == "LinDA"].copy()
    li["nsig"] = li.taxon.map(nsig).fillna(0).astype(int)
    li["y"] = -np.log10(li.p_raw)
    cols = {0: "#c9c9c9", 1: "#fdae61", 2: "#d7301f", 3: "#7f0000"}
    li["c"] = li.nsig.clip(upper=3).map(cols)
    li = li.sort_values("nsig")
    ax.scatter(li.effect_log2, li.y, s=26, c=li.c, edgecolor=NAVY, linewidth=0.4, zorder=3)
    passed = li[li.q_within < ALPHA]
    if len(passed):
        pcut = passed.p_raw.max()
        ax.axhline(-np.log10(pcut), color=NAVY, ls="--", lw=0.9)
        ax.text(0.98, -np.log10(pcut) + 0.06, "LinDA BH q = 0.05", transform=ax.get_yaxis_transform(),
                fontsize=8.5, color=NAVY, va="bottom", ha="right")
    ax.axvline(0, color="#bbbbbb", lw=0.8, zorder=1)
    r = li[li.taxon == "Prevotella_copri"].iloc[0]
    ax.annotate("$\\it{P.\\ copri}$", (r.effect_log2, r.y), xytext=(8, -12),
                textcoords="offset points", fontsize=9.5, color=NAVY)
    ax.set_xlabel("log$_2$ fold change (LinDA)")
    ax.set_ylabel("$-$log$_{10}$ p (raw)")
    ax.spines[["top", "right"]].set_visible(False)
    handles = [Line2D([], [], marker="o", ls="", markersize=6, markerfacecolor=c,
                      markeredgecolor=NAVY, markeredgewidth=0.4,
                      label=(f"{k}" if k < 3 else "≥3") + " method" + ("" if k == 1 else "s"))
               for k, c in cols.items()]
    ax.legend(handles=handles, title="q < 0.05 in", frameon=False, loc="center right",
              title_fontsize=9)
    ax.set_title("All sensitized vs matched controls\n(species, 5% prevalence)", fontsize=9.5)
    return len(li)


# ── D: consensus-taxon stability ──────────────────────────────────────────
def panel_D(ax):
    base = t[(t.taxon == "Prevotella_copri") & (t.partition == "any_ige") & (t.level == "species")]
    se = pd.read_csv(T / "da_consensus_long.csv", low_memory=False,
                     usecols=["taxon", "method", "level", "prevalence", "partition", "se", "stat"])
    se = se[(se.taxon == "Prevotella_copri") & (se.partition == "any_ige") & (se.level == "species")]
    ext = pd.concat([pd.read_csv(DAEXT / f"daext_long_species_{pv}_any_ige.csv").assign(prevalence=pv)
                     for pv in ["prev10", "prev05", "prev01"]])
    ext = ext[ext.taxon == "Prevotella_copri"][["method", "prevalence", "se", "stat"]]
    sef = pd.concat([se[["method", "prevalence", "se", "stat"]], ext])
    d = base.merge(sef, on=["method", "prevalence"], how="left")
    lim = d.method == "limma_voom"
    d.loc[lim, "se"] = (d.loc[lim, "effect_log2"] / d.loc[lim, "stat"]).abs()  # SE = logFC / t
    PREV = [("prev10", "10%", -0.22), ("prev05", "5%", 0.0), ("prev01", "1%", 0.22)]
    PCOL = {"prev10": "#08519c", "prev05": "#3182bd", "prev01": "#9ecae1"}
    for i, m in enumerate(METHODS):
        for pv, plab, off in PREV:
            r = d[(d.method == m) & (d.prevalence == pv)]
            if not len(r):
                continue
            r = r.iloc[0]; y = i + off
            if pd.notna(r.se):
                ax.errorbar(r.effect_log2, y, xerr=1.96 * r.se, fmt="none",
                            ecolor=PCOL[pv], elinewidth=1.4, capsize=2.5, zorder=2)
            filled = r.q_within < ALPHA
            ax.plot(r.effect_log2, y, "o", ms=7, zorder=3,
                    markerfacecolor=PCOL[pv] if filled else "white",
                    markeredgecolor=PCOL[pv], markeredgewidth=1.6)
    ax.axvline(0, color=NAVY, lw=0.9, ls="--")
    ax.set_yticks(range(len(METHODS)))
    ax.set_yticklabels([MLAB[m] + ("$^{*}$" if m == "ZicoSeq" else "") for m in METHODS])
    ax.invert_yaxis()
    ax.set_xlabel("effect (log$_2$, method-specific scale)")
    ax.set_title("$\\it{P.\\ copri}$ aggregate across methods\nand prevalence filters", fontsize=9.5)
    ax.spines[["top", "right"]].set_visible(False)
    h = [Line2D([], [], marker="o", ls="-", color=PCOL[pv], markerfacecolor=PCOL[pv],
                markeredgecolor=PCOL[pv], label=f"{plab} prevalence") for pv, plab, _ in PREV]
    h += [Line2D([], [], marker="o", ls="", color="#555555", markerfacecolor="#555555",
                 label="q < 0.05"),
          Line2D([], [], marker="o", ls="", color="#555555", markerfacecolor="white",
                 label="q ≥ 0.05")]
    ax.legend(handles=h, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13),
              ncol=3, fontsize=8.5)


# ── individual panels ─────────────────────────────────────────────────────
print("=== Figure 3B-D ===")
fig = plt.figure(figsize=(5.2, 6.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 5], width_ratios=[20, 1], hspace=0.06, wspace=0.05)
ax_bar = fig.add_subplot(gs[0, 0]); ax = fig.add_subplot(gs[1, 0])
im = panel_B(ax_bar, ax)
cb = fig.colorbar(im, cax=fig.add_subplot(gs[1, 1])); cb.set_label("% of tested taxa with q < 0.05")
save(fig, "fig3B_method_heatmap"); plt.close(fig)

fig, ax = plt.subplots(figsize=(5.4, 4.6)); n = panel_C(ax)
save(fig, "fig3C_volcano_corrected"); plt.close(fig)

fig, ax = plt.subplots(figsize=(5.4, 4.6)); panel_D(ax)
save(fig, "fig3D_consensus_stability"); plt.close(fig)

# ── combined B-D, for placing beside the author's Figure 3A ───────────────
fig = plt.figure(figsize=(13.5, 6.9))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 5], width_ratios=[10, 11, 11],
                      hspace=0.08, wspace=0.55)
ax_bar = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[1, 0])
im = panel_B(ax_bar, axB)
# inset colorbar: does not shrink axB, so the bars above stay aligned with the columns
cb = fig.colorbar(im, cax=axB.inset_axes([1.04, 0.30, 0.045, 0.45])); cb.set_label("% significant\n(q < 0.05)")
axC = fig.add_subplot(gs[:, 1]); panel_C(axC)
axD = fig.add_subplot(gs[:, 2]); panel_D(axD)
for a, lab, xy in [(ax_bar, "B", (-0.22, 1.15)), (axC, "C", (-0.16, 1.10)), (axD, "D", (-0.30, 1.10))]:
    a.text(*xy, lab, transform=a.transAxes, fontsize=16, fontweight="bold", va="bottom")
save(fig, "fig3BCD"); plt.close(fig)
print(f"  volcano: {n} species tested (5% prevalence)")
