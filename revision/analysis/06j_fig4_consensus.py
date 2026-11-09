#!/usr/bin/env python
"""
Module 06j - revised Figure 4B-D: the consensus story in three simple panels.
Panel A (alpha diversity) is kept as in the manuscript.

Data: full cohort, covariate-adjusted design only (all sensitized vs all 233
IgE-negative, ~ IgE status + sex + age + BMI; 02f + --design full), five
methods: LinDA, ANCOM-BC2, ALDEx2, ZicoSeq, limma-voom. Significant = BH
q_within < 0.05 (within one analysis).

  B  Methods differ: % of all tests called significant per method, pooled over
     the 180 analyses (4 ranks x 3 prevalence filters x 15 IgE groups).
     Denominator = features tested (46,080 per method), the same for every
     method; the appendix heatmap shows the per-analysis breakdown.
  C  Agreement is rare: every taxon-analysis significant in >= 1 method,
     split by how many methods called it.
  D  One volcano per method at FAMILY level for the primary contrast (all
     sensitized vs controls, 5% prevalence filter): x = the method's own effect
     (log2 scale; LinDA/ANCOM-BC2/limma-voom fold change, ALDEx2/ZicoSeq CLR
     difference), y = -log10 raw p, dashed line = BH q = 0.05 for that method.
     The four families of the submitted paper and Prevotellaceae keep their
     original colours in every panel; black edge = q < 0.05.
     (An earlier consensus-map version of D was replaced: hard to read.)

Runs in scont_plots (cnsplots), Roboto font.
Output (results/figures/review/):
  fig4B_methods_differ.*  fig4C_agreement.*  fig4D_method_volcanoes.*  fig4BCD_consensus.*
  (each .png/.svg/.pdf)
Source data: results/tables/fig4BC_counts.csv, fig4D_family_volcano.csv,
             fig4_consensus_primary.csv
Usage:  python revision/analysis/06j_fig4_consensus.py
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
ALPHA = 0.05
METHODS = ["ALDEx2", "ZicoSeq", "LinDA", "limma_voom", "ANCOMBC2"]
MLAB = {"LinDA": "LinDA", "ANCOMBC2": "ANCOM-BC2", "ALDEx2": "ALDEx2",
        "ZicoSeq": "ZicoSeq", "limma_voom": "limma-voom"}
FAM = {"Lachnospiraceae": "#fcc853", "Bacteroidaceae": "#078a07",
       "Oscillospiraceae": "#1f77b4", "Veillonellaceae": "#ff7f0e",
       "Prevotellaceae": "#c2185b"}
NAVY, GREY = "#2a3f5f", "#c9ccd3"

# ── data ──────────────────────────────────────────────────────────────────
files = sorted((RES / "da_full").glob("da_long_*.csv")) + \
        sorted((RES / "da_extended_full").glob("daext_long_*.csv"))
cols = ["taxon", "method", "effect_log2", "q_within", "level", "prevalence", "partition"]
d = pd.concat([pd.read_csv(f, usecols=cols) for f in files], ignore_index=True)
d = d[d.method.isin(METHODS)]
d["sig"] = d.q_within < ALPHA
n_analyses = len(d[["level", "prevalence", "partition"]].drop_duplicates())
assert n_analyses == 180, n_analyses

# denominators = features tested per analysis, identical across methods
n_feat = d[d.method == "LinDA"].groupby(["level", "prevalence", "partition"]).size().sum()
rate = d.groupby("method").sig.sum().reindex(METHODS)
rate = pd.DataFrame({"n_sig": rate, "pct": 100 * rate / n_feat})

key = ["level", "prevalence", "partition", "taxon"]
g = d.groupby(key).agg(n_sig=("sig", "sum"),
                       up=("effect_log2", lambda e: int((e > 0).sum())),
                       dn=("effect_log2", lambda e: int((e < 0).sum())),
                       methods=("method", lambda m: ""),).reset_index()
sig_methods = d[d.sig].groupby(key).method.apply(lambda m: ";".join(sorted(m)))
g = g.drop(columns="methods").merge(sig_methods.rename("sig_methods").reset_index(),
                                    on=key, how="left")
g["agree"] = g.up - g.dn
g["consensus"] = (g.n_sig >= 2) & ((g.up == 0) | (g.dn == 0))
hits = g[g.n_sig >= 1]
agree_n = hits.n_sig.clip(upper=4).value_counts().reindex([1, 2, 3, 4], fill_value=0)

pd.concat([rate.assign(panel="B").rename_axis("item").reset_index(),
           pd.DataFrame({"item": ["1 method", "2 methods", "3 methods", "4-5 methods"],
                         "n_sig": agree_n.values,
                         "pct": 100 * agree_n.values / len(hits), "panel": "C"})]) \
  .to_csv(RES / "tables" / "fig4BC_counts.csv", index=False)

prim = g[(g.prevalence == "prev05") & (g.partition == "any_ige")].copy()
prim.to_csv(RES / "tables" / "fig4_consensus_primary.csv", index=False)

# ── style ─────────────────────────────────────────────────────────────────
cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "savefig.transparent": False,
                     "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "font.size": 9.5, "axes.labelsize": 10,
                     "xtick.labelsize": 9, "ytick.labelsize": 9.5, "svg.fonttype": "none"})


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)


def panel_b(ax):
    y = np.arange(len(METHODS))
    ax.barh(y, rate.pct, color=NAVY, height=0.62)
    for yi, (p, n) in enumerate(zip(rate.pct, rate.n_sig)):
        ax.text(p + 0.25, yi, f"{p:.1f}%" if p >= 0.1 or p == 0 else "<0.1%",
                va="center", fontsize=9)
    ax.set_yticks(y); ax.set_yticklabels([MLAB[m] for m in METHODS])
    ax.set_xlim(0, rate.pct.max() * 1.25)
    ax.set_xlabel("% of tests significant (q < 0.05)")
    ax.set_title("Methods differ", fontsize=10.5, color=NAVY, loc="left")
    clean(ax)


def panel_c(ax):
    labs = ["1 method", "2 methods", "3 methods", "4–5 methods"]
    y = np.arange(4)
    pct = 100 * agree_n.values / len(hits)
    ax.barh(y, pct, color=["#9aa5b8", NAVY, NAVY, NAVY], height=0.62)
    for yi, (p, n) in enumerate(zip(pct, agree_n.values)):
        txt = f"{p:.1f}%  (n = {n:,})" if p >= 0.1 or n == 0 else f"<0.1%  (n = {n})"
        ax.text(p + 1.5, yi, txt, va="center", fontsize=9)
    ax.set_yticks(y); ax.set_yticklabels(labs)
    ax.invert_yaxis()
    ax.set_xlim(0, 135)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel(f"% of significant calls (n = {len(hits):,})")
    ax.set_title("Methods rarely agree", fontsize=10.5, color=NAVY, loc="left")
    clean(ax)


VMETH = ["LinDA", "ANCOMBC2", "ALDEx2", "ZicoSeq", "limma_voom"]
fam_tab = pd.concat([pd.read_csv(RES / "da_full" / "da_long_family_prev05_any_ige.csv"),
                     pd.read_csv(RES / "da_extended_full" / "daext_long_family_prev05_any_ige.csv")])
fam_tab = fam_tab[fam_tab.method.isin(VMETH)]
fam_tab.to_csv(RES / "tables" / "fig4D_family_volcano.csv", index=False)


def panel_d(axes):
    """One family-level volcano per method, primary contrast, full cohort."""
    for k, (ax, m) in enumerate(zip(axes, VMETH)):
        t = fam_tab[fam_tab.method == m].copy()
        t["y"] = -np.log10(t.p_raw)
        sig = t.q_within < ALPHA
        hl = t.taxon.isin(FAM)
        rest = t[~hl]
        ax.scatter(rest.effect_log2, rest.y, s=14, lw=0, zorder=2,
                   c=np.where(rest.q_within < ALPHA, NAVY, GREY))
        for _, r in t[hl].iterrows():
            ax.scatter(r.effect_log2, r.y, s=58, c=FAM[r.taxon], zorder=4,
                       ec="black" if r.q_within < ALPHA else "white",
                       lw=1.1 if r.q_within < ALPHA else 0.6)
        if sig.any():   # BH q = 0.05 drawn at the largest raw p that passes
            ax.axhline(-np.log10(t.p_raw[sig].max()), ls="--", lw=0.8, color="#d62728", zorder=1)
        ax.axvline(0, color="#e0e0e0", lw=0.8, zorder=0)
        ax.set_title(f"{MLAB[m]}  ({int(sig.sum())} of {len(t)})", fontsize=10, color=NAVY)
        ax.set_xlabel("effect (log$_2$)")
        if k == 0:
            ax.set_ylabel("$-$log$_{10}$ p")
        lo, hi = ax.get_xlim(); r = max(abs(lo), abs(hi)); ax.set_xlim(-r, r)
        ax.set_ylim(0, max(t.y.max() * 1.08, 2.0))
        clean(ax)


def d_legend(fig, y):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", mfc=c, mec="white", ms=8, label=f)
         for f, c in FAM.items()]
    h += [Line2D([], [], marker="o", ls="", mfc=NAVY, mec=NAVY, ms=5, label="other family, q < 0.05"),
          Line2D([], [], marker="o", ls="", mfc=GREY, mec=GREY, ms=5, label="other family, n.s."),
          Line2D([], [], marker="o", ls="", mfc="white", mec="black", mew=1.1, ms=8,
                 label="highlighted, q < 0.05"),
          Line2D([], [], ls="--", color="#d62728", label="BH q = 0.05")]
    fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, y), ncol=5,
               frameon=False, fontsize=9, handletextpad=0.3, columnspacing=1.4)


def save(fig, stem):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=300, bbox_inches="tight",
                    transparent=False, facecolor="white")
    print("  wrote", stem)



OUT.mkdir(parents=True, exist_ok=True)
for fn, stem, size in [(panel_b, "fig4B_methods_differ", (4.2, 3.0)),
                       (panel_c, "fig4C_agreement", (4.2, 3.0))]:
    fig, ax = plt.subplots(figsize=size)
    fn(ax)
    save(fig, stem)
    plt.close(fig)

fig, axes = plt.subplots(1, 5, figsize=(15, 3.4), gridspec_kw={"wspace": 0.3})
panel_d(axes)
d_legend(fig, -0.08)
save(fig, "fig4D_method_volcanoes")
plt.close(fig)

fig = plt.figure(figsize=(15, 8.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 0.95], hspace=0.55)
top = gs[0].subgridspec(1, 2, wspace=0.55)
for i, (fn, tag) in enumerate([(panel_b, "B"), (panel_c, "C")]):
    ax = fig.add_subplot(top[0, i])
    fn(ax)
    ax.text(-0.25, 1.12, tag, transform=ax.transAxes, fontsize=16, fontweight="bold")
bottom = gs[1].subgridspec(1, 5, wspace=0.3)
axes = [fig.add_subplot(bottom[0, j]) for j in range(5)]
panel_d(axes)
axes[0].text(-0.42, 1.16, "D", transform=axes[0].transAxes, fontsize=16, fontweight="bold")
d_legend(fig, 0.035)
save(fig, "fig4BCD_consensus")

print(rate.round(2))
print("agreement:", dict(zip(["1", "2", "3", "4-5"], agree_n.values)), "of", len(hits))
print("primary contrast, consensus taxa:")
print(prim[prim.consensus][["level", "taxon", "n_sig", "agree", "sig_methods"]].to_string(index=False))
print(prim[(prim.level == "family") & prim.taxon.isin(FAM)][["taxon", "n_sig", "agree"]].to_string(index=False))
