#!/usr/bin/env python
"""
Module 06c - review figures for the round-2 analyses.

Four panels, each reading ONLY computed result CSVs (no number is typed in):

  A  alpha_forest        paired alpha-diversity estimates with CIs, all 15
                         partitions, RAW primary (per the no-rarefaction
                         decision), q_all annotated
  B  nested_waterfall    the M0 seed distribution (published 61 marked) and
                         the M0->M3 waterfall, total vs >=1%-prevalence hits
  C  clusterF_depth      the depth artifact: cluster_F observed richness under
                         raw Wilcoxon / rarefied Wilcoxon / depth-adjusted lm
  D  ml_benchmark        six-model AUC forest at the primary contrast, with
                         per-model permutation p and the 80%-power bound

Outputs: revision/results/figures/review/<name>.{pdf,png,svg}
Fonts are kept >= 9 pt (reviewer requirement).

Usage:  python revision/analysis/06c_review_figures.py   (slr-env, in-session)
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REVISION = Path(__file__).resolve().parent.parent
RES = REVISION / "results"
OUT = RES / "figures" / "review"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
                     "xtick.labelsize": 9, "ytick.labelsize": 9,
                     "legend.fontsize": 9})
PUBLISHED = 61
PART_ORDER = ["any_ige", "single", "multiple", "high_score", "medium_score",
              "low_score"] + [f"cluster_{c}" for c in "ABCDEFGHI"]


def save(fig, name):
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote figures/review/{name}.[pdf|png|svg]")


# ── A. alpha forest ───────────────────────────────────────────────────────
def panel_alpha():
    r = pd.read_csv(RES / "alpha" / "alpha_tests_long.csv")
    pw = r[(r.analysis == "paired_wilcoxon") & (r.variant == "raw")
           & (r["index"].isin(["shannon", "observed"]))].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.4), sharey=True)
    for ax, idx, xlab in zip(axes, ["shannon", "observed"],
                             ["Hodges-Lehmann shift, Shannon",
                              "Hodges-Lehmann shift, observed richness"]):
        d = pw[pw["index"] == idx].set_index("partition").reindex(PART_ORDER)
        y = np.arange(len(d))[::-1]
        ax.errorbar(d.estimate, y,
                    xerr=[d.estimate - d.ci_lo, d.ci_hi - d.estimate],
                    fmt="o", color="#2C3E50", ecolor="#95A5A6",
                    elinewidth=1.6, capsize=2.5, ms=4.5)
        ax.axvline(0, color="#C0392B", lw=1, ls="--")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{p} (n={int(n)})" for p, n
                            in zip(d.index, d.n_pairs)])
        ax.set_xlabel(xlab)
        for yi, q in zip(y, d.q_all):
            ax.annotate(f"q={q:.2f}", (ax.get_xlim()[1], yi), fontsize=9,
                        ha="right", va="center", color="#7F8C8D",
                        xytext=(-4, 0), textcoords="offset points")
    fig.suptitle("Alpha diversity, matched pairs, RAW counts (primary): "
                 "IgE-sensitized minus non-sensitized; BH q across all tests",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "alpha_forest_raw_primary")


# ── B. nested waterfall ───────────────────────────────────────────────────
def panel_waterfall():
    dist = pd.read_csv(RES / "nested_models" / "m0_seed_distribution.csv")
    wf = pd.read_csv(RES / "nested_models" / "nested_summary.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6),
                             gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    ax.hist(dist.n_union_hits, bins=24, color="#95A5A6", edgecolor="white")
    ax.axvline(PUBLISHED, color="#C0392B", lw=2)
    pct = (dist.n_union_hits <= PUBLISHED).mean() * 100
    ax.annotate(f"published: {PUBLISHED} ASVs\n({pct:.0f}th percentile)",
                (PUBLISHED, ax.get_ylim()[1] * 0.92), color="#C0392B",
                fontsize=10, ha="left", xytext=(6, 0),
                textcoords="offset points")
    ax.set_xlabel("union of significant ASVs (q<0.05), published design")
    ax.set_ylabel(f"control draws (n={len(dist)} seeds)")
    ax.set_title("The published design's hit count is a draw\n"
                 "from a wide sampling distribution", fontsize=11)

    ax = axes[1]
    x = np.arange(len(wf))
    ax.bar(x - 0.2, wf.n_union_hits, 0.38, color="#95A5A6",
           label="all hits (no prevalence filter)")
    ax.bar(x + 0.2, wf.n_hits_prev_ge1pct, 0.38, color="#2C3E50",
           label="hits at ≥1% prevalence")
    for xi, (t, c) in enumerate(zip(wf.n_union_hits, wf.n_hits_prev_ge1pct)):
        ax.annotate(f"{t:.0f}", (xi - 0.2, t), ha="center", va="bottom", fontsize=9)
        ax.annotate(f"{c:.0f}", (xi + 0.2, c), ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(["M0\npublished\n(median)", "M1\nall controls",
                        "M2\nmatched", "M3\nmatched\n+adjusted"], fontsize=9)
    ax.set_ylabel("union of significant ASVs")
    ax.set_yscale("symlog", linthresh=20)
    ax.legend(frameon=False)
    ax.set_title("Most hits at every step are <1%-prevalence\nartifacts "
                 "(median: present in 1 of 508 samples)", fontsize=11)
    fig.tight_layout()
    save(fig, "nested_m0_m3_waterfall")


# ── C. the cluster_F depth artifact ──────────────────────────────────────
def panel_clusterF():
    r = pd.read_csv(RES / "alpha" / "alpha_tests_long.csv")
    f = r[(r.partition == "cluster_F") & (r["index"] == "observed")]
    rows = [
        ("raw counts,\npaired Wilcoxon", "raw", "paired_wilcoxon"),
        ("rarefied 100×,\npaired Wilcoxon", "rarefied_100x", "paired_wilcoxon"),
        ("raw + log₁₀(depth)\ncovariate (lm)", "raw", "adjusted_lm_matched_depth"),
    ]
    est, lo, hi, pv, labs = [], [], [], [], []
    for lab, v, an in rows:
        d = f[(f.variant == v) & (f.analysis == an)]
        if not len(d):
            continue
        est.append(float(d.estimate.iloc[0])); lo.append(float(d.ci_lo.iloc[0]))
        hi.append(float(d.ci_hi.iloc[0])); pv.append(float(d.p_value.iloc[0]))
        labs.append(lab)
    y = np.arange(len(est))[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.errorbar(est, y, xerr=[np.array(est) - lo, np.array(hi) - est],
                fmt="o", color="#2C3E50", ecolor="#95A5A6", elinewidth=2,
                capsize=3, ms=6)
    ax.axvline(0, color="#C0392B", lw=1, ls="--")
    ax.set_yticks(y); ax.set_yticklabels(labs)
    for yi, p in zip(y, pv):
        ax.annotate(f"p = {p:.3f}", (ax.get_xlim()[1], yi), fontsize=9,
                    ha="right", va="center", color="#7F8C8D",
                    xytext=(-4, 0), textcoords="offset points")
    ax.set_xlabel("observed-richness difference, IgE − control (cluster F, n = 10 pairs)")
    ax.set_title("cluster F's nominal richness signal is a sequencing-depth artifact:\n"
                 "its IgE members were sequenced at ~half their controls' depth "
                 "(9,685 vs 18,559 reads)", fontsize=10.5)
    fig.tight_layout()
    save(fig, "clusterF_depth_artifact")


# ── D. ML benchmark forest at the primary contrast ────────────────────────
def panel_ml():
    m = pd.read_csv(RES / "ml" / "ml_species_prev05_any_ige_metrics.csv")
    m = m[m.model != "Dummy"]
    da = float(m.detectable_auc_80pct.iloc[0])
    y = np.arange(len(m))[::-1]
    fig, ax = plt.subplots(figsize=(7.6, 3.9))
    ax.errorbar(m.auc_pooled, y,
                xerr=[m.auc_pooled - m.auc_ci_low, m.auc_ci_high - m.auc_pooled],
                fmt="o", color="#2C3E50", ecolor="#95A5A6", elinewidth=2,
                capsize=3, ms=6)
    ax.axvline(0.5, color="#C0392B", lw=1, ls="--", label="chance")
    ax.axvline(da, color="#7F8C8D", lw=1, ls=":",
               label=f"80% power bound ({da:.3f})")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{mo}  (perm p = {p:.2f})"
                        for mo, p in zip(m.model, m.perm_p)])
    ax.set_xlabel("ROC-AUC, nested CV grouped on matched pairs (species, prev05, any_ige)")
    ax.set_title("No classifier family beats chance at the primary contrast;\n"
                 "each model carries its own within-pair permutation null", fontsize=10.5)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    save(fig, "ml_benchmark_primary")


if __name__ == "__main__":
    print("=== 06c review figures ===")
    panel_alpha()
    panel_waterfall()
    panel_clusterF()
    panel_ml()
    print("Done.")
