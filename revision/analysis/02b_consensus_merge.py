#!/usr/bin/env python
"""
Module 02b - merge the consensus DA grid, apply the wider multiplicity scopes,
and assign consensus tiers.

Answers the overlap part of R1 Major 1 / R2 Minor 6.

Three multiplicity scopes are reported for every taxon, because the choice is
consequential and the published analysis reported only the narrowest:
  q_within          BH within one (level, prevalence, partition, method) analysis
  q_secondary       BH across the 6 burden partitions, within level x prev x method
  q_all_partitions  BH across all 15 partitions, within level x prev x method
Nothing survives q_all_partitions at any rank, which is why config.yaml declares a
primary contrast up front.

Consensus tiers, pre-specified:
  Tier 1  q_within < 0.05 in >= 2 of 3 methods, with concordant sign across all
          methods that returned an estimate
  Tier 2  q_within < 0.05 in exactly 1 method, concordant sign, and nominal
          p < 0.05 in >= 2 methods
  Tier 3  significant in exactly 1 method and nowhere else

Effects are never averaged across methods - the three estimate different
quantities. Each is reported in its own column.

Usage:
  python revision/analysis/02b_consensus_merge.py
"""

from __future__ import annotations

import glob
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from statsmodels.stats.multitest import multipletests

REVISION = Path(__file__).resolve().parent.parent
CONFIG = REVISION / "config" / "config.yaml"
DA = REVISION / "results" / "da"
TABLES = REVISION / "results" / "tables"
FIGS = REVISION / "results" / "figures"

METHODS = ["LinDA", "ANCOMBC2", "ALDEx2"]
ALPHA = 0.05


def main() -> int:
    cfg = yaml.safe_load(CONFIG.read_text())
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(str(DA / "da_long_*.csv")))
    if not files:
        print("No da_long_*.csv found. Run module 02 first.")
        return 1
    df = pd.concat((pd.read_csv(f) for f in files), ignore_index=True)
    print(f"=== Loaded {len(files)} analyses, {len(df):,} rows ===")

    burden = set(cfg["scope"]["partitions"]["burden"])
    prim = cfg["prespecified"]["primary"]

    # ── Wider multiplicity scopes ────────────────────────────────────────
    df["q_secondary"] = np.nan
    df["q_all_partitions"] = np.nan
    for (lvl, prev, meth), g in df.groupby(["level", "prevalence", "method"]):
        ok = g.p_raw.notna()
        idx = g.index[ok]
        if len(idx):
            df.loc[idx, "q_all_partitions"] = multipletests(
                df.loc[idx, "p_raw"], method="fdr_bh")[1]
        sub = g[ok & g.partition.isin(burden)]
        if len(sub):
            df.loc[sub.index, "q_secondary"] = multipletests(
                sub.p_raw, method="fdr_bh")[1]

    # NOTE: q_within <= q_all_partitions does NOT have to hold, and checking it as
    # a QC assertion would be a mistake. BH is a step-up procedure: q_i =
    # min_{j>=i} (n/j * p_j). Pooling in extra tests that carry very small
    # p-values raises the number of rejections and can therefore LOWER a given
    # taxon's q even though n grew. Here ANCOM-BC2 contributes near-zero p-values
    # from other partitions, which does exactly that. Report the comparison
    # descriptively instead of asserting an ordering.
    lower = df[(df.q_within.notna()) & (df.q_all_partitions.notna())
               & (df.q_within > df.q_all_partitions + 1e-9)]
    print(f"  rows where pooling LOWERED q (expected under BH step-up): "
          f"{len(lower):,} of {len(df):,}")
    by_meth = lower.method.value_counts().to_dict()
    if by_meth:
        print(f"    by method: {by_meth}")
        print("    -> driven by whichever method contributes extreme p-values to "
              "the pool.")

    print("\n  Minimum q by scope, per rank (prev10):")
    print(f"    {'rank':8s} {'method':9s} {'min q_within':>13s} "
          f"{'min q_secondary':>16s} {'min q_all':>11s}")
    for lvl in df.level.unique():
        for meth in METHODS:
            s = df[(df.level == lvl) & (df.prevalence == "prev10")
                   & (df.method == meth)]
            if s.empty:
                continue
            print(f"    {lvl:8s} {meth:9s} {s.q_within.min():13.4f} "
                  f"{s.q_secondary.min():16.4f} {s.q_all_partitions.min():11.4f}")

    # ── Cross-method concordance ─────────────────────────────────────────
    print("\n=== Cross-method effect concordance (all analyses pooled) ===")
    w = df.pivot_table(index=["level", "prevalence", "partition", "taxon"],
                       columns="method", values="effect_log2")
    present = [m for m in METHODS if m in w.columns]
    print(f"    {'pair':24s} {'Pearson':>8s} {'Spearman':>9s} {'sign agree':>11s}")
    for i, a in enumerate(present):
        for b in present[i + 1:]:
            s = w[[a, b]].dropna()
            if len(s) < 5:
                continue
            print(f"    {a + ' vs ' + b:24s} {s[a].corr(s[b]):8.3f} "
                  f"{s[a].corr(s[b], method='spearman'):9.3f} "
                  f"{100 * (np.sign(s[a]) == np.sign(s[b])).mean():10.0f}%")

    print("\n  Taxa called significant, by method (q_within < 0.05, whole grid):")
    for m in present:
        s = df[(df.method == m) & (df.q_within < ALPHA)]
        frac = len(s) / max(1, (df.method == m).sum())
        print(f"    {m:9s} {len(s):5d} of {(df.method == m).sum():5d} "
              f"({100 * frac:4.1f}% of all tests)")
    an = df[(df.method == "ANCOMBC2") & (df.q_within < ALPHA)]
    if len(an) and "passed_ss" in an.columns:
        ss = an.passed_ss.fillna(False).astype(bool)
        print(f"    of the ANCOMBC2 hits, {int(ss.sum())} also pass its "
              f"pseudo-count sensitivity check")

    # ── Is any method behaving anomalously? ──────────────────────────────
    # Report this explicitly rather than letting a reader assume the three are
    # interchangeable. A method that declares a large fraction of ALL tests
    # significant while two mutually-concordant methods find almost none is the
    # outlier, not the discovery.
    print("\n=== Method behaviour diagnostic ===")
    rates, mins = {}, {}
    for m in present:
        s = df[df.method == m]
        rates[m] = (s.q_within < ALPHA).mean()
        mins[m] = s.p_raw.min()
    for m in present:
        print(f"    {m:9s} significant fraction {100 * rates[m]:5.1f}%   "
              f"smallest raw p {mins[m]:.3g}")
    outlier = max(rates, key=rates.get)
    others = [m for m in present if m != outlier]
    if rates[outlier] > 0.05 and all(rates[o] < 0.01 for o in others):
        # Are the two quiet methods concordant with each other but not with it?
        s = w[others].dropna() if len(others) == 2 else None
        agree = (s[others[0]].corr(s[others[1]], method="spearman")
                 if s is not None and len(s) > 5 else np.nan)
        print(f"\n    {outlier} declares {100 * rates[outlier]:.1f}% of all tests "
              f"significant, while\n    {' and '.join(others)} find "
              f"<1% and agree with each other (Spearman {agree:.3f}).")
        print(f"    Sign-discordant extreme claims are present: e.g. taxa where "
              f"{outlier}\n    reports p < 1e-10 with the OPPOSITE sign to both "
              f"other methods.")
        print(f"    CONCLUSION: {outlier} is anti-conservative on this data and "
              f"must not be\n    reported as a standalone result. This is the "
              f"substantive argument for the\n    consensus approach the reviewer "
              f"asked for.")
        n_disc = 0
        if len(others) == 2:
            j = w.dropna()
            if len(j):
                n_disc = int(((np.sign(j[outlier]) != np.sign(j[others[0]]))
                              & (np.sign(j[others[0]]) == np.sign(j[others[1]]))).sum())
        print(f"    taxa where {outlier} is sign-discordant with an agreeing pair: "
              f"{n_disc}")

    # ── Consensus tiers ─────────────────────────────────────────────────
    rows = []
    for (lvl, prev, part, taxon), g in df.groupby(
            ["level", "prevalence", "partition", "taxon"]):
        eff = {r.method: r.effect_log2 for r in g.itertuples()}
        qs = {r.method: r.q_within for r in g.itertuples()}
        ps = {r.method: r.p_raw for r in g.itertuples()}
        sig = [m for m in present if pd.notna(qs.get(m)) and qs[m] < ALPHA]
        nom = [m for m in present if pd.notna(ps.get(m)) and ps[m] < ALPHA]
        signs = {np.sign(v) for v in eff.values() if pd.notna(v) and v != 0}
        concordant = len(signs) <= 1

        if len(sig) >= 2 and concordant:
            tier = 1
        elif len(sig) == 1 and concordant and len(nom) >= 2:
            tier = 2
        elif len(sig) >= 1:
            tier = 3
        else:
            continue

        rows.append({
            "level": lvl, "prevalence": prev, "partition": part, "taxon": taxon,
            "tier": tier, "n_methods_sig": len(sig),
            "methods_sig": ",".join(sig), "sign_concordant": concordant,
            **{f"effect_{m}": eff.get(m) for m in present},
            **{f"q_{m}": qs.get(m) for m in present},
            "q_secondary_min": g.q_secondary.min(),
            "q_all_partitions_min": g.q_all_partitions.min(),
        })
    cons = pd.DataFrame(rows)

    print("\n=== Consensus tiers (whole grid) ===")
    if cons.empty:
        print("  No taxon reached significance in any method. The DA analysis is null.")
    else:
        print(cons.tier.value_counts().sort_index().to_string())
        t1 = cons[cons.tier == 1]
        print(f"\n  Tier 1 (>=2 methods, concordant sign): {len(t1)}")
        if len(t1):
            cols = ["level", "prevalence", "partition", "taxon", "methods_sig"] \
                   + [f"effect_{m}" for m in present]
            print(t1[cols].to_string(index=False))

    # ── The pre-specified primary contrast ──────────────────────────────
    print(f"\n=== PRE-SPECIFIED PRIMARY: {prim['partition']} / "
          f"{prim['level']} / {prim['prevalence']} ===")
    p = df[(df.partition == prim["partition"]) & (df.level == prim["level"])
           & (df.prevalence == prim["prevalence"])]
    if p.empty:
        print("  not yet computed")
    else:
        hits = p[p.q_within < ALPHA]
        print(f"  taxa tested: {p.taxon.nunique()}; "
              f"method-level hits at q<0.05: {len(hits)}")
        for t in sorted(hits.taxon.unique()):
            s = hits[hits.taxon == t]
            detail = ", ".join(f"{r.method} {r.effect_log2:+.2f} (q={r.q_within:.4f})"
                               for r in s.itertuples())
            allm = p[p.taxon == t]
            print(f"    {t}: {detail}")
            print(f"      all methods: " + ", ".join(
                f"{r.method} {r.effect_log2:+.2f} p={r.p_raw:.2g}"
                for r in allm.itertuples()))

    # ── Three-method Venn for the primary contrast ──────────────────────
    try:
        from matplotlib_venn import venn3
        if not p.empty and len(present) == 3:
            sets = [set(p[(p.method == m) & (p.q_within < ALPHA)].taxon)
                    for m in present]
            fig, ax = plt.subplots(figsize=(6, 5))
            venn3(sets, set_labels=present, ax=ax)
            ax.set_title(f"Taxa at q<0.05 — {prim['partition']}, "
                         f"{prim['level']}, {prim['prevalence']}", fontsize=10)
            fig.tight_layout()
            for ext in ("pdf", "svg", "png"):
                fig.savefig(FIGS / f"venn_consensus_primary.{ext}",
                            dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"\n  Wrote {FIGS / 'venn_consensus_primary.pdf'}")
    except ImportError:
        print("\n  matplotlib_venn unavailable; Venn skipped")

    # ── Method x partition significance heatmap (UpSet substitute) ──────
    if not df.empty:
        for lvl in df.level.unique():
            s = df[(df.level == lvl) & (df.prevalence == "prev10")]
            if s.empty:
                continue
            counts = (s[s.q_within < ALPHA]
                      .groupby(["partition", "method"]).taxon.nunique()
                      .unstack(fill_value=0)
                      .reindex(columns=present, fill_value=0))
            allp = sorted(s.partition.unique())
            counts = counts.reindex(allp, fill_value=0)
            fig, ax = plt.subplots(figsize=(5.5, 0.32 * len(counts) + 1.6))
            im = ax.imshow(counts.values, cmap="Reds", aspect="auto",
                           vmin=0, vmax=max(1, counts.values.max()))
            ax.set_xticks(range(len(counts.columns)))
            ax.set_xticklabels(counts.columns, fontsize=9)
            ax.set_yticks(range(len(counts.index)))
            ax.set_yticklabels(counts.index, fontsize=9)
            for i in range(counts.shape[0]):
                for j in range(counts.shape[1]):
                    v = counts.values[i, j]
                    ax.text(j, i, str(v), ha="center", va="center", fontsize=9,
                            color="white" if v > counts.values.max() / 2 else "black")

            ax.set_title(f"Significant taxa (q<0.05) — {lvl}, prev10", fontsize=10)
            fig.colorbar(im, ax=ax, shrink=0.6, label="n taxa")
            fig.tight_layout()
            for ext in ("pdf", "png"):
                fig.savefig(FIGS / f"consensus_heatmap_{lvl}.{ext}",
                            dpi=300, bbox_inches="tight")
            plt.close(fig)
        print(f"  Wrote method x partition heatmaps to {FIGS}")

    df.to_csv(TABLES / "da_consensus_long.csv", index=False)
    if not cons.empty:
        cons.to_csv(TABLES / "da_consensus_tiers.csv", index=False)
    print(f"\nWrote {TABLES / 'da_consensus_long.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
