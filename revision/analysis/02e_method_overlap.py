#!/usr/bin/env python
"""
Module 02e - merge all six DA methods and quantify their overlap.

This is the direct answer to "use more than one differential-abundance method,
and show whether they agree."

  module 02   LinDA, ANCOM-BC2, ALDEx2                  (log-ratio family)
  module 02d  ZicoSeq, Wilcoxon-CLR paired, limma-voom  (permutation /
                                                         non-parametric / count)

WHAT THIS DELIBERATELY DOES NOT DO
  It does not average effect sizes. The six are on different scales - log2 fold
  change of relative abundance, natural-log LFC, CLR difference, log2 CPM fold
  change - and averaging them would produce a number that means nothing.
  Agreement is measured two ways instead: shared significance calls, and rank
  correlation of effect sizes.

  It also compares every method at the SAME multiplicity procedure (BH within
  analysis, the `q_within` column). ZicoSeq additionally reports its own
  permutation FDR in `q_native`; using that for ZicoSeq while using BH for the
  others would make the overlap partly an artifact of the correction rather
  than of the test. `q_native` is reported alongside, never mixed in.

Usage:
  python revision/analysis/02e_method_overlap.py
"""

from __future__ import annotations

import itertools
import time
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

REVISION = Path(__file__).resolve().parent.parent
DA = REVISION / "results" / "da"
DAEXT = REVISION / "results" / "da_extended"
OUT = REVISION / "results" / "tables"
ALPHA = 0.05

LOGRATIO = {"LinDA", "ANCOMBC2", "ALDEx2", "ZicoSeq", "Wilcoxon_CLR_paired"}
COUNT = {"limma_voom"}


def load(folder: Path, pattern: str) -> pd.DataFrame:
    files = sorted(folder.glob(pattern))
    if not files:
        return pd.DataFrame()
    return pd.concat([pd.read_csv(f, low_memory=False) for f in files],
                     ignore_index=True)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    base = load(DA, "da_long_*.csv")
    ext = load(DAEXT, "daext_long_*.csv")
    if base.empty and ext.empty:
        print("No DA output found. Run modules 02 and 02d first.")
        return 1

    # ── Guard: is the input grid complete and settled? ───────────────────
    # An earlier run read results/da_extended while the array was still
    # rewriting it and reported ZicoSeq coverage of 26,693 where the settled
    # data gives 29,762. Nothing warned; the numbers simply differed. Silence is
    # not acceptable for a table that goes into a response letter.
    EXPECT_EXT, EXPECT_ASV = 180, 45
    ext_files = sorted(DAEXT.glob("daext_long_*.csv"))
    asv_files = sorted(DA.glob("da_long_asv_*.csv"))
    now = time.time()
    fresh = [f for f in ext_files + asv_files if now - f.stat().st_mtime < 120]
    problems = []
    if len(ext_files) != EXPECT_EXT:
        problems.append(f"extended DA: {len(ext_files)} files, expected {EXPECT_EXT}")
    if len(asv_files) != EXPECT_ASV:
        problems.append(f"ASV consensus DA: {len(asv_files)} files, expected {EXPECT_ASV}")
    if fresh:
        problems.append(f"{len(fresh)} file(s) modified in the last 2 min - "
                        "the grid may still be being written")
    if problems:
        print("\n  ⚠ INPUT GRID IS NOT SETTLED - numbers below are PROVISIONAL:")
        for p in problems:
            print(f"      - {p}")
        print()
    else:
        print(f"  input grid complete and settled: {len(ext_files)} extended "
              f"+ {len(asv_files)} ASV files")

    keep = ["taxon", "method", "effect_log2", "effect_type", "p_raw",
            "q_within", "level", "prevalence", "partition", "n_taxa"]
    for df in (base, ext):
        for c in keep:
            if not df.empty and c not in df.columns:
                df[c] = pd.NA
    parts = [d[keep] for d in (base, ext) if not d.empty]
    t = pd.concat(parts, ignore_index=True)
    t["key"] = (t.level.astype(str) + "|" + t.prevalence.astype(str) + "|"
                + t.partition.astype(str) + "|" + t.taxon.astype(str))

    methods = sorted(t.method.dropna().unique())
    print("=== Six-method DA overlap ===")
    print(f"  methods : {', '.join(methods)}")
    print(f"  tests   : {len(t):,}")
    print(f"  levels  : {sorted(t.level.dropna().unique())}")
    print(f"  analyses: {t.groupby(['level','prevalence','partition']).ngroups}")

    # ── Per-method significance rate ─────────────────────────────────────
    print(f"\n  {'method':22s} {'tests':>8s} {'sig':>7s} {'%':>7s} {'min q':>10s}")
    rows = []
    for m in methods:
        g = t[t.method == m]
        s = int((g.q_within < ALPHA).sum())
        print(f"  {m:22s} {len(g):8,} {s:7,} {100*s/max(len(g),1):6.2f}% "
              f"{g.q_within.min():10.3g}")
        rows.append({"method": m, "n_tests": len(g), "n_sig": s,
                     "pct_sig": 100 * s / max(len(g), 1),
                     "min_q": g.q_within.min()})
    pd.DataFrame(rows).to_csv(OUT / "da_method_rates.csv", index=False)

    # ── Overlap of significance calls ────────────────────────────────────
    # `evaluated` is NOT the same as `all keys`. ZicoSeq drops features that are
    # constant within a partition, and ANCOM-BC2 applies internal filtering
    # (at ASV rank it returned 898 of 1459 features where LinDA returned 1459).
    # A taxon a method never tested must not be counted as "that method found
    # it non-significant" - that would silently inflate apparent disagreement.
    # Every pairwise statistic below is therefore computed on the intersection
    # of what BOTH methods actually evaluated.
    evaluated = {m: set(t.loc[(t.method == m) & t.q_within.notna(), "key"])
                 for m in methods}
    sig = {m: set(t.loc[(t.method == m) & (t.q_within < ALPHA), "key"])
           for m in methods}

    print("\n  === Feature coverage (how much each method actually tested) ===")
    total_keys = t.key.nunique()
    for m in methods:
        print(f"    {m:22s} evaluated {len(evaluated[m]):7,} of {total_keys:,} "
              f"taxon-analyses ({100*len(evaluated[m])/max(total_keys,1):5.1f}%)")

    print("\n  === Pairwise overlap of significant calls (q_within < 0.05) ===")
    print("      restricted to taxon-analyses BOTH methods evaluated")
    pw = []
    for a, b in itertools.combinations(methods, 2):
        common = evaluated[a] & evaluated[b]
        sa, sb = sig[a] & common, sig[b] & common
        inter = len(sa & sb)
        union = len(sa | sb)
        j = inter / union if union else float("nan")
        if union:
            print(f"    {a:22s} & {b:22s}  shared {inter:4d} of {union:4d} "
                  f"(Jaccard {j:.3f})   common features {len(common):,}")
        pw.append({"method_a": a, "method_b": b, "n_common_evaluated": len(common),
                   "n_sig_a_in_common": len(sa), "n_sig_b_in_common": len(sb),
                   "n_shared": inter, "n_union": union, "jaccard": j})
    pd.DataFrame(pw).to_csv(OUT / "da_method_overlap_pairwise.csv", index=False)

    # ── How many methods called each taxon-analysis? (UpSet input) ───────
    allsig = sorted(set().union(*sig.values())) if sig else []
    upset = []
    for k in allsig:
        callers = sorted(m for m in methods if k in sig[m])
        # n_tested is the honest denominator for "how many methods agreed":
        # 2-of-3 that tested it is stronger evidence than 2-of-6 where three
        # methods never evaluated the feature at all.
        tested = sorted(m for m in methods if k in evaluated[m])
        lev, prev, part, tax = k.split("|", 3)
        upset.append({"level": lev, "prevalence": prev, "partition": part,
                      "taxon": tax, "n_methods": len(callers),
                      "n_tested": len(tested),
                      "methods": ",".join(callers),
                      "methods_not_tested": ",".join(
                          m for m in methods if m not in tested),
                      "n_logratio": sum(m in LOGRATIO for m in callers),
                      "n_count": sum(m in COUNT for m in callers)})
    up = pd.DataFrame(upset)
    if not up.empty:
        up = up.sort_values(["n_methods", "level", "partition"],
                            ascending=[False, True, True])
        up.to_csv(OUT / "da_method_upset.csv", index=False)
        print(f"\n  === Taxon-analyses called by >=1 method: {len(up):,} ===")
        print(up.n_methods.value_counts().sort_index(ascending=False)
              .rename("n_taxon_analyses").to_string())

        multi = up[up.n_methods >= 2]
        print(f"\n  === CONSENSUS (>=2 of {len(methods)} methods): {len(multi)} ===")
        if len(multi):
            print(multi.head(40).to_string(index=False))

        # The reviewer-facing question: does the consensus hold at each rank?
        print("\n  === Consensus by taxonomic rank ===")
        for lev in ["family", "genus", "species", "asv"]:
            sub = up[up.level == lev]
            if not len(sub):
                continue
            print(f"    {lev:8s} any-method {len(sub):5d}   >=2 methods "
                  f"{int((sub.n_methods >= 2).sum()):4d}   >=3 methods "
                  f"{int((sub.n_methods >= 3).sum()):4d}")

    # ── Effect-size agreement, the fairer measure ────────────────────────
    # Two methods can rank taxa near-identically yet share no significance
    # calls if one is simply more conservative. Reporting only the overlap of
    # calls would hide that, and would overstate genuine disagreement.
    w = t.pivot_table(index="key", columns="method", values="effect_log2",
                      aggfunc="first")
    # ⚠ Not every pair is an independent comparison. Three of the six report a
    # CLR-based difference, and for a BALANCED MATCHED design the mean of paired
    # differences is algebraically equal to the difference of means - so
    # ZicoSeq's and Wilcoxon's effect columns are the same number, and correlate
    # at rho = 1.000 by construction, not by agreement. ALDEx2's CLR difference
    # is the same quantity computed a third way.
    #
    # ZicoSeq's column is especially not a ZicoSeq statistic: its own R2 is
    # unsigned, so module 02d substitutes the CLR mean difference purely to
    # recover direction. Presenting any of these as cross-method validation
    # would be circular, so each pair is labelled with whether the two methods
    # estimate the SAME quantity.
    CLR_FAMILY = {"ALDEx2", "ZicoSeq", "Wilcoxon_CLR_paired"}

    def estimand(m: str) -> str:
        if m in CLR_FAMILY:
            return "CLR difference"
        return {"LinDA": "log2FC relative abundance",
                "ANCOMBC2": "log2FC absolute abundance",
                "limma_voom": "log2FC logCPM"}.get(m, "unknown")

    print("\n  === Effect-size rank correlation (Spearman, complete cases) ===")
    print("      'same estimand' pairs agree partly BY CONSTRUCTION - see note below")
    cor = []
    for a, b in itertools.combinations(methods, 2):
        if a not in w or b not in w:
            continue
        d = w[[a, b]].dropna()
        if len(d) < 30:
            continue
        r, _ = spearmanr(d[a], d[b])
        same = float((d[a] * d[b] > 0).mean() * 100)
        same_est = estimand(a) == estimand(b)
        flag = "  [SAME ESTIMAND - not independent]" if same_est else ""
        print(f"    {a:22s} & {b:22s}  rho {r:+.3f}  n={len(d):6,}  "
              f"same-sign {same:5.1f}%{flag}")
        cor.append({"method_a": a, "method_b": b, "spearman": r,
                    "n": len(d), "pct_same_sign": same,
                    "estimand_a": estimand(a), "estimand_b": estimand(b),
                    "same_estimand": same_est,
                    "independent_comparison": not same_est})

    ind = [c for c in cor if c["independent_comparison"]]
    if ind:
        print("\n    Independent comparisons only (different estimands):")
        for c in sorted(ind, key=lambda x: -x["spearman"]):
            print(f"      {c['method_a']:22s} & {c['method_b']:22s}  "
                  f"rho {c['spearman']:+.3f}  same-sign {c['pct_same_sign']:5.1f}%")
    pd.DataFrame(cor).to_csv(OUT / "da_method_effect_correlation.csv",
                             index=False)

    t.drop(columns=["key"]).to_csv(OUT / "da_all_methods_long.csv", index=False)
    print(f"\nWrote {OUT}/da_method_rates.csv")
    print(f"Wrote {OUT}/da_method_overlap_pairwise.csv")
    print(f"Wrote {OUT}/da_method_upset.csv")
    print(f"Wrote {OUT}/da_method_effect_correlation.csv")
    print(f"Wrote {OUT}/da_all_methods_long.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
