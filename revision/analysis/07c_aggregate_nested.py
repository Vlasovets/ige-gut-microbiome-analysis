#!/usr/bin/env python
"""
Module 07c - aggregate the M0 -> M3 nested-model chain into the waterfall.

Reads the per-task outputs of 07b and answers, with receipts:
  1. How many ASV hits does the PUBLISHED design (M0) produce, as a
     distribution over 100 down-sampling seeds - and where does the published
     61 sit in it?
  2. How does the hit count fall as each design element is fixed:
     M0 (down-sampled, unadjusted) -> M1 (all controls) -> M2 (matched)
     -> M3 (matched + sex/age/BMI)?
  3. Which of the seed-42 M0 hits survive each later step (the hit ledger)?

Definitions locked to the published counting rule: a "hit" is an ASV with
BH-within-partition q < 0.05; the headline number is the UNION of hits across
the 15 partitions. M0 seed tasks skip any_ige (it is seed-invariant: the
published design used all 233 controls there), so every seed's union merges in
the M1 any_ige hit set.

Usage:  python revision/analysis/07c_aggregate_nested.py   (slr-env, in-session)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REVISION = Path(__file__).resolve().parent.parent
NM = REVISION / "results" / "nested_models"
OUT = REVISION / "results" / "tables"
PUBLISHED = 61
SEEDS = list(range(42, 142))
PREV_CUT = 0.01     # the stratification boundary; also the smallest cutoff in
                    # the revision's prev01/05/10 sensitivity grid


def _prevalence() -> pd.Series:
    """ASV prevalence across the 508 analytic samples (share of samples > 0)."""
    import yaml
    cfg = yaml.safe_load((REVISION / "config" / "config.yaml").read_text())
    kora = Path(cfg["paths"]["kora_data"])
    ft = pd.read_csv(kora / "feature_table.tsv", sep="\t", index_col=0)
    ige = pd.read_csv(kora / "ige_kora_subset_all.csv", index_col=0)
    cols = [c for c in ft.columns if c in set(map(str, ige.index))]
    return (ft[cols] > 0).mean(axis=1)


PREV = None


def hits(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["partition", "asv"])
    df = pd.read_csv(path)
    return df if len(df) else pd.DataFrame(columns=df.columns)


def n_common(asvs) -> int:
    """Hits at >= PREV_CUT prevalence - the scale on which biology is claimed.

    WHY THIS STRATIFICATION EXISTS: measured on the computed results, 84% of
    the reproduced published-design hits (M0 seed 42), 93% of the M1 explosion
    and 100% of the M2 any_ige hits are ASVs below 1% prevalence - features
    present in a handful of samples, median ONE of 508, with |lfc| ~ 0.02,
    where LinDA's normal approximation has no support. The published pipeline
    applied NO prevalence filter (utils/preprocessing freq_threshold = 0), so
    its hit count is dominated by exactly this artifact class. Reporting the
    waterfall without the split would launder those artifacts into a design
    story.
    """
    return int((PREV.reindex(list(asvs)) >= PREV_CUT).sum())


def main() -> int:
    global PREV
    PREV = _prevalence()
    m1 = hits(NM / "M1" / "M1_hits.csv")
    m2 = hits(NM / "M2" / "M2_hits.csv")
    m3 = hits(NM / "M3" / "M3_hits.csv")
    m1_any = set(m1.loc[m1.partition == "any_ige", "asv"]) if len(m1) else set()

    # ── M0 distribution over seeds ────────────────────────────────────────
    rows, missing = [], []
    for s in SEEDS:
        f = NM / "M0" / f"M0_seed{s}_hits.csv"
        if not f.exists():
            missing.append(s)
            continue
        h = hits(f)
        union = set(h.asv) | m1_any            # merge the seed-invariant any_ige
        rows.append({"seed": s, "n_union_hits": len(union),
                     "n_union_hits_prev1pct": n_common(union),
                     "n_hits_excl_any_ige": h.asv.nunique(),
                     "n_partitions_with_hits": h.partition.nunique()})
    if missing:
        print(f"⚠ {len(missing)} M0 seed file(s) missing (array incomplete?): "
              f"{missing[:8]}{' ...' if len(missing) > 8 else ''}")
    m0 = pd.DataFrame(rows)
    if m0.empty or not len(m1):
        print("M0/M1 outputs incomplete - run after the 07b array finishes.")
        return 1
    m0.to_csv(NM / "m0_seed_distribution.csv", index=False)

    q = m0.n_union_hits.quantile
    print("=== M0 (published design), union hit count over "
          f"{len(m0)} seeds ===")
    print(f"  median {m0.n_union_hits.median():.0f}  IQR [{q(.25):.0f}, {q(.75):.0f}]"
          f"  range [{m0.n_union_hits.min()}, {m0.n_union_hits.max()}]")
    pct = (m0.n_union_hits <= PUBLISHED).mean() * 100
    seed42 = m0.loc[m0.seed == 42, "n_union_hits"]
    print(f"  published value {PUBLISHED} sits at the {pct:.0f}th percentile; "
          f"seed 42 here gives {int(seed42.iloc[0]) if len(seed42) else 'NA'}")

    # ── The waterfall, stratified by prevalence ───────────────────────────
    u1, u2, u3 = set(m1.asv), set(m2.asv), set(m3.asv)
    qp = m0.n_union_hits_prev1pct.quantile
    wf = pd.DataFrame([
        {"model": "M0", "design": "published: per-partition down-sampled controls (seed median), unadjusted",
         "n_union_hits": float(m0.n_union_hits.median()),
         "n_union_hits_iqr": f"[{q(.25):.0f}, {q(.75):.0f}]",
         "n_hits_prev_ge1pct": float(m0.n_union_hits_prev1pct.median()),
         "prev_ge1pct_iqr": f"[{qp(.25):.0f}, {qp(.75):.0f}]"},
        {"model": "M1", "design": "all 233 controls in every partition, unadjusted",
         "n_union_hits": len(u1), "n_union_hits_iqr": "",
         "n_hits_prev_ge1pct": n_common(u1), "prev_ge1pct_iqr": ""},
        {"model": "M2", "design": "matched pairs (sex exact, age ±5, BMI ±4), unadjusted",
         "n_union_hits": len(u2), "n_union_hits_iqr": "",
         "n_hits_prev_ge1pct": n_common(u2), "prev_ge1pct_iqr": ""},
        {"model": "M3", "design": "matched pairs + sex/age/BMI in the model",
         "n_union_hits": len(u3), "n_union_hits_iqr": "",
         "n_hits_prev_ge1pct": n_common(u3), "prev_ge1pct_iqr": ""},
    ])
    wf.to_csv(NM / "nested_summary.csv", index=False)
    print("\n=== waterfall ===")
    print(wf.to_string(index=False))

    # ── Hit ledger: what happens to the seed-42 M0 hits ───────────────────
    f42 = NM / "M0" / "M0_seed42_hits.csv"
    if f42.exists():
        h42 = set(hits(f42).asv) | m1_any
        ledger = pd.DataFrame({"asv": sorted(h42)})
        ledger["prevalence_508"] = PREV.reindex(ledger.asv).values
        ledger["in_M1"] = ledger.asv.isin(set(m1.asv))
        ledger["in_M2"] = ledger.asv.isin(set(m2.asv))
        ledger["in_M3"] = ledger.asv.isin(set(m3.asv))
        ledger.to_csv(NM / "hit_waterfall.csv", index=False)
        print(f"\n=== seed-42 hit ledger ({len(ledger)} ASVs) ===")
        print(f"  survive M1: {ledger.in_M1.sum()}  |  survive M2: "
              f"{ledger.in_M2.sum()}  |  survive M3: {ledger.in_M3.sum()}")
        if len(m3):
            print(f"  M3 hits: {', '.join(sorted(set(m3.asv))[:6])}")

    # ── QC gate ───────────────────────────────────────────────────────────
    lo, hi = m0.n_union_hits.quantile([0.05, 0.95])
    print("\n=== QC gate: does the M0 distribution cover the published 61? ===")
    if lo <= PUBLISHED <= hi:
        print(f"  PASS: 61 lies inside the seed 5-95% band [{lo:.0f}, {hi:.0f}] - "
              "the published design is reproduced.")
    else:
        print(f"  ATTENTION: 61 lies OUTSIDE [{lo:.0f}, {hi:.0f}].")
        print("  The re-implementation may differ from the notebook (control-file")
        print("  row order changes each seed draw; see 07a header). The waterfall")
        print("  is then still valid as a 're-implemented legacy design', but the")
        print("  letter wording must say so - do NOT claim exact reproduction.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
