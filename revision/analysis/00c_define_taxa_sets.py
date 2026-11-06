#!/usr/bin/env python
"""
Stage 0: freeze the network taxon set and build the genus -> phylum lookup.

Two prerequisites for module 03 that do not exist upstream:

1. "Top 50 most abundant genera cohort-wide" must be a reproducible artifact, not
   a phrase. The published network analysis used filtTax = "highestVar" with
   filtTaxPar = list(highestVar = 50), which (a) selects on VARIANCE, not
   abundance, and (b) was INERT - only 41-43 family-level taxa existed against a
   threshold of 50, so nothing was filtered. Worse, because the filter ran inside
   each netConstruct call, the control network was re-estimated on a DIFFERENT
   taxa subset in every comparison. Freezing one manifest here removes all three
   problems by construction.

2. 00b_aggregate_otu.py writes tax_table.csv only when level == "family"
   (see its line ~154), so genus-level networks have no phylum lookup for node
   colouring.

Ranking metric: mean of WITHIN-SAMPLE proportions, mean_j(x_ij / sum_i x_ij).
Not colSums/total, which is dominated by high-depth samples.

"Cohort-wide" = the 508 IgE-typed individuals, not the 2034 sequenced and not the
442 matched. Using the full typed cohort makes the taxa set independent of any
partition or of the matching, so every network is built on identical nodes.

Usage:
  python revision/analysis/00c_define_taxa_sets.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

REVISION = Path(__file__).resolve().parent.parent
CONFIG = REVISION / "config" / "config.yaml"
OUT = REVISION / "results" / "taxa_sets"

UNKNOWN = ("unknown", "unclassif", "uncult")
PREFIX = {"family": "f__", "genus": "g__", "species": "s__"}


def clean_label(raw, prefix: str):
    """Mirror 00b_aggregate_otu.py::_clean_label exactly, so labels match."""
    if pd.isna(raw):
        return None
    lab = str(raw).replace(prefix, "").strip().rstrip(";").strip()
    if not lab or lab in ("nan", "NA"):
        return None
    if any(tok in lab.lower() for tok in UNKNOWN):
        return None
    return lab


def main() -> int:
    cfg = yaml.safe_load(CONFIG.read_text())
    kora = Path(cfg["paths"]["kora_data"])
    net = cfg["networks"]
    OUT.mkdir(parents=True, exist_ok=True)

    print("=== Freezing the network taxon set ===")
    print(f"  level={net['level']}  top_n={net['top_n']}  "
          f"metric={net['abundance_metric']}")

    feat = pd.read_csv(kora / "feature_table.tsv", sep="\t", index_col=0)
    feat.columns = feat.columns.astype(str)
    tax = pd.read_csv(kora / "taxonomy_clean.csv", index_col=0)
    tax.index = tax.index.astype(str)
    print(f"  feature table: {feat.shape[0]} ASVs x {feat.shape[1]} samples")

    # Cohort = the IgE-typed individuals present in the feature table.
    ige = pd.read_csv(kora / "ige_kora_subset_all.csv", index_col="u3_16s_id",
                      low_memory=False)
    cohort = [s for s in ige.index.astype(str) if s in feat.columns]
    print(f"  cohort (IgE-typed, with sequencing): {len(cohort)}")
    feat = feat[cohort]

    level = net["level"]
    lab = {a: clean_label(tax.at[a, level] if level in tax.columns else None,
                          PREFIX[level]) for a in feat.index.astype(str)}
    feat = feat.copy()
    feat["_label"] = [lab.get(a) for a in feat.index.astype(str)]
    n_before = len(feat)
    feat = feat.dropna(subset=["_label"])
    print(f"  ASVs with a valid {level} label: {len(feat)} / {n_before}")

    agg = feat.groupby("_label").sum()          # genus x samples
    print(f"  aggregated {level}-level features: {agg.shape[0]}")

    # Prevalence on the cohort.
    prev = (agg > 0).mean(axis=1)
    thresh = {"prev10": 0.10, "prev05": 0.05, "prev01": 0.01}[net["prevalence"]]
    keep = prev[prev >= thresh].index
    agg = agg.loc[keep]
    print(f"  after prevalence >= {thresh:.0%}: {agg.shape[0]}")

    # Mean of within-sample proportions.
    depth = agg.sum(axis=0)
    prop = agg.div(depth.replace(0, np.nan), axis=1)
    ranked = pd.DataFrame({
        "mean_rel_abundance": prop.mean(axis=1),
        "sd_rel_abundance": prop.std(axis=1),
        "prevalence": prev.reindex(agg.index),
        "n_asvs_collapsed": feat.groupby("_label").size().reindex(agg.index),
    }).sort_values(["mean_rel_abundance", "prevalence"],
                   ascending=[False, False])
    # Deterministic tie-break so the manifest is byte-reproducible.
    ranked = ranked.sort_values(
        ["mean_rel_abundance"], ascending=False, kind="mergesort")
    ranked.insert(0, "rank", range(1, len(ranked) + 1))
    ranked.index.name = level

    ranked.to_csv(OUT / f"{level}_all_{net['prevalence']}_ranked.csv")

    sets = {}
    for n in sorted({net["top_n"], *net["robustness_top_n"]}):
        sel = ranked.head(n)
        sets[n] = list(sel.index)
        sel.to_csv(OUT / f"{level}_top{n}_cohort.csv")
        cum = sel.mean_rel_abundance.sum()
        print(f"  top {n:3d}: cumulative mean relative abundance = {cum:.3f}")

    # Nesting must hold, or the robustness comparison is meaningless.
    ns = sorted(sets)
    for a, b in zip(ns, ns[1:]):
        ok = set(sets[a]).issubset(sets[b])
        print(f"  QC top{a} subset of top{b}: {'PASS' if ok else 'FAIL'}")
        if not ok:
            return 1

    top = sets[net["top_n"]]
    print(f"\n  Top 10 by abundance: {', '.join(top[:10])}")

    # ── genus -> phylum lookup, via one representative ASV per label ──────
    rep = {}
    for a, l in lab.items():
        if l is not None and l not in rep:
            rep[l] = a
    rows = []
    for g in ranked.index:
        a = rep.get(g)
        phy = ord_ = ""
        if a is not None and a in tax.index:
            phy = str(tax.at[a, "phylum"]).replace("p__", "").strip().rstrip(";") \
                if "phylum" in tax.columns else ""
            ord_ = str(tax.at[a, "order"]).replace("o__", "").strip().rstrip(";") \
                if "order" in tax.columns else ""
        rows.append({level.capitalize(): g, "Phylum": phy, "Order": ord_,
                     "representative_asv": a})
    tt = pd.DataFrame(rows).set_index(level.capitalize())
    tt.to_csv(OUT / f"tax_table_{level}.csv")
    n_unk = int((tt.loc[top, "Phylum"].isin(["", "nan"])).sum())
    print(f"  genus -> phylum table: {len(tt)} rows "
          f"({n_unk} of the top {len(top)} lack a phylum)")
    print(f"  phyla among the top {len(top)}: "
          f"{sorted(set(tt.loc[top, 'Phylum']) - {'', 'nan'})}")

    meta = {
        "level": level, "prevalence": net["prevalence"],
        "abundance_metric": net["abundance_metric"],
        "cohort_n": len(cohort),
        "n_features_after_prevalence": int(agg.shape[0]),
        "sets": {str(k): v for k, v in sets.items()},
    }
    (OUT / "manifest.json").write_text(json.dumps(meta, indent=2))
    print(f"\nWrote manifest and taxa sets to {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
