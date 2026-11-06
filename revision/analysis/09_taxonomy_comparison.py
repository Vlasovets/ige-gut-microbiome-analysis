#!/usr/bin/env python
"""
Module 09 - compare the original taxonomy against a documented RESCRIPt reference.

`taxonomy_clean.csv` came from a RESCRIPt classification whose reference
construction is recorded nowhere: no RESCRIPt call appears in either preprocessing
notebook. Rather than guess the original settings, module 00e builds a reference
with published best-practice parameters - SILVA 138.2 SSU NR99, culled,
length-filtered, dereplicated, and **extracted to the V3-V4 window bounded by
341F/785R** - and reclassifies the same ASVs.

This is a SEPARATE analysis, not a reproduction. Two things differ from the
original at once (reference construction and QIIME 2 version: 2026.7 vs 2021.11),
so disagreement cannot be attributed to either alone. It is reported as an
independent reference point, and `taxonomy_clean.csv` remains the record of what
the manuscript reports.

The question it answers concretely: only 12.4% of KORA ASVs currently carry a
species label. If region extraction lifts that substantially, the species-level
framing in the revision needs revisiting. If it does not, the sparsity is
inherent to a V3-V4 amplicon against SILVA, which is a stronger and more
defensible statement than the one currently in the draft.

Usage:
  python revision/analysis/09_taxonomy_comparison.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

REVISION = Path(__file__).resolve().parent.parent
CONFIG = REVISION / "config" / "config.yaml"
OUT = REVISION / "results" / "tables"
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PREFIX = {r: f"{r[0]}__" for r in RANKS}
PREFIX["domain"] = "d__"
UNKNOWN = ("unknown", "unclassif", "uncult", "nan", "")

# SILVA carries placeholder "species" labels that are not species at all -
# environmental or metagenome-assembly stand-ins. Counting them as assignments
# inflates the species-level rate, so they are tracked separately.
PLACEHOLDER_SPECIES = {
    "metagenome", "gut_metagenome", "human_gut", "unidentified",
    "uncultured_bacterium", "uncultured_organism", "bacterium",
    "unidentified_bacterium", "human_gut_metagenome", "uncultured",
}

# SILVA 138.2 adopts the ICNP-standardised phylum names and the Prevotella
# split. These are RENAMES, not disagreements; scoring them as mismatches makes
# the two taxonomies look far more discordant than they are.
SYNONYMS = {
    ("Firmicutes", "Bacillota"), ("Proteobacteria", "Pseudomonadota"),
    ("Actinobacteriota", "Actinomycetota"), ("Cyanobacteria", "Cyanobacteriota"),
    ("Desulfobacterota", "Thermodesulfobacteriota"),
    ("Euryarchaeota", "Methanobacteriota"),
    ("Campilobacterota", "Campylobacterota"), ("Chloroflexi", "Chloroflexota"),
    ("Halobacterota", "Halobacteriota"), ("Bacteroidota", "Bacteroidota"),
    # Genus-level reclassifications (Hitch et al.)
    ("Prevotella", "Segatella"), ("Prevotella", "Leyella"),
    ("Butyricicoccus", "Agathobaculum"),
}


def same_taxon(a, b) -> bool:
    """Equal, or a documented rename between the two reference versions."""
    if a == b:
        return True
    return (a, b) in SYNONYMS


def clean(v, pre):
    s = str(v).replace(pre, "").strip().rstrip(";").strip()
    return None if (not s or s.lower() in UNKNOWN
                    or any(u in s.lower() for u in UNKNOWN[:3])) else s


def parse_qiime_taxonomy(path: Path) -> pd.DataFrame:
    """Split a QIIME 2 'Taxon' column into fixed ranks."""
    t = pd.read_csv(path, sep="\t")
    idcol = t.columns[0]
    t = t.set_index(idcol)
    rows = {}
    for fid, tax in t["Taxon"].items():
        parts = [p.strip() for p in str(tax).split(";")]
        d = {}
        for r in RANKS:
            pre = PREFIX[r]
            hit = next((p for p in parts if p.startswith(pre)), None)
            d[r] = clean(hit, pre) if hit else None
        rows[fid] = d
    out = pd.DataFrame.from_dict(rows, orient="index")
    if "Confidence" in t.columns:
        out["confidence"] = t["Confidence"]
    return out


def main() -> int:
    cfg = yaml.safe_load(CONFIG.read_text())
    kora = Path(cfg["paths"]["kora_data"])
    new_path = kora / "q2_rescript_taxonomy" / "kora_taxonomy_rescript_v3v4.tsv"
    OUT.mkdir(parents=True, exist_ok=True)

    if not new_path.exists():
        print(f"RESCRIPt taxonomy not found: {new_path}")
        print("Run revision/slurm/00e_rescript_taxonomy.sh first.")
        return 1

    old = pd.read_csv(kora / "taxonomy_clean.csv", index_col=0)
    for r in RANKS:
        old[r] = old[r].map(lambda v, p=PREFIX[r]: clean(v, p))
    new = parse_qiime_taxonomy(new_path)

    shared = old.index.intersection(new.index)
    print("=== Taxonomy comparison: original vs documented RESCRIPt V3-V4 ===")
    print(f"  original ASVs : {len(old):,}")
    print(f"  RESCRIPt ASVs : {len(new):,}")
    print(f"  shared        : {len(shared):,}")
    o, n = old.loc[shared], new.loc[shared]

    # ── Assignment rate by rank ──────────────────────────────────────────
    print(f"\n  {'rank':9s} {'original':>10s} {'RESCRIPt':>10s} {'change':>9s}")
    rows = []
    for r in RANKS:
        a = o[r].notna().mean() * 100
        b = n[r].notna().mean() * 100
        print(f"  {r:9s} {a:9.1f}% {b:9.1f}% {b - a:+8.1f}pp")
        rows.append({"rank": r, "assigned_original_pct": a,
                     "assigned_rescript_pct": b, "change_pp": b - a})

    # ── Placeholder species labels ───────────────────────────────────────
    ph_new = n["species"].isin(PLACEHOLDER_SPECIES)
    ph_old = o["species"].isin(PLACEHOLDER_SPECIES)
    inf_new = (n["species"].notna() & ~ph_new).mean() * 100
    inf_old = (o["species"].notna() & ~ph_old).mean() * 100
    print(f"\n  Species labels that are placeholders, not species "
          f"(metagenome, human_gut, ...):")
    print(f"    original : {int(ph_old.sum()):,}   RESCRIPt: {int(ph_new.sum()):,}")
    print(f"  INFORMATIVE species assignment (placeholders excluded):")
    print(f"    original {inf_old:.1f}%  ->  RESCRIPt {inf_new:.1f}%  "
          f"({inf_new - inf_old:+.1f}pp)")

    # ── Agreement where both assign ──────────────────────────────────────
    # Reported two ways: raw string equality, and after collapsing documented
    # renames. The raw figure alone is misleading - SILVA 138.2 renamed most
    # phyla, so raw phylum "agreement" looks catastrophic when nothing moved.
    print(f"\n  {'rank':9s} {'both':>8s} {'raw agree':>11s} {'w/ renames':>12s}")
    for i, r in enumerate(RANKS):
        both = o[r].notna() & n[r].notna()
        raw = (o.loc[both, r] == n.loc[both, r]).sum()
        adj = sum(same_taxon(a, b)
                  for a, b in zip(o.loc[both, r], n.loc[both, r]))
        praw = 100 * raw / both.sum() if both.sum() else float("nan")
        padj = 100 * adj / both.sum() if both.sum() else float("nan")
        print(f"  {r:9s} {both.sum():8,} {praw:10.1f}% {padj:11.1f}%")
        rows[i]["n_both_assigned"] = int(both.sum())
        rows[i]["n_agree_raw"] = int(raw)
        rows[i]["agreement_raw_pct"] = praw
        rows[i]["n_agree_with_renames"] = int(adj)
        rows[i]["agreement_with_renames_pct"] = padj
    pd.DataFrame(rows).to_csv(OUT / "taxonomy_comparison_by_rank.csv", index=False)

    # ── The headline taxon ───────────────────────────────────────────────
    print("\n=== What happens to the headline species signal ===")
    pc_old = o.index[o["species"].astype(str).str.contains("copri", case=False,
                                                           na=False)]
    print(f"  ASVs labelled *copri* originally : {len(pc_old)}")
    if len(pc_old):
        relabel = n.loc[pc_old, "species"].fillna("<unassigned>").value_counts()
        print("  their RESCRIPt species labels:")
        for lab, cnt in relabel.head(8).items():
            print(f"    {cnt:5d}  {lab}")
        g = n.loc[pc_old, "genus"].fillna("<unassigned>").value_counts()
        print("  their RESCRIPt genus labels:")
        for lab, cnt in g.head(5).items():
            print(f"    {cnt:5d}  {lab}")

    pc_new = n.index[n["species"].astype(str).str.contains("copri", case=False,
                                                           na=False)]
    print(f"\n  ASVs labelled *copri* by RESCRIPt: {len(pc_new)}")
    print(f"  overlap with the original set    : {len(set(pc_old) & set(pc_new))}")

    # Is the genus itself stable, even if the species label moves?
    for gname in ("Prevotella", "Segatella"):
        oc = (o["genus"] == gname).sum()
        nc = (n["genus"] == gname).sum()
        if oc or nc:
            print(f"  genus {gname:11s} original {oc:5d} | RESCRIPt {nc:5d}")

    # ── Species-level detail, for the framing decision ───────────────────
    print("\n=== Most frequent species labels ===")
    cmp = pd.DataFrame({
        "original": o["species"].value_counts().head(8),
        "rescript": n["species"].value_counts().head(8),
    })
    print(cmp.to_string())

    both_sp = o["species"].notna() & n["species"].notna()
    if both_sp.sum():
        disagree = o.loc[both_sp & (o["species"] != n["species"])]
        print(f"\n  ASVs with a species label in BOTH: {both_sp.sum():,}")
        print(f"  of those, label differs           : {len(disagree):,} "
              f"({100 * len(disagree) / both_sp.sum():.1f}%)")

    full = o.join(n, lsuffix="_original", rsuffix="_rescript")
    full.to_csv(OUT / "taxonomy_comparison_per_asv.csv")
    print(f"\nWrote {OUT / 'taxonomy_comparison_by_rank.csv'}")
    print(f"Wrote {OUT / 'taxonomy_comparison_per_asv.csv'}")

    # ── Verdict for the species-framing decision ─────────────────────────
    sp_old = inf_old
    sp_new = inf_new
    print("\n=== Verdict ===")
    if sp_new - sp_old > 15:
        print(f"  Species assignment rises {sp_new - sp_old:.1f}pp "
              f"({sp_old:.1f}% -> {sp_new:.1f}%). Region extraction matters "
              "materially;\n  the species-level framing should be revisited and the "
              "reference choice reported.")
    elif sp_new - sp_old > 5:
        print(f"  Species assignment rises modestly ({sp_old:.1f}% -> {sp_new:.1f}%). "
              "Worth reporting,\n  but does not overturn the provisional framing.")
    else:
        print(f"  Species assignment barely moves ({sp_old:.1f}% -> {sp_new:.1f}%). "
              "The sparsity is\n  inherent to a V3-V4 amplicon against SILVA, not an "
              "artifact of reference choice.\n  That STRENGTHENS the provisional "
              "species-level framing - it now rests on a\n  tested alternative rather "
              "than an assumption.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
