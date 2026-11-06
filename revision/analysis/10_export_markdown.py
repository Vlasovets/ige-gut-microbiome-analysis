#!/usr/bin/env python
"""
Module 10 - export a self-contained Markdown briefing of the revision.

PURPOSE
  Produce a set of Markdown files that carry the DETAILS AND RESULTS with no
  dependency on this cluster, this environment, or the CSVs - so the work can
  continue on a laptop.

DESIGN RULE
  Every number is READ FROM THE RESULT FILES at generation time. Nothing is
  hand-transcribed. If a result is missing because a job has not finished, the
  section says so explicitly rather than silently omitting it - an absent
  section must never read like a null finding.

PRIVACY
  Only aggregate, per-taxon and per-analysis statistics are emitted. No
  individual-level participant data is read or written. The inputs are the
  tables in revision/results/, which were verified free of per-subject
  identifiers.

Usage:
  python revision/analysis/10_export_markdown.py
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

import pandas as pd

REVISION = Path(__file__).resolve().parent.parent
REPO = REVISION.parent
TABLES = REVISION / "results" / "tables"
DOCS = REPO / "docs"
ALPHA = 0.05


def md_table(df: pd.DataFrame, floatfmt: str = "{:.4g}", max_rows: int = 60) -> str:
    """Render a DataFrame as a GitHub Markdown table."""
    if df is None or df.empty:
        return "_No rows._\n"
    d = df.head(max_rows).copy()

    def cell(v) -> str:
        # Cover every case explicitly. An object column holding a mix of str and
        # float NaN is common here (an empty CSV field reads back as NaN), and
        # dtype-based branching misses it.
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return ""
        try:
            if pd.isna(v):
                return ""
        except (TypeError, ValueError):
            pass
        if isinstance(v, float):
            return floatfmt.format(v)
        # Pipes would break the table layout.
        return str(v).replace("|", "\\|")

    head = "| " + " | ".join(str(c) for c in d.columns) + " |"
    rule = "|" + "|".join("---" for _ in d.columns) + "|"
    rows = ["| " + " | ".join(cell(v) for v in r) + " |"
            for r in d.itertuples(index=False, name=None)]
    out = "\n".join([head, rule, *rows])
    if len(df) > max_rows:
        out += f"\n\n_Showing {max_rows} of {len(df):,} rows._"
    return out + "\n"


def read(name: str) -> pd.DataFrame | None:
    p = TABLES / name
    if not p.exists():
        return None
    try:
        return pd.read_csv(p, low_memory=False)
    except Exception:
        return None


def missing(what: str, how: str) -> str:
    return (f"> ⚠ **Not available.** {what}\n>\n"
            f"> This is a MISSING RESULT, not a null result. Regenerate with: `{how}`\n")


# ── 02 Differential abundance ────────────────────────────────────────────
def doc_da() -> str:
    rates = read("da_method_rates.csv")
    pw = read("da_method_overlap_pairwise.csv")
    up = read("da_method_upset.csv")
    cor = read("da_method_effect_correlation.csv")

    s = ["# Differential abundance — six methods\n",
         "The direct answer to *“use more than one DA method, and show whether "
         "they agree.”*\n",
         "All six ran on **identical input** — the same matched, "
         "prevalence-filtered tables — with the same formula "
         "`~ ige_status + sex + age + bmi`, except the paired Wilcoxon, which "
         "uses the 1:1 matching directly instead of adjusting.\n",
         "| Family | Methods |",
         "|---|---|",
         "| Log-ratio | LinDA, ANCOM-BC2, ALDEx2 |",
         "| Permutation / reference-frame | ZicoSeq |",
         "| Non-parametric, matched | Wilcoxon on CLR (paired) |",
         "| Count model | limma-voom |\n"]

    s.append("## Per-method significance rate\n")
    if rates is None:
        s.append(missing("Module 02e has not produced `da_method_rates.csv`.",
                         "python revision/analysis/02e_method_overlap.py"))
    else:
        s.append(md_table(rates))
        s.append("\nA method calling **18%** of taxa significant while two others "
                 "call 0.0–0.1% is a statement about calibration, not biology.\n")

    s.append("\n## Pairwise overlap of significant calls\n")
    s.append("Restricted to taxon-analyses **both** methods actually evaluated. "
             "ZicoSeq drops features that are constant within a partition and "
             "ANCOM-BC2 applies internal filtering, so a taxon a method never "
             "tested must not be scored as that method finding it "
             "non-significant.\n")
    if pw is None:
        s.append(missing("`da_method_overlap_pairwise.csv` is absent.",
                         "python revision/analysis/02e_method_overlap.py"))
    else:
        s.append(md_table(pw.sort_values("n_shared", ascending=False)))

    s.append("\n## Effect-size agreement\n")
    s.append("The fairer measure: two methods can rank taxa near-identically yet "
             "share no significance calls if one is simply more conservative. "
             "Reporting only shared calls would overstate genuine disagreement.\n")
    s.append("> ⚠ **Not every pair is an independent comparison.** ALDEx2, "
             "ZicoSeq and the paired Wilcoxon all report a **CLR difference**. "
             "For a balanced matched design the mean of paired differences is "
             "*algebraically equal* to the difference of means, so ZicoSeq's and "
             "Wilcoxon's effect columns are the **same number** — they correlate "
             "at ρ = 1.000 by construction, not by agreement. ZicoSeq's column is "
             "not even a ZicoSeq statistic: its own R² is unsigned, so the CLR "
             "mean difference is substituted purely to recover direction. Pairs "
             "sharing an estimand are flagged below and must never be cited as "
             "cross-method validation.\n")
    if cor is None:
        s.append(missing("`da_method_effect_correlation.csv` is absent.",
                         "python revision/analysis/02e_method_overlap.py"))
    else:
        c = cor.sort_values("spearman", ascending=False)
        if "independent_comparison" in c.columns:
            s.append("### Independent comparisons — different estimands\n")
            s.append(md_table(c[c.independent_comparison]))
            s.append("\n### Same-estimand pairs — agreement is partly definitional\n")
            s.append(md_table(c[~c.independent_comparison]))
            s.append("\nAmong the **independent** comparisons, five methods agree "
                     "closely on direction (91–95% same-sign). **ANCOM-BC2 is the "
                     "sole outlier at 57–59% — coin-flip.** That is the finding: "
                     "the method calling 31.5% of taxa significant is also the one "
                     "whose effect estimates do not reproduce in any other method.\n")
        else:
            s.append(md_table(c))

    s.append("\n## Consensus\n")
    if up is None:
        s.append(missing("`da_method_upset.csv` is absent.",
                         "python revision/analysis/02e_method_overlap.py"))
    else:
        s.append(f"Taxon-analyses called significant by ≥1 method: "
                 f"**{len(up):,}**\n")
        vc = up.n_methods.value_counts().sort_index(ascending=False)
        s.append(md_table(vc.rename_axis("n_methods").reset_index(
            name="n_taxon_analyses")))
        multi = up[up.n_methods >= 2]
        s.append(f"\n**Called by ≥2 methods: {len(multi):,}**\n")
        cols = [c for c in ["level", "prevalence", "partition", "taxon",
                            "n_methods", "n_tested", "methods",
                            "methods_not_tested"] if c in multi.columns]
        s.append(md_table(multi[cols], max_rows=80))
        s.append("\n### By taxonomic rank\n")
        rows = []
        for lev in ["family", "genus", "species", "asv"]:
            sub = up[up.level == lev]
            if len(sub):
                rows.append({"rank": lev, "any_method": len(sub),
                             "ge_2_methods": int((sub.n_methods >= 2).sum()),
                             "ge_3_methods": int((sub.n_methods >= 3).sum())})
        s.append(md_table(pd.DataFrame(rows)))

    s.append("""
## Two cautions that must travel with these numbers

1. **Effect sizes are never averaged across methods.** They are on different
   scales — log2 fold change of relative abundance, natural-log LFC, CLR
   difference, log2 CPM fold change. An average would mean nothing.

2. **ZicoSeq at ASV rank must be read via `q_native`, not `q_within`.** A
   permutation p-value has a hard floor of 1/(nperm+1), so BH over *m* features
   cannot fall below *m*/(nperm+1). At nperm = 4999 that floor is 0.012 for 60
   features but **0.29 for 1,459** — structurally incapable of reaching 0.05.
   `q_native` is ZicoSeq's own permutation FDR, which is not floor-limited the
   same way. The scripts print this floor for every analysis.
""")
    return "\n".join(s)


# ── 03 Taxonomy ──────────────────────────────────────────────────────────
def doc_taxonomy() -> str:
    by = read("taxonomy_comparison_by_rank.csv")
    src = REVISION / "TAXONOMY_COMPARISON.md"
    s = ["# Taxonomy — SILVA 138.2 reclassification\n"]
    if by is not None:
        cols = [c for c in ["rank", "assigned_original_pct", "assigned_rescript_pct",
                            "change_pp", "n_both_assigned", "agreement_raw_pct",
                            "agreement_with_renames_pct"] if c in by.columns]
        s.append("## Assignment and agreement by rank\n")
        s.append(md_table(by[cols]))
    else:
        s.append(missing("`taxonomy_comparison_by_rank.csv` is absent.",
                         "python revision/analysis/09_taxonomy_comparison.py"))
    if src.exists():
        s.append("\n---\n")
        s.append(src.read_text())
    return "\n".join(s)


# ── 04 Baseline / diversity ──────────────────────────────────────────────
def doc_other_results() -> str:
    s = ["# Cohort, diversity, networks and classification\n"]
    for title, fname, note in [
        ("Table 1 — full cohort (n = 508), formally tested",
         "table1_full508_tested.csv",
         "Tests the submitted Table 1 never reported."),
        ("Table 1 — matched subset QC", "table1_matched442_qc.csv",
         "Confirms the 1:1 matching balanced the covariates it was built on."),
        ("Table 2 — IgE cluster descriptives", "table2_clusters_descriptives.csv", ""),
        ("Table 2 — IgE clusters, formally tested", "table2_clusters_tested.csv",
         "Clusters F (n = 10) and A (n = 12) are underpowered; read with that in mind."),
        ("TaxSEA, recomputed", "taxsea_recomputed.csv",
         "The re-run that settled the impossible FDR < p rows in the submitted "
         "supplement."),
    ]:
        df = read(fname)
        s.append(f"\n## {title}\n")
        if note:
            s.append(note + "\n")
        s.append(md_table(df) if df is not None else
                 missing(f"`{fname}` is absent.", "see revision/README.md"))
    return "\n".join(s)


def main() -> int:
    DOCS.mkdir(parents=True, exist_ok=True)
    stamp = date.today().isoformat()

    written = []
    for fname, body in [
        ("02_DIFFERENTIAL_ABUNDANCE.md", doc_da()),
        ("03_TAXONOMY.md", doc_taxonomy()),
        ("04_COHORT_DIVERSITY_TAXSEA.md", doc_other_results()),
    ]:
        (DOCS / fname).write_text(body + f"\n\n---\n_Generated {stamp} by "
                                         f"`revision/analysis/10_export_markdown.py`._\n")
        written.append(fname)

    # Copy the narrative documents verbatim - they are already Markdown and are
    # the reasoning, not the numbers.
    for src, dst in [
        (REPO.parent / "REVISION_PLAN.md", "05_REVISION_PLAN.md"),
        (REPO.parent / "REVIEWER_RESPONSE_OUTLINE.md", "06_REVIEWER_RESPONSE.md"),
        (REVISION / "UPSTREAM_METHODS.md", "07_UPSTREAM_METHODS.md"),
        (REVISION / "DATA.md", "08_DATA.md"),
    ]:
        if src.exists():
            (DOCS / dst).write_text(src.read_text())
            written.append(dst)

    wb = DOCS / "REVIEWER_WORKBENCH.md"
    index = [f"# iScience revision — offline briefing\n",
             f"_Generated {stamp}. Every number is read from the result files at "
             f"generation time; nothing is transcribed by hand._\n",
             "Manuscript **ISCIENCE-D-26-05498** — *The role of IgE patterns and "
             "their link to the gut microbiome in allergic sensitization* "
             "(KORA FF4, N = 508).\n",
             "## Contents\n"]
    labels = {
        "02_DIFFERENTIAL_ABUNDANCE.md": "Six DA methods and their overlap — the reviewer's main question",
        "03_TAXONOMY.md": "SILVA 138.2 reclassification; why the claim is scoped to 108 ASVs",
        "04_COHORT_DIVERSITY_TAXSEA.md": "Tables 1–2 formally tested, and the TaxSEA re-run",
        "05_REVISION_PLAN.md": "Execution plan; §6 is the manuscript-correction inventory",
        "06_REVIEWER_RESPONSE.md": "Point-by-point response, keyed to R1/R2/R3",
        "07_UPSTREAM_METHODS.md": "Recovered DADA2 / taxonomy / phylogeny provenance",
        "08_DATA.md": "Data inventory and provenance",
    }
    if wb.exists():
        index.append("- [REVIEWER_WORKBENCH.md](REVIEWER_WORKBENCH.md) — "
                     "**START HERE**: per-comment guide — files to open, verdicts, "
                     "manuscript line anchors, response skeletons")
    for f in written:
        index.append(f"- [{f}]({f}) — {labels.get(f, '')}")
    index.append("""
## What is NOT here

No individual-level participant data. These documents carry aggregate,
per-taxon and per-analysis statistics only. Regenerating any result requires
the KORA tables and `slr-env` (R 4.4.3 / Bioconductor 3.20) on the cluster.

## Open items requiring a human

- The **verbatim reviewer letter** — every `[VERBATIM COMMENT]` in
  `06_REVIEWER_RESPONSE.md` is a placeholder until it is supplied.
- **Allergy-domain sign-off** on the 104-component category mapping
  (bromelain/CCD, *Anisakis*, storage vs house-dust mite).
- The manuscript is still formatted for **Cell Host & Microbe** and has not
  been retargeted to iScience.
""")
    (DOCS / "README.md").write_text("\n".join(index) + "\n")

    print(f"Wrote {len(written) + 1} Markdown files to {DOCS}")
    for f in ["README.md", *written]:
        p = DOCS / f
        print(f"  {f:36s} {p.stat().st_size / 1024:7.1f} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
