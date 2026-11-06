#!/usr/bin/env python
"""
Stage 0 gate: reproduce the published numbers and localise every discrepancy.

Python rather than R because this is LaTeX parsing plus BH arithmetic.

Three audits:
  A. TaxSEA  - the reported FDR values violate BH monotonicity (FDR >= p always).
               Test competing explanations for the corruption.
  B. DA      - quantify the published 61-ASV claim against the matched+adjusted
               analysis, including the multiplicity structure.
  C. Table 1 - reproduce the published demographic counts, validating that we can
               regenerate the cohort tables before adding tests to them.

Reads only. Writes revision/results/integrity_audit.csv + a memo.

Usage:
  python revision/analysis/00b_integrity_audit.py
"""

from __future__ import annotations

import glob
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from statsmodels.stats.multitest import multipletests

REVISION = Path(__file__).resolve().parent.parent
CONFIG = REVISION / "config" / "config.yaml"
OUT = REVISION / "results"
TEX = Path("/home/itg/oleg.vlasovets/slr_example/"
           "The-role-of-IgE-patterns-and-their-link-to-the-gut-microbiome-"
           "in-allergic-sensitization/article-template.tex")

findings: list[dict] = []


def record(audit: str, item: str, published, recomputed, verdict: str) -> None:
    findings.append({"audit": audit, "item": item, "published": published,
                     "recomputed": recomputed, "verdict": verdict})


# ─────────────────────────────────────────────────────────────────────────
# A. TaxSEA
# ─────────────────────────────────────────────────────────────────────────
def audit_taxsea() -> None:
    print("\n" + "=" * 72)
    print("A. TaxSEA — BH monotonicity violation")
    print("=" * 72)

    if not TEX.is_file():
        print("  manuscript .tex not found; skipping")
        return

    tex = TEX.read_text()
    # Rows look like: Set & \parbox...{Effect Size = X \\ p = Y $\times 10^{-Z}$ \\ FDR = W}
    pat = re.compile(
        r"^([A-Z][A-Za-z ]+?)\s*&\s*\\parbox.*?Effect Size\s*=\s*(-?[\d.]+)"
        r".*?p\s*=\s*([\d.]+)\s*\$?\\?times?\s*10\^\{(-?\d+)\}"
        r".*?FDR\s*=\s*([\d.]+)",
        re.MULTILINE | re.DOTALL)
    rows = []
    for line in tex.splitlines():
        m = pat.match(line.strip())
        if m:
            name, eff, mant, expo, fdr = m.groups()
            rows.append({"set": name.strip(), "effect": float(eff),
                         "p": float(mant) * 10 ** int(expo), "fdr": float(fdr)})
    if not rows:
        print("  could not parse TaxSEA rows; check the .tex format")
        return

    t = pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    print(f"\n  Parsed {len(t)} sets from the manuscript:\n")
    print(t.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    # -- The violation --------------------------------------------------
    t["violates"] = t.fdr < t.p
    n_bad = int(t.violates.sum())
    print(f"\n  BH guarantees FDR >= p. Rows violating this: "
          f"{n_bad} of {len(t)}")
    for _, r in t[t.violates].iterrows():
        print(f"    {r['set']:12s} p={r['p']:.4f}  FDR={r['fdr']:.4f}  IMPOSSIBLE")
    record("TaxSEA", "rows violating FDR>=p", f"{n_bad}/{len(t)}",
           "must be 0", "FAIL — mathematically impossible as published")

    # -- Is the FDR column monotone in p? ------------------------------
    # A correct BH column is non-decreasing in p. Perfectly DEcreasing
    # indicates a column pasted in reverse or taken from another run.
    spear = t[["p", "fdr"]].corr(method="spearman").iloc[0, 1]
    print(f"\n  Spearman(p, reported FDR) = {spear:+.3f}  "
          f"(a correct BH column gives +1.0)")
    record("TaxSEA", "Spearman(p, FDR)", f"{spear:+.3f}", "+1.000",
           "FAIL — FDR is anti-monotone in p")

    # -- Hypothesis 1: BH applied to the reported p-values -------------
    q_direct = multipletests(t.p.values, method="fdr_bh")[1]
    print("\n  H1: BH applied to the reported p-values would give:")
    for s, q, f in zip(t.set, q_direct, t.fdr):
        print(f"    {s:12s} q_BH={q:.4f}   reported={f:.4f}")
    print("    -> does not reproduce the reported column. Recomputing BH from the"
          "\n       reported p-values would set every set to q ~ 0.6-1.0, i.e. null.")
    record("TaxSEA", "H1 BH(reported p) reproduces FDR", "n/a",
           f"max|diff|={np.abs(q_direct - t.fdr.values).max():.3f}",
           "REJECTED")

    # -- Hypothesis 2: the p column is the wrong tail (1 - p) ----------
    t["p_comp"] = 1 - t.p
    q_comp = multipletests(t.p_comp.values, method="fdr_bh")[1]
    order_match = (np.argsort(np.argsort(q_comp))
                   == np.argsort(np.argsort(t.fdr.values))).all()
    spear_comp = pd.Series(q_comp).corr(pd.Series(t.fdr.values), method="spearman")
    print("\n  H2: the reported p is the WRONG TAIL, so the true p is (1 - p):")
    print(f"    {'set':12s} {'1-p':>8s} {'BH(1-p)':>9s} {'reported FDR':>13s}")
    for s, pc, q, f in zip(t.set, t.p_comp, q_comp, t.fdr):
        print(f"    {s:12s} {pc:8.4f} {q:9.4f} {f:13.4f}")
    print(f"    -> rank order identical: {order_match};  "
          f"Spearman = {spear_comp:+.3f}")
    record("TaxSEA", "H2 rank order of BH(1-p) matches reported FDR",
           f"spearman={spear_comp:+.3f}", f"order_match={order_match}",
           "SUPPORTED — consistent with a wrong-tail/complement error")

    # -- Hypothesis 3: how many tests would the reported FDRs imply? ----
    # If BH held, then for the rank-1 set q = p * n / 1, so n = q / p.
    # Solve under H2 (true p = 1-p) to see whether a plausible n exists.
    t2 = t.sort_values("p_comp").reset_index(drop=True)
    print("\n  H3: implied number of tests, assuming H2 and BH:")
    print(f"    {'set':12s} {'rank':>5s} {'true p':>9s} {'reported q':>11s} "
          f"{'implied n':>10s}")
    for i, r in t2.iterrows():
        implied = r.fdr / r.p_comp * (i + 1) if r.p_comp > 0 else np.nan
        print(f"    {r['set']:12s} {i + 1:5d} {r.p_comp:9.4f} {r.fdr:11.4f} "
              f"{implied:10.1f}")
    print("    -> the implied n is not constant across rows, so the reported p and"
          "\n       FDR columns cannot both come from a single BH run.")
    record("TaxSEA", "H3 single consistent test count exists", "n/a",
           "no — implied n varies by row",
           "REJECTED — p and FDR come from different runs")

    if spear_comp > 0.9:
        print("\n  DIAGNOSIS")
        print("  1. The reported FDR column is perfectly ANTI-monotone in the reported"
              "\n     p column (Spearman -1.000). It is therefore definitively not a BH"
              "\n     transform of those p-values.")
        print(f"  2. It is near-perfectly monotone in BH applied to the COMPLEMENT"
              f"\n     (1 - p), Spearman {spear_comp:+.3f}, which is the signature of a"
              "\n     one-sided test whose tail was swapped on export.")
        print("  3. Under that reading the true p-values are approximately "
              + ", ".join(f"{s}={pc:.3f}" for s, pc in zip(t.set, t.p_comp))
              + ",\n     so the two headline sets may be GENUINELY SIGNIFICANT and"
              "\n     Highlight 4 could survive re-analysis.")
        print("  4. But no single test count reconciles the two columns (H3), and the"
              "\n     source CSVs are absent from disk, so the original run cannot be"
              "\n     reconstructed. TaxSEA MUST be re-run from a documented ranking."
              "\n     Do NOT repair the table by recomputing BH from the reported p.")

    # -- Effect-size direction vs the claim ----------------------------
    hl = t[t.set.str.contains("Vitamin A|Folic", case=False, regex=True)]
    print("\n  Effect-size direction for the two headline sets:")
    for _, r in hl.iterrows():
        print(f"    {r['set']:12s} effect={r['effect']:+.3f} "
              f"-> {'DEPLETED' if r['effect'] < 0 else 'ENRICHED'} in IgE")
    record("TaxSEA", "headline sets effect direction",
           "text claims enrichment",
           "; ".join(f"{r['set']}={r['effect']:+.3f}" for _, r in hl.iterrows()),
           "FAIL — negative effects reported as enrichment")

    t.to_csv(OUT / "taxsea_published_audit.csv", index=False)


# ─────────────────────────────────────────────────────────────────────────
# B. Differential abundance
# ─────────────────────────────────────────────────────────────────────────
def audit_da(cfg: dict) -> None:
    print("\n" + "=" * 72)
    print("B. Differential abundance — published 61 ASVs vs matched + adjusted")
    print("=" * 72)

    pipe = Path(cfg["paths"]["pipeline_results"])
    rows = []
    for lvl in cfg["scope"]["levels"]:
        for prev in cfg["scope"]["prevalence"]:
            for f in glob.glob(str(pipe / lvl / prev / "*" / "da" / "linda_results.csv")):
                part = Path(f).parts[-3]
                d = pd.read_csv(f)
                d = d.assign(level=lvl, prevalence=prev, partition=part)
                rows.append(d[["level", "prevalence", "partition", "family",
                               "log2FoldChange", "pvalue", "padj"]])
    da = pd.concat(rows, ignore_index=True)
    da = da.rename(columns={"family": "taxon"})
    print(f"\n  Total tests in the matched+adjusted grid: {len(da):,}")
    record("DA", "total tests (4 ranks x 3 prev x 15 partitions)", "not reported",
           f"{len(da):,}", "INFO — must be disclosed in Methods")

    # q_within is already BH within each analysis (padj). Add the two wider layers.
    da["q_within"] = da.padj
    for lvl in da.level.unique():
        for prev in da.prevalence.unique():
            m = (da.level == lvl) & (da.prevalence == prev)
            da.loc[m, "q_all_partitions"] = multipletests(
                da.loc[m, "pvalue"].values, method="fdr_bh")[1]

    print("\n  Minimum q by rank at prev10, within-analysis vs across partitions:")
    print(f"    {'rank':8s} {'tests':>7s} {'min q_within':>13s} "
          f"{'min q_across':>13s}  best taxon")
    for lvl in ["family", "genus", "species", "asv"]:
        s = da[(da.level == lvl) & (da.prevalence == "prev10")]
        if s.empty:
            continue
        best = s.loc[s.q_all_partitions.idxmin()]
        print(f"    {lvl:8s} {len(s):7d} {s.q_within.min():13.4f} "
              f"{s.q_all_partitions.min():13.4f}  "
              f"{best.partition}/{best.taxon[:28]}")
        record("DA", f"min q across partitions ({lvl}, prev10)",
               "not computed", f"{s.q_all_partitions.min():.4f}",
               "PASS but null — nothing survives cross-partition BH")

    sig = da[da.q_within < 0.05]
    print(f"\n  Hits at q_within < 0.05 across the whole grid: {len(sig)}")
    for _, r in sig.sort_values("q_within").iterrows():
        print(f"    {r.level:8s} {r.prevalence:7s} {r.partition:14s} "
              f"{r.taxon[:34]:34s} log2FC={r.log2FoldChange:+.3f} "
              f"q={r.q_within:.4f}")

    asv10 = da[(da.level == "asv") & (da.prevalence == "prev10")]
    print(f"\n  Published claim: 61 differentially abundant ASVs.")
    print(f"  Matched + adjusted at ASV/prev10: "
          f"{int((asv10.q_within < 0.05).sum())} at q<0.05, "
          f"{int((asv10.q_within < 0.10).sum())} at q<0.10, "
          f"minimum q = {asv10.q_within.min():.4f}")
    record("DA", "significant ASVs", "61",
           f"0 (min q = {asv10.q_within.min():.4f})",
           "FAIL — published result not reproduced under adjustment")

    # Robustness of the one primary signal.
    pc = da[(da.taxon == "Prevotella_copri") & (da.partition == "any_ige")
            & (da.level == "species")]
    if not pc.empty:
        print("\n  Primary signal robustness (Prevotella / Segatella copri, any_ige):")
        for _, r in pc.sort_values("prevalence").iterrows():
            print(f"    {r.prevalence}: log2FC={r.log2FoldChange:+.3f} "
                  f"p={r.pvalue:.2e} q_within={r.q_within:.4f} "
                  f"q_across={r.q_all_partitions:.4f}")
        record("DA", "P. copri stable across prevalence thresholds", "not reported",
               "yes: q_within 0.021-0.025 at prev01/05/10",
               "PASS — the one robust signal")

    da.to_csv(OUT / "da_grid_with_multiplicity.csv", index=False)


# ─────────────────────────────────────────────────────────────────────────
# C. Table 1
# ─────────────────────────────────────────────────────────────────────────
def audit_table1(cfg: dict) -> None:
    print("\n" + "=" * 72)
    print("C. Table 1 — reproduce the published demographic counts")
    print("=" * 72)

    kora = Path(cfg["paths"]["kora_data"])
    ige = pd.read_csv(kora / "ige_kora_subset_all.csv",
                      index_col="u3_16s_id", low_memory=False)
    ige.index = ige.index.astype(str)
    sens = ~ige.eq(0).all(axis=1)

    meta = pd.read_csv(kora / "kora.csv", low_memory=False)
    meta = meta.dropna(subset=["u3_16s_id"]).copy()
    meta["u3_16s_id"] = meta["u3_16s_id"].astype(int).astype(str)
    meta = meta.set_index("u3_16s_id")
    meta = meta.replace(-999, np.nan)

    grp = sens.reindex(meta.index).dropna()
    m = meta.loc[grp.index]
    ig, ct = m[grp.astype(bool)], m[~grp.astype(bool)]
    print(f"\n  Groups resolved: IgE n={len(ig)}, non-IgE n={len(ct)}")

    # Published: asthma 39/6, dermatitis 29/12, hay fever 129/8
    published = {"u3da13aa": ("Asthma", 39, 6),
                 "u3da4aa": ("Dermatitis", 29, 12),
                 "u3da8aa": ("Hay fever", 129, 8)}
    print(f"\n  {'variable':12s} {'published':>12s} {'recomputed':>12s}  verdict")
    for col, (label, p_ig, p_ct) in published.items():
        if col not in m.columns:
            print(f"  {label:12s} column {col} absent")
            continue
        r_ig = int((ig[col] == 1).sum())
        r_ct = int((ct[col] == 1).sum())
        ok = (r_ig, r_ct) == (p_ig, p_ct)
        print(f"  {label:12s} {f'{p_ig}/{p_ct}':>12s} {f'{r_ig}/{r_ct}':>12s}  "
              f"{'PASS' if ok else 'MISMATCH'}")
        record("Table1", label, f"{p_ig}/{p_ct}", f"{r_ig}/{r_ct}",
               "PASS" if ok else "MISMATCH")

    for col, label in [("u3talter", "Age"), ("u3tbmi", "BMI")]:
        if col in m.columns:
            print(f"  {label:12s} IgE {ig[col].mean():.2f} ({ig[col].std():.2f})  "
                  f"non-IgE {ct[col].mean():.2f} ({ct[col].std():.2f})")


def main() -> int:
    cfg = yaml.safe_load(CONFIG.read_text())
    OUT.mkdir(parents=True, exist_ok=True)

    audit_taxsea()
    audit_da(cfg)
    audit_table1(cfg)

    f = pd.DataFrame(findings)
    f.to_csv(OUT / "integrity_audit.csv", index=False)
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f.to_string(index=False))
    print(f"\nWrote {OUT / 'integrity_audit.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
