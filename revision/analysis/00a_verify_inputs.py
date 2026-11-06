#!/usr/bin/env python
"""
Stage 0: verify and hash every input the revision depends on.

Purpose: make every downstream number traceable, and fail loudly now rather than
silently later if the read-only input archive has changed. Writes a manifest of
checksums plus the assertions that define the analysis scope.

Reads only. Writes only into revision/results/.

Usage:
  python revision/analysis/00a_verify_inputs.py
  python revision/analysis/00a_verify_inputs.py --quick   # skip the 124 MB hash
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd
import yaml

REVISION = Path(__file__).resolve().parent.parent
CONFIG = REVISION / "config" / "config.yaml"
OUT = REVISION / "results"

# Expected values, established by audit. A mismatch means an input changed.
EXPECTED_PAIRS = {
    "cluster_A": 12, "cluster_B": 23, "cluster_C": 25, "cluster_D": 45,
    "cluster_E": 22, "cluster_F": 10, "cluster_G": 21, "cluster_H": 40,
    "cluster_I": 70, "any_ige": 221, "single": 44, "multiple": 196,
    "high_score": 57, "medium_score": 107, "low_score": 104,
}
EXPECTED_FEATURES_PREV10 = {"family": 45, "genus": 112, "species": 44, "asv": 288}
EXPECTED_IGE_SHAPE = (508, 104)          # NOT 112 components
EXPECTED_IGE_SPLIT = (275, 233)          # under "any non-zero"

failures: list[str] = []
notes: list[str] = []


def check(cond: bool, msg: str) -> None:
    if cond:
        print(f"  PASS  {msg}")
    else:
        print(f"  FAIL  {msg}")
        failures.append(msg)


def sha256(path: Path, quick: bool = False) -> str:
    if quick and path.stat().st_size > 50_000_000:
        return "skipped-quick"
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="skip hashing files larger than 50 MB")
    args = ap.parse_args()

    cfg = yaml.safe_load(CONFIG.read_text())
    kora = Path(cfg["paths"]["kora_data"])
    pipe = Path(cfg["paths"]["pipeline_results"])
    OUT.mkdir(parents=True, exist_ok=True)

    print(f"=== Stage 0: input verification ===")
    print(f"  kora_data        {kora}")
    print(f"  pipeline_results {pipe}")

    # ---- 1. Raw inputs exist, and hash them ------------------------------
    print("\n[1] Raw inputs")
    raw = ["feature_table.tsv", "taxonomy_clean.csv", "kora.csv",
           "ige_kora_subset_all.csv", "cluster_df_latent.csv"]
    manifest = []
    for name in raw:
        p = kora / name
        check(p.is_file(), f"{name} present")
        if p.is_file():
            manifest.append({
                "path": str(p), "sha256": sha256(p, args.quick),
                "bytes": p.stat().st_size,
            })

    # ---- 2. Vendored utils source ---------------------------------------
    print("\n[2] utils source (CLR + pair matching implementations)")
    us = Path(cfg["paths"]["utils_source"])
    for mod in ["helper.py", "pair_matching.py", "preprocessing.py"]:
        check((us / mod).is_file(), f"utils/{mod} available to vendor")

    # ---- 3. IgE panel shape and group split -----------------------------
    print("\n[3] IgE panel")
    ige = pd.read_csv(kora / "ige_kora_subset_all.csv",
                      index_col="u3_16s_id", low_memory=False)
    check(ige.shape == EXPECTED_IGE_SHAPE,
          f"panel is {ige.shape} (expected {EXPECTED_IGE_SHAPE}; "
          f"manuscript says 112 components -- that is the platform, not the matrix)")

    sens = ~ige.eq(0).all(axis=1)
    split = (int(sens.sum()), int((~sens).sum()))
    check(split == EXPECTED_IGE_SPLIT,
          f"'any non-zero' split is {split} (reproduces the published 275/233)")

    # The threshold the Methods claim, which inverts the groups.
    for thr in cfg["ige_definition"]["sensitivity_thresholds"]:
        s = int((ige > thr).any(axis=1).sum())
        notes.append(f"threshold >={thr} ISU would give {s} sensitized / "
                     f"{len(ige) - s} non-sensitized")
        print(f"  NOTE  >= {thr} ISU -> {s} / {len(ige) - s}")

    nz = ige.values[ige.values > 0]
    print(f"  NOTE  {len(nz)} non-zero values, min {nz.min():.4f} ISU-E, "
          f"{(nz < 0.3).sum()} ({100 * (nz < 0.3).mean():.1f}%) below 0.3")

    ccd = cfg["ige_definition"]["ccd_marker"]
    if ccd in ige.columns:
        others = ige.drop(columns=[ccd])
        only = int(((ige[ccd] > 0) & others.eq(0).all(axis=1)).sum())
        print(f"  NOTE  {ccd}: {int((ige[ccd] > 0).sum())} positive, "
              f"{only} sensitized solely via it (CCD marker, not an allergen)")

    # ---- 4. Cluster assignments ----------------------------------------
    print("\n[4] Cluster assignments")
    cl = pd.read_csv(kora / "cluster_df_latent.csv", index_col=0)
    check(len(cl) == 275, f"cluster file has {len(cl)} individuals (expected 275)")
    sizes = cl["cluster"].value_counts().to_dict()
    check(sum(sizes.values()) == 275, "cluster sizes sum to 275")

    # ---- 5. Matched pair counts ----------------------------------------
    print("\n[5] Matched pair counts (read-only archive)")
    msum = pipe / "matching_summary.csv"
    check(msum.is_file(), "matching_summary.csv present")
    if msum.is_file():
        ms = pd.read_csv(msum).set_index("partition")["n_matched_pairs"].to_dict()
        bad = {k: (ms.get(k), v) for k, v in EXPECTED_PAIRS.items()
               if ms.get(k) != v}
        check(not bad, f"all 15 partitions match expected pair counts"
                       + (f" -- MISMATCH {bad}" if bad else ""))

    # ---- 6. Feature counts and result-tree completeness -----------------
    print("\n[6] Result tree completeness")
    levels = cfg["scope"]["levels"]
    prevs = cfg["scope"]["prevalence"]
    partitions = (cfg["scope"]["partitions"]["burden"]
                  + cfg["scope"]["partitions"]["clusters"])

    missing_src, missing_da, scope_rows = [], [], []
    for lvl in levels:
        for prev in prevs:
            for part in partitions:
                base = pipe / lvl / prev / part
                src = base / "source" / "filter_summary.csv"
                da = base / "da" / "linda_results.csv"
                if not src.is_file():
                    missing_src.append(str(src.relative_to(pipe)))
                    continue
                fs = pd.read_csv(src).iloc[0]
                if not da.is_file():
                    missing_da.append(str(da.relative_to(pipe)))
                scope_rows.append({
                    "level": lvl, "prevalence": prev, "partition": part,
                    "n_features": int(fs["n_features"]),
                    "n_ige": int(fs["n_ige"]), "n_ctrl": int(fs["n_ctrl"]),
                    "has_da": da.is_file(),
                })

    n_expected = len(levels) * len(prevs) * len(partitions)
    check(not missing_src,
          f"all {n_expected} source/ dirs present"
          + (f" -- {len(missing_src)} missing" if missing_src else ""))
    check(not missing_da,
          f"all {n_expected} linda_results.csv present"
          + (f" -- {len(missing_da)} missing" if missing_da else ""))

    scope = pd.DataFrame(scope_rows)
    if not scope.empty:
        got = (scope[scope.prevalence == "prev10"]
               .groupby("level")["n_features"].first().to_dict())
        bad = {k: (got.get(k), v) for k, v in EXPECTED_FEATURES_PREV10.items()
               if got.get(k) != v}
        check(not bad, "prev10 feature counts as expected (45/112/44/288)"
                       + (f" -- MISMATCH {bad}" if bad else ""))

        # Every partition must carry n_ige == n_ctrl (1:1 matching).
        unbal = scope[scope.n_ige != scope.n_ctrl]
        check(unbal.empty,
              "every partition is 1:1 balanced"
              + (f" -- {len(unbal)} unbalanced" if not unbal.empty else ""))

    # ---- 7. Covariates present with no missingness ---------------------
    print("\n[7] Covariates in the matched any_ige set")
    meta = pd.read_csv(kora / "kora.csv", low_memory=False)
    want = (cfg["covariates"]["base"] + cfg["covariates"]["extended"]
            + cfg["covariates"]["medication"])
    absent = [c for c in want if c not in meta.columns]
    check(not absent, f"all {len(want)} covariates present in kora.csv"
                      + (f" -- absent: {absent}" if absent else ""))

    ige_ids = pipe / "matched" / "any_ige" / "ige_ids.txt"
    ctrl_ids = pipe / "matched" / "any_ige" / "ctrl_ids.txt"
    if ige_ids.is_file() and ctrl_ids.is_file():
        ids = (ige_ids.read_text().split() + ctrl_ids.read_text().split())
        meta = meta.dropna(subset=["u3_16s_id"]).copy()
        meta["u3_16s_id"] = meta["u3_16s_id"].astype(int).astype(str)
        m = meta.set_index("u3_16s_id").reindex(ids)
        check(len(m) == 442, f"matched any_ige set has {len(m)} samples (expected 442)")
        have = [c for c in want if c in m.columns]
        miss = m[have].isna().sum()
        nonzero = miss[miss > 0].to_dict()
        check(not nonzero, "zero missingness across all covariates (442/442)"
                           + (f" -- missing: {nonzero}" if nonzero else ""))

    # ---- 8. Confirm no diet variable, no phylogeny ----------------------
    print("\n[8] Known blockers")
    d3 = "u3su_d3a1"
    if d3 in pd.read_csv(kora / "kora.csv", nrows=0, low_memory=False).columns:
        col = pd.read_csv(kora / "kora.csv", usecols=[d3],
                          low_memory=False)[d3].dropna().astype(str)
        top = col.value_counts().head(4).index.tolist()
        print(f"  NOTE  {d3} free text: {top} -> an intolerance/allergy "
              f"questionnaire, NOT dietary intake. No diet adjustment possible.")
    trees = [p for p in kora.iterdir() if p.suffix in (".nwk", ".qza", ".tre")]
    print(f"  NOTE  phylogeny files in kora_data: {len(trees)} "
          f"-> UniFrac and the published Faith's PD are blocked.")

    # ---- Write outputs -------------------------------------------------
    pd.DataFrame(manifest).to_csv(OUT / "manifest.csv", index=False)
    if not scope.empty:
        scope.to_csv(OUT / "scope.csv", index=False)
    print(f"\nWrote {OUT / 'manifest.csv'} and {OUT / 'scope.csv'}")

    print(f"\n=== {len(failures)} failure(s) ===")
    for f in failures:
        print(f"  - {f}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
