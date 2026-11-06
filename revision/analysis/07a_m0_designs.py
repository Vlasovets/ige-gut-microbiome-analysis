#!/usr/bin/env python
"""
Module 07a - reproduce the PUBLISHED (M0) sample designs, across 100 seeds.

Runs under revision/env/legacy-venv (pandas 2.2.3), because the published
notebook (notebooks/03_DA_analysis.ipynb) was written against pandas 2-era
semantics and the reproduction must not silently inherit pandas 3 behaviour.

WHAT THE PUBLISHED DESIGN ACTUALLY WAS (read from the notebook, cell by cell):
  * Sensitized (275) = rows of ige_kora_subset_all.csv with ANY non-zero
    component; controls (233) = rows with none.
  * Partitions: any_ige/single/multiple from non-zero counts; high/medium/low
    from the burden score (count > 15, (5, 15], <= 5) computed on the 275;
    clusters A-I from cluster_df_latent.csv.
  * Down-sampling: for every partition EXCEPT any_ige, controls were
    down-sampled to the partition size with pandas .sample(n, random_state=42)
    - a single draw. any_ige used all 233 controls.
  * Feature filtering happened AFTER sample selection, per comparison:
    keep ASVs non-zero within the selected samples (freq_threshold=0 in
    utils/preprocessing.filter_and_process_asv_table - the docstring mentions
    1% but the code compares prevalence > freq_threshold, which is 0), then
    drop samples with zero mean.
  * LinDA ~w, prev_cut=0, lib_cut=1; BH per partition x rank at alpha 0.05.
    The published "61 ASVs" is the UNION of ASV-rank hits across partitions.

REPRODUCTION CAVEAT, recorded rather than hidden: the notebook drew controls
from control_samples.csv, whose ROW ORDER we do not have; .sample() with a
fixed seed depends on row order. We use the row order of ige_kora_subset_all,
so the seed-42 draw here is A realization of the published design, not
necessarily THE published draw. That is why 100 seeds are run: the claim rests
on the distribution of hit counts, with the published 61 as one realization.

Outputs (revision/results/nested_models/designs/):
  base.csv               partition membership + full control list (M1)
  seed_<s>.csv           per-seed down-sampled control draw per partition (M0)
  validation.txt         partition-size checks against the published counts

Usage:  revision/env/legacy-venv/bin/python revision/analysis/07a_m0_designs.py
"""

import sys
from pathlib import Path

import pandas as pd
import yaml

REVISION = Path(__file__).resolve().parent.parent
cfg = yaml.safe_load((REVISION / "config" / "config.yaml").read_text())
KORA = Path(cfg["paths"]["kora_data"])
OUT = REVISION / "results" / "nested_models" / "designs"
OUT.mkdir(parents=True, exist_ok=True)

SEEDS = list(range(42, 142))          # 100 seeds; 42 first = the published seed
EXPECT = {"any_ige": 275, "single": 47, "multiple": 228,
          "high_score": 57, "medium_score": 110, "low_score": 108}

print("=== 07a: M0 design extraction ===")
print(f"  pandas {pd.__version__} (legacy-venv)")

ige = pd.read_csv(KORA / "ige_kora_subset_all.csv", index_col=0, low_memory=False)
ige.index = ige.index.astype(str)
nz = (ige != 0).sum(axis=1)
sens = ige.index[nz >= 1].tolist()
ctrl = ige.index[nz == 0].tolist()          # row order of ige_kora_subset_all

# burden score = COUNT of positive components, computed on the sensitized
score = nz[nz >= 1]
partitions = {
    "any_ige": sens,
    "single": nz.index[nz == 1].tolist(),
    "multiple": nz.index[nz > 1].tolist(),
    "high_score": score.index[score > 15].tolist(),
    "medium_score": score.index[(score <= 15) & (score > 5)].tolist(),
    "low_score": score.index[score <= 5].tolist(),
}
cl = pd.read_csv(KORA / "cluster_df_latent.csv", index_col=0, low_memory=False)
cl.index = cl.index.astype(str)
for c in sorted(cl["cluster"].unique()):
    partitions[f"cluster_{c}"] = cl.index[cl["cluster"] == c].tolist()

# ── validation against the published sizes ────────────────────────────────
lines = [f"controls: {len(ctrl)} (expect 233)"]
ok = len(ctrl) == 233
for p, exp in EXPECT.items():
    got = len(partitions[p])
    ok &= (got == exp)
    lines.append(f"{p}: {got} (expect {exp}) {'OK' if got == exp else 'MISMATCH'}")
for p in sorted(k for k in partitions if k.startswith("cluster_")):
    lines.append(f"{p}: {len(partitions[p])}")
print("  " + "\n  ".join(lines))
(OUT / "validation.txt").write_text("\n".join(lines) + "\n")
if not ok:
    print("FATAL: partition sizes do not reproduce the published design")
    sys.exit(1)

# ── base design (memberships + the M1 full-control comparison) ────────────
rows = [{"partition": p, "sample_id": s, "role": "case"}
        for p, m in partitions.items() for s in m]
rows += [{"partition": "__controls__", "sample_id": s, "role": "control_pool"}
         for s in ctrl]
pd.DataFrame(rows).to_csv(OUT / "base.csv", index=False)

# ── per-seed M0 control draws ─────────────────────────────────────────────
# EXACTLY the notebook's operation: ctlr.sample(n=class_size, random_state=seed),
# applied per partition, skipped for any_ige. `ctrl_df` mirrors the notebook's
# ctlr_samples (a DataFrame indexed by sample id, in control-file row order).
ctrl_df = pd.DataFrame(index=pd.Index(ctrl, name="u3_16s_id"))
for s in SEEDS:
    rec = []
    for p, members in partitions.items():
        if p == "any_ige":
            continue                      # published design: no down-sampling
        draw = ctrl_df.sample(n=len(members), random_state=s).index.tolist()
        rec += [{"partition": p, "sample_id": x, "role": "control"} for x in draw]
    pd.DataFrame(rec).to_csv(OUT / f"seed_{s}.csv", index=False)
print(f"  wrote base.csv + {len(SEEDS)} seed files to {OUT}")
print("Done.")
