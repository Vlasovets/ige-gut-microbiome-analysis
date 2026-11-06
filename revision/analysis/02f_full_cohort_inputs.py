#!/usr/bin/env python
"""
Module 02f - inputs for the FULL-COHORT, covariate-adjusted DA sensitivity run.

WHY
  Modules 02/02d test the 1:1 matched subset (221 pairs for any_ige) with
  ~ ige_status + sex + age + bmi. The alternative way to control the same
  confounders is adjustment alone on every available sample: all sensitized
  members of a partition vs all 233 controls. It keeps the 54+ sensitized
  individuals that matching drops, so it is the better-powered of the two.
  Running both and comparing is the sensitivity analysis; neither design is
  unadjusted.

WHAT IS HELD FIXED (so the only change is the sample set)
  * Feature set: aggregation and prevalence filter use the pipeline's own
    functions on the full feature table, label-blind, exactly as
    IGE_Causal_Pipeline/analysis/00b_aggregate_otu.py does. Asserted equal to
    the matched source columns.
  * Counts: for every matched sample, the rebuilt row must equal the pipeline's
    matched source row cell for cell.
  * Membership: results/nested_models/designs/base.csv (module 07a, validated
    against the published partition sizes).

OUTPUT (per-participant rows; gitignored, regenerated in ~1 min)
  results/full_cohort/inputs/<level>/<prev>/<partition>/source/
      otu_ige.csv  otu_non_ige.csv  sample_data_ige.csv  sample_data_non_ige.csv
  sample_data carries only u3csex, u3talter, u3tbmi - what the DA models use.
  results/full_cohort/inputs_summary.csv  (counts only, tracked)

Usage:  python revision/analysis/02f_full_cohort_inputs.py
        python revision/analysis/02f_full_cohort_inputs.py --unfiltered   # prev00, any_ige only
"""

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import yaml

REVISION = Path(__file__).resolve().parent.parent
cfg = yaml.safe_load((REVISION / "config" / "config.yaml").read_text())
KORA = Path(cfg["paths"]["kora_data"])
PIPE = Path(cfg["paths"]["pipeline_results"])
PIPE_AN = Path(cfg["paths"]["pipeline_analysis"])
OUT = REVISION / "results" / "full_cohort"
LEVELS = cfg["scope"]["levels"]
PREVS = cfg["scope"]["prevalence"]
PARTS = cfg["scope"]["partitions"]["burden"] + cfg["scope"]["partitions"]["clusters"]
COVS = ["u3csex", "u3talter", "u3tbmi"]

# Import the pipeline's aggregation helpers instead of re-implementing them.
spec = importlib.util.spec_from_file_location("agg", PIPE_AN / "00b_aggregate_otu.py")
agg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(agg)

print("=== 02f: full-cohort DA inputs ===")
feat = pd.read_csv(KORA / "feature_table.tsv", sep="\t", index_col=0)
feat.columns = feat.columns.astype(str)
feat.index = feat.index.astype(str)
tax = pd.read_csv(KORA / "taxonomy_clean.csv", index_col=0)
tax.index = tax.index.astype(str)

base = pd.read_csv(REVISION / "results" / "nested_models" / "designs" / "base.csv",
                   dtype=str)
ctrl = base.loc[base.partition == "__controls__", "sample_id"].tolist()
assert len(ctrl) == 233, len(ctrl)

kora = pd.read_csv(KORA / "kora.csv", usecols=["u3_16s_id"] + COVS, low_memory=False)
kora["u3_16s_id"] = kora["u3_16s_id"].astype("Int64").astype(str)
kora = kora.drop_duplicates("u3_16s_id").set_index("u3_16s_id")
cohort = set(base.sample_id)
missing = sorted(cohort - set(kora.index))
assert not missing, f"{len(missing)} cohort samples absent from kora.csv"
assert cohort <= set(feat.columns), "cohort samples absent from feature table"
print(f"  cohort {len(cohort)} samples; covariates missing: "
      + ", ".join(f"{c}={int(kora.loc[list(cohort), c].isna().sum())}" for c in COVS))

# --unfiltered: the submitted manuscript's feature set ("prev00") = every feature
# non-zero in at least one sample of the comparison, no prevalence filter.
# Only the primary contrast is built; the pipeline has no prev00 source to
# compare against, so the anchor is the published count of 6,958 ASVs.
if "--unfiltered" in sys.argv:
    part = "any_ige"
    cases = base.loc[(base.partition == part) & (base.role == "case"), "sample_id"].tolist()
    rows = []
    for level in LEVELS:
        lab = agg._build_label_map(tax, level)
        f = feat[cases + ctrl].copy()
        f["_label"] = [lab.get(a) for a in f.index]
        f = f.dropna(subset=["_label"]).groupby("_label").sum()
        fa = f.loc[(f > 0).any(axis=1)].T
        if level == "asv":
            assert fa.shape[1] == 6958, f"unfiltered ASVs = {fa.shape[1]}, published 6,958"
        d = OUT / "inputs" / level / "prev00" / part / "source"
        d.mkdir(parents=True, exist_ok=True)
        fa.loc[cases].to_csv(d / "otu_ige.csv")
        fa.loc[ctrl].to_csv(d / "otu_non_ige.csv")
        kora.loc[cases, COVS].rename_axis("u3_16s_id").to_csv(d / "sample_data_ige.csv")
        kora.loc[ctrl, COVS].rename_axis("u3_16s_id").to_csv(d / "sample_data_non_ige.csv")
        rows.append({"level": level, "prevalence": "prev00", "partition": part,
                     "n_features": fa.shape[1], "n_ige": len(cases), "n_ctrl": len(ctrl)})
        print(f"  {level}/prev00/{part}: {fa.shape[1]} features")
    pd.DataFrame(rows).to_csv(OUT / "inputs_summary_unfiltered.csv", index=False)
    sys.exit(0)

summary = []
for level in LEVELS:
    lab = agg._build_label_map(tax, level)
    f = feat.copy()
    f["_label"] = [lab.get(a) for a in f.index]
    f = f.dropna(subset=["_label"]).groupby("_label").sum()
    for prev in PREVS:
        thr = agg.PREV_VALUES[prev]
        fa = f.loc[(f > 0).mean(axis=1) >= thr].T          # samples x features
        for part in PARTS:
            cases = base.loc[(base.partition == part) & (base.role == "case"),
                             "sample_id"].tolist()
            # QC 1+2: identical features, identical counts on the matched samples
            src = PIPE / level / prev / part / "source"
            ref_ige = pd.read_csv(src / "otu_ige.csv", index_col=0)
            ref = pd.concat([ref_ige, pd.read_csv(src / "otu_non_ige.csv", index_col=0)])
            ref.index = ref.index.astype(str)
            assert list(ref.columns) == list(fa.columns), f"{level}/{prev}: feature set differs"
            assert set(ref.index) <= set(cases) | set(ctrl), f"{part}: matched sample outside cohort"
            assert (fa.loc[ref.index, ref.columns].values == ref.values).all(), \
                f"{level}/{prev}/{part}: counts differ from pipeline source"

            d = OUT / "inputs" / level / prev / part / "source"
            d.mkdir(parents=True, exist_ok=True)
            fa.loc[cases].to_csv(d / "otu_ige.csv")
            fa.loc[ctrl].to_csv(d / "otu_non_ige.csv")
            kora.loc[cases, COVS].rename_axis("u3_16s_id").to_csv(d / "sample_data_ige.csv")
            kora.loc[ctrl, COVS].rename_axis("u3_16s_id").to_csv(d / "sample_data_non_ige.csv")
            summary.append({"level": level, "prevalence": prev, "partition": part,
                            "n_features": fa.shape[1], "n_ige": len(cases),
                            "n_ctrl": len(ctrl), "n_matched_ige": len(ref_ige)})
        print(f"  {level}/{prev}: {fa.shape[1]} features, 15 partitions written, QC pass")

pd.DataFrame(summary).to_csv(OUT / "inputs_summary.csv", index=False)
print(f"Wrote {len(summary)} input sets + {OUT / 'inputs_summary.csv'}")
