# Data inputs and provenance

Individual-level data is **not** in this repository and must not be committed to it
(KORA participant consent). This file records what every analysis reads, where it
comes from, and which previously-referenced files do not exist.

## Primary inputs — read-only

Location: `../KORA_data/` (outside this repository), set as `paths.kora_data` in
`revision/config/config.yaml`.

| File | Shape | Contents |
|---|---|---|
| `feature_table.tsv` | 15,170 ASVs × 2,034 samples | ASV count table |
| `taxonomy_clean.csv` | 15,170 rows | ASV taxonomy (domain → species); a *cleaned* derivative whose cleaning script is not in the repository |
| `kora.csv` | 76 columns | KORA FF4 metadata |
| `ige_kora_subset_all.csv` | 508 × **104** | Allergen-specific IgE, ISU-E |
| `cluster_df_latent.csv` | 275 rows | IgE cluster assignments A–I |

**The IgE panel has 104 components, not 112.** 112 is the ISAC platform's component
count; 104 is what remains after total IgE/IgG and all-zero components are dropped, and
104 is what every analysis uses.

`revision/results/manifest.csv` records a SHA-256 for each of these, written by
`00a_verify_inputs.py`. Any downstream number can be traced to a specific input state.

## Derived inputs — read-only archive

Location: `../IGE_Causal_Pipeline/results/`, set as `paths.pipeline_results`.

Matched sample sets and per-partition aggregated tables for
4 taxonomic ranks × 3 prevalence thresholds × 15 partitions = 180 combinations.

⚠ **This directory is not under version control** (no `.git`) and is the provenance of
the numbers the revision reports. It is treated as read-only: no revision script writes
into it. That is also why the revision code lives here rather than there — code in an
unversioned directory cannot be cited.

## Vendored code

`utils/` was imported by the notebooks but absent from the submitted repository. It is
now vendored from `../Causal_Sparse_Low_Rank_Microbiome_Tutorial/utils/`; see
`utils/PROVENANCE.md`. It contains the CLR/mCLR transform and the pair-matching
implementation, both methodologically load-bearing and previously unavailable to
reviewers.

## Files the notebooks read but nothing writes

These were read by `notebooks/01`–`06` and are not produced by any committed code, so
the published pipeline could not be re-run end to end. They are listed here rather than
quietly dropped.

| File | Read by | Status |
|---|---|---|
| `data/ige_clustered_275.csv` | 01 | not produced by any notebook |
| `data/latent_allergy_components.csv` | 01 | read in cell 17, but the PCA that would create it runs later (cells 37/43) and never saves |
| `data/cluster_df_latent.csv` | 01, 03, 04, 05 | recovered — present in `KORA_data/` |
| `data/control_samples.csv` | 01, 04, 05 | not produced; recomputable as "all allergens = 0" |
| `data/taxa_dict_ige.pickle` | 02 | not produced |
| `data/asv_ige_508.csv`, `data/54subset_latent.csv`, `data/co_occurrence_matrix.csv` | 05 | not produced |
| `data/phylum_labels_only_taxa.csv` | 06 (Rmd) | not produced; replaced by `revision/results/taxa_sets/tax_table_genus.csv` |
| `data/TaxSEA_Metabolite_Producers_{genus,species}.csv` | 03 | **not produced and not on disk** — the enrichment analysis was unreproducible, which is why module 07 re-runs it from a documented ranking |
| `data/q2/tree/rooted_tree_15170/tree.nwk` | 02 | **absent** — see below |

## Not available, and what follows

**~~No phylogeny exists~~ — being rebuilt.** The source reads were located at
`/lustre/groups/bds01/datasets/KORA/paired-end-demux.qza` (8.45 GB, UUID `c5a66171`,
`SampleData[PairedEndSequencesWithQuality]`), and the preprocessing notebooks record the
exact commands. `revision/slurm/00d_rebuild_phylogeny.sh` re-runs DADA2 under
**qiime2-2021.11** — the version the manuscript states — and rebuilds the tree with
`align-to-tree-mafft-fasttree`.

Matching the QIIME 2 version is the point: an ASV id is the MD5 of its sequence, so
identical sequences give identical ids, and the rebuilt tree's tips then line up with the
existing `feature_table.tsv`. The job verifies that overlap and declines to endorse the
tree if coverage is poor. Once it lands, **UniFrac and Faith's PD become computable**.
Full provenance in `revision/UPSTREAM_METHODS.md`.

**No dietary variable exists.** `u3su_d3*` is a self-reported intolerance/allergy
questionnaire — its free-text fields contain "Lactose", "Gluten", "Penicillin" — not
dietary intake, and `-999` is a skip-logic sentinel. Diet cannot be adjusted for. Do not
construct a diet index from these columns.

**~~DADA2 parameters are not recorded~~ — RECOVERED.** The preprocessing notebooks
(`KORA_data/preprocessing/`) record every setting; see `revision/UPSTREAM_METHODS.md`.
In brief: `--p-trim-left-f 17 --p-trim-left-r 21 --p-trunc-len-f 245 --p-trunc-len-r 245
--p-max-ee-f 2 --p-max-ee-r 2 --p-trunc-q 2`, chimera method `consensus`, classifier
`silva-138-99-nb-classifier.qza` from the QIIME 2 2021.4 resources.

Two caveats that follow from it, both now stated in the Methods rather than left
implicit: primers were removed by **positional trimming, not cutadapt** (the trim
lengths are exactly the primer lengths), and the classifier is the **full-length**
SILVA 138 artifact rather than one extracted to V3–V4 — which is the mechanism behind
96% of ASVs being unassigned at species level.

**Still open:** the notebook holds both a `classify-sklearn` and a RESCRIPt
classification output, and which produced `taxonomy_clean.csv` cannot be determined from
it. The exact SILVA 138 point release is likewise not recorded in the artifact name.

## Outputs

Written under `revision/results/`. Aggregate tables and figures are tracked; anything
individual-level is excluded by `revision/.gitignore`.
