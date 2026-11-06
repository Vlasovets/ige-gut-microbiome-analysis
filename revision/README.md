# Revision analyses — ISCIENCE-D-26-05498

Code for the iScience revision of *The role of IgE patterns and their link to the gut
microbiome in allergic sensitization*.

**The notebooks in `../notebooks/` are frozen and byte-identical to the submitted
version.** They are the record of what was published, including analyses this revision
supersedes. Nothing there is deleted or edited; the table below states what each
revision module replaces and why.

## What supersedes what

| Revision module | Supersedes | Why |
|---|---|---|
| `02_consensus_da.R` + `02b/02d/02e` | `notebooks/03_DA_analysis.ipynb` | The published analysis set `pair_matching = False` and randomly down-sampled controls, so age, sex and BMI never entered the model. The replacement is 1:1 matched and covariate-adjusted, and runs SIX DA methods on identical input (LinDA, ANCOM-BC2, ALDEx2, ZicoSeq, paired Wilcoxon-CLR, limma-voom) with a fixed consensus rule. |
| `05_beta_diversity.R` | — (new) | The submitted manuscript contains no beta-diversity analysis at all, despite framing itself against a prior beta-diversity result. |
| `05b_phylo_diversity.R` | — (new) | Faith's PD + UniFrac on the rebuilt phylogeny (reproduces the published ASV set 15,170/15,170). |
| `05c_alpha_diversity.R` | `notebooks/02_diversity.ipynb` | The published alpha diversity used an unmatched, unadjusted Mann-Whitney test with no multiplicity correction and no effort control. Replacement: matched pairs, BH within and across, RAW counts primary (no rarefaction, per McMurdie & Holmes 2014) with a depth-balance test and a log10(library-size) covariate model; repeated-rarefaction kept as labelled sensitivity only. |
| `03_genus_networks.R` + `03b_network_overlay.py` + `06q_network_core_overlay.py` + `06b_networks_interactive.R` | `notebooks/06_Network_analysis_with_NetCoMi.Rmd` | The published differential networks used `diffMethod = "threshold", diffThresh = 0, adjust = "none"`, which applies no statistical test, and compared clusters of n = 22–24 against the 233-control network, confounding edge loss with sample size. Rebuilt with per-edge permutation tests (nPerm 5,000, adaptBH) on partitions with n ≥ 74. Primary design `--design full` (all sensitized vs all 233 controls, the same samples as the DA; `results/networks_full`); `--design matched` (1:1 pairs, `results/networks`) is the sensitivity analysis. In the full design `medium_score`, `low_score` and `cluster_I` are unbalanced, so only per-edge tests, not edge counts, are compared there. 03b annotates nodes with the prevalence core (≥ 80% of the cohort) and the paper's DA (prev05) and tests core centrality; 06q draws the core-overlay figure; 06b is the RStudio exploration script (loads stored permutation statistics; never recomputes them). |
| `04_predictive_ml.py` + `04b/04c` wrappers | `notebooks/04_classification.ipynb` | The published evaluation filtered by prevalence before splitting and allowed matched pairs to straddle CV folds. Replacement: nested CV grouped on pairs; six model families incl. Random Forest and XGBoost, each with its own within-pair permutation null. |
| `07_taxsea_recompute.R` | external TaxSEA CSVs | TaxSEA was never run inside the repository, the input ranking was undocumented, and the exported table reports FDR < p in four of five rows. |
| `07a/07b/07c` (nested models) | — (new) | Re-runs the PUBLISHED design (M0) over 100 control draws and walks it to matched+adjusted (M3): the published 61-ASV count is the 55th percentile of its own sampling distribution (range 6–122), and most hits at every unadjusted step are sub-1%-prevalence artifacts. |
| `09_taxonomy_comparison.py` | — (new) | Independent, documented RESCRIPt SILVA 138.2 V3-V4 reference; species sparsity shown inherent to the amplicon. |
| `06c_review_figures.py`, `08_star_methods.py`, `10_export_markdown.py`, `11_workbench.py` | — (support) | Review figures; live KEY_RESOURCES/versions; the docs/ bundle; the reviewer workbench. |

## Environment

`slr-env` is **never mutated** — it produced the published results. Revision-only R
packages install into `env/rlibs/`, which each script prepends to `.libPaths()`.

⚠ **Do not use `conda create --clone` on this cluster — it deadlocks.**

```bash
conda activate slr-env
sbatch revision/slurm/00_install_r_extras.sh     # ANCOMBC, ALDEx2, TaxSEA, UpSetR
```

## Running

All paths and every pre-specified choice live in `config/config.yaml`. No script may
override it.

```bash
python revision/analysis/00a_verify_inputs.py        # manifest + scope assertions
python  revision/analysis/00b_integrity_audit.py     # gate: what is actually wrong
Rscript revision/analysis/01_baseline_stats.R        # Tables 1 & 2 with tests
```

Execution order and SLURM shapes are documented in
`../../REVISION_PLAN.md` §8.

## Pre-specified scope — read before interpreting anything

BH correction across all 15 partitions leaves **nothing** significant at any rank
(minimum q = 0.2367 species, 0.4170 genus, 0.6073 family, 0.9809 ASV). A primary
contrast is therefore declared in `config/config.yaml` *before* results are inspected:
`any_ige`, species rank, 5% prevalence, 221 matched pairs. Every other partition, rank
and threshold is secondary or exploratory, and every results table carries
`q_within`, `q_secondary` and `q_all_partitions` so the multiplicity is visible.

## Data

Inputs are read-only and live outside this repository; see `DATA.md`.
`../../IGE_Causal_Pipeline/results/` is an **unversioned input archive** and is never
written to.
