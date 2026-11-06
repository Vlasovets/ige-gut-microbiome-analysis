# Offline review bundle: which file answers which comment

Open `00_start_here/REVIEWER_WORKBENCH.md` first: per category it gives the
reading key for each file, the computed verdict, manuscript line anchors and a
response skeleton. This index only tells you where things are.

No per-participant data is included (no per-sample alpha values, no matched
count tables, no sample IDs). Everything here is aggregate. Anything that needs
re-computation must be done on the cluster.

## Major comments

| Reviewer comment | Folder | Key files |
|---|---|---|
| R3: formal tests for Tables 1 and 2 | `01_cohort_tables` | `table1_full508_tested.csv`, `table2_clusters_tested.csv`; final LaTeX in `latex/` |
| R2: Table 1 p-values | `01_cohort_tables` | `latex/table1_main.tex` (as pasted in the manuscript) |
| R2: Table 2 tests with type I error control | `01_cohort_tables` | `latex/table2_main.tex`, `latex/tableS_cluster_tests.tex` |
| R3: confounders, matching | `01_cohort_tables` | `table1_matched442_qc.csv`, `latex/tableS_matched_cohort.tex` |
| R3: confounders, alpha diversity | `02_diversity/alpha` | `alpha_tests_long.csv` (raw counts primary; `analysis` column gives paired Wilcoxon / adjusted lm / depth-adjusted lm / depth balance), `alpha_forest_raw_primary.pdf`, `clusterF_depth_artifact.pdf` |
| R3: confounders, beta diversity | `02_diversity/beta` | `beta_<level>_prev10_<partition>.csv` (read the `ige_status` rows; M1 = sex+age+BMI, M2 adds lifestyle, M3 adds PPI), `betadisper_*.csv`, `pcoa_*.pdf` |
| R3: Faith's PD / UniFrac (supports alpha/beta) | `02_diversity/phylo_faithpd_unifrac` | `faithpd_*.csv`, `unifrac_*.csv` |
| R1: more DA methods (ANCOM-BC2, ALDEx2) | `03_differential_abundance` | `da_consensus_tiers.csv` (tier 1 = consensus), `da_method_rates.csv`, `da_method_upset.csv`, `da_method_effect_correlation.csv`, `consensus_heatmap_*.pdf`, `venn_consensus_primary.pdf` |
| R1 / R2: why the 61 published ASVs do not survive | `03_differential_abundance/why_61_became_1` | `nested_summary.csv`, `m0_seed_distribution.csv`, `hit_waterfall.csv`, `nested_m0_m3_waterfall.pdf` |
| R1: top-50 genera SPIEC-EASI networks | `04_networks` | `genus_top50_cohort.csv`, `per_partition/*/edges_*.csv`, rendered plots in `interactive/` |
| R1: overlay consensus taxa on networks | `04_networks` | `network_overlay_annotation.csv` (genus CONTAINING consensus-DA ASVs) |
| R2: networks underreported | `04_networks` | `per_partition/*/params_*.txt` (every NetCoMi argument), `network_edges_supplementary.csv` (all edges with permutation p and q) |
| R3: Random Forest benchmark | `05_classification` | `metrics/ml_species_prev05_any_ige_metrics.csv` (primary contrast, six models, per-model permutation p), `metrics/*_xgb_metrics.csv`, `ml_benchmark_primary.pdf` |
| R2: small subgroups, FDR control | across folders | every result table carries q-values; cluster-network partitions B/E were excluded by design; `clusterF_depth_artifact.pdf` |
| R1: species-level claims, V3-V4 | `06_discussion_support` | `TAXONOMY_COMPARISON.md`, `taxonomy_comparison_by_rank.csv` |
| R3: justify filtering cutoffs | `03_.../why_61_became_1`, `06_discussion_support` | `hit_waterfall.csv` (prevalence of each published hit), `genus_all_prev10_ranked.csv` |
| R3: vitamin-producing taxa in Discussion | `06_discussion_support` | `taxsea_recomputed.csv` (TaxSEA re-run: null) |
| R3: exact software versions | `07_methods_versions` | `KEY_RESOURCES.md`, `versions.json`, `slr-env.yml`, `*-requirements.txt` |

## Minor comments

| Comment | Where | Status |
|---|---|---|
| R2: Fig 1A all allergen categories | `06_discussion_support/fig1A_allergen_categories.*`, `lib/allergen_categories.csv` | figure exists; needs allergy co-author sign-off on category calls |
| R2: Fig 2B cluster colours | `06_discussion_support/lib/cluster_palette.csv` | palette only; figure NOT regenerated |
| R2: beta diversity association to IgE subgroups | `02_diversity/beta` | done, all 15 partitions |
| R2: Fig 3C volcano plateau | `03_differential_abundance/da_consensus_long.csv` | NO new volcano yet. The plateau came from plotting BH q; replot `-log10(p_raw)` from this file (filter `method == "LinDA"`) |
| R2: Fig 3D pie charts | `03_differential_abundance/fig3D_effect_distribution.*` | boxplot replacement exists |
| R2: burden cutoffs (page 200) | none | no analysis; needs a text justification |
| R2: "selected as predictors" (pages 329-337) | `05_classification/metrics/*_importance.csv` | text answer; importances are small (largest mean AUC drop 0.04; 1 of 60 species above 2 SD), consistent with chance AUC. *P. copri* ranks second, in line with the DA consensus |
| R2: Lachnospiraceae-Butyricicoccaceae claim (page 387) | `04_networks` | the old family-level network is superseded by the genus networks with permutation statistics |
| R2: DADA2 / taxonomy settings (6.4.3) | `07_methods_versions/UPSTREAM_METHODS.md` | recovered and verified (15,170/15,170 ASVs reproduced) |
| R2: SPIEC-EASI / NetCoMi defaults | `04_networks/per_partition/*/params_*.txt` | all arguments explicit |
| R3: reference style | none | answered in letter (Cell Press numbered style) |

## Also here

`08_plans_and_letter_drafts/`: the revision plan, the August response-letter draft (check its claims against the workbench chapter 7 before reusing any sentence), and the full `docs/` briefing. `09_code/`: every analysis module and SLURM wrapper, for the Methods text.
