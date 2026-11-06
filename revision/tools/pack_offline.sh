#!/bin/bash
# Build an offline review bundle, organised by reviewer topic.
#
# INCLUDED: aggregate results (per taxon, per test, per analysis), figures,
# LaTeX tables, docs, analysis code and environment pins.
# EXCLUDED, deliberately: anything with one row per participant
#   - results/alpha/alpha_persample.csv
#   - results/phylo/*_persample.csv
#   - results/nested_models/designs/  (sample IDs per draw)
#   - results/networks/*/netprops_*.rds  (embed matched per-sample counts)
#   - the KORA source tables
# plus regenerable bulk: networks/*/assoPerm.bmat (~200 MB each) and the
# per-analysis DA files (their merged form is in tables/).
#
# Usage:  bash revision/tools/pack_offline.sh
set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
WS=/home/itg/oleg.vlasovets/slr_example
R=${REPO}/revision/results
STAMP=$(date '+%Y%m%d')
NAME="ige_revision_offline_${STAMP}"
STAGE="/ictstr01/scratch/users/${USER}/pack_${STAMP}/${NAME}"
OUTDIR="/ictstr01/scratch/users/${USER}/revision_archive"
rm -rf "${STAGE}"; mkdir -p "${STAGE}" "${OUTDIR}"
cd "${REPO}"

c() { mkdir -p "${STAGE}/$1"; shift_dir="$1"; shift; for f in "$@"; do [ -e "$f" ] && cp -r "$f" "${STAGE}/${shift_dir}/"; done; }

# 00 start here
c 00_start_here docs/REVIEWER_WORKBENCH.md docs/README.md
# 01 cohort tables (R3 major 1, R2 Table 1 / Table 2)
c 01_cohort_tables ${R}/tables/table1_full508_tested.csv ${R}/tables/table1_matched442_qc.csv \
  ${R}/tables/table2_clusters_tested.csv ${R}/tables/table2_clusters_descriptives.csv \
  ${R}/tables/tableS1_partitions_tested.csv ${R}/tables/group_descriptives_pairs_topIgE.csv \
  ${R}/tables/latex
# 02 confounders: alpha + beta diversity (R3 major 2, R2 beta minor)
c 02_diversity/alpha ${R}/alpha/alpha_tests_long.csv ${R}/alpha/alpha_summary_tableS5.csv \
  ${R}/figures/review/alpha_forest_raw_primary.* ${R}/figures/review/clusterF_depth_artifact.*
c 02_diversity/phylo_faithpd_unifrac $(ls ${R}/phylo/*.csv | grep -v persample)
c 02_diversity/beta ${R}/beta ${R}/figures/pcoa_*
# 03 differential abundance (R1 major; R2 Fig 3C/3D minor)
c 03_differential_abundance ${R}/tables/da_consensus_tiers.csv ${R}/tables/da_method_rates.csv \
  ${R}/tables/da_method_upset.csv ${R}/tables/da_method_overlap_pairwise.csv \
  ${R}/tables/da_method_effect_correlation.csv ${R}/tables/da_method_overlap_report.txt \
  ${R}/tables/da_consensus_long.csv ${R}/tables/da_all_methods_long.csv \
  ${R}/figures/consensus_heatmap_* ${R}/figures/venn_consensus_primary.* \
  ${R}/figures/fig3D_effect_distribution.*
c 03_differential_abundance/why_61_became_1 ${R}/nested_models/nested_summary.csv \
  ${R}/nested_models/m0_seed_distribution.csv ${R}/nested_models/hit_waterfall.csv \
  ${R}/nested_models/M1 ${R}/nested_models/M2 ${R}/nested_models/M3 \
  ${R}/nested_models/designs/validation.txt ${R}/figures/review/nested_m0_m3_waterfall.*
# 04 networks (R1 major; R2 underreported)
c 04_networks ${R}/networks/network_edges_supplementary.csv ${R}/networks/network_overlay_annotation.csv \
  ${R}/taxa_sets/genus_top50_cohort.csv ${R}/taxa_sets/manifest.json ${R}/figures/interactive
for p in any_ige multiple medium_score low_score cluster_I; do
  c 04_networks/per_partition/${p} ${R}/networks/${p}/edges_${p}.csv ${R}/networks/${p}/nodes_${p}.csv \
    ${R}/networks/${p}/params_${p}.txt
done
# 05 predictive modelling (R3 major; R2 pages 329-337)
c 05_classification ${R}/figures/review/ml_benchmark_primary.*
c 05_classification/metrics $(ls ${R}/ml/*.csv)
c 05_classification/figures_per_cell $(ls ${R}/figures/ml_*)
# 06 discussion: species claims, cutoffs, TaxSEA, Fig 1A (R1, R2, R3)
c 06_discussion_support revision/TAXONOMY_COMPARISON.md ${R}/tables/taxonomy_comparison_by_rank.csv \
  ${R}/tables/taxonomy_comparison_per_asv.csv ${R}/tables/taxsea_recomputed.csv \
  ${R}/tables/taxsea_recomputed_params.txt ${R}/taxa_sets/genus_all_prev10_ranked.csv \
  ${R}/figures/fig1A_allergen_categories.* revision/lib
# 07 methods and versions (R3 versions; R2 DADA2 / SPIEC-EASI settings)
c 07_methods_versions ${R}/tables/KEY_RESOURCES.md revision/UPSTREAM_METHODS.md revision/DATA.md \
  revision/env/versions.json revision/env/slr-env.yml revision/env/xgb-requirements.txt \
  revision/env/legacy-requirements.txt revision/config/config.yaml
# 08 planning docs and code
c 08_plans_and_letter_drafts "${WS}/REVISION_PLAN.md" "${WS}/REVIEWER_RESPONSE_OUTLINE.md" docs
c 09_code revision/analysis revision/slurm revision/tools revision/README.md
rm -rf "${STAGE}/09_code/slurm/logs"

cp "${REPO}/revision/tools/OFFLINE_INDEX.md" "${STAGE}/INDEX.md"

# Guard: nothing per-participant may be staged.
if grep -rlE '^"?(sample|u3_16s_id|sample_id)"?,' "${STAGE}" --include='*.csv' | grep -q .; then
  echo "FATAL: per-sample file staged:"; grep -rlE '^"?(sample|u3_16s_id|sample_id)"?,' "${STAGE}" --include='*.csv'
  exit 1
fi
find "${STAGE}" -name '*.rds' -o -name '*.bmat' | grep -q . && { echo "FATAL: rds/bmat staged"; exit 1; }

cd "$(dirname "${STAGE}")"
rm -f "${OUTDIR}/${NAME}.zip"
zip -qr "${OUTDIR}/${NAME}.zip" "${NAME}"
rm -rf "$(dirname "${STAGE}")"
echo "archive: ${OUTDIR}/${NAME}.zip"
ls -lh "${OUTDIR}/${NAME}.zip"
unzip -l "${OUTDIR}/${NAME}.zip" | tail -1
