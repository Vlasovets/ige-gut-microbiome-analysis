#!/bin/bash
#SBATCH --job-name=rev_beta_asv
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment='agent:claude-code; ASV-level PERMANOVA/betadisper, completes the beta table'
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=0-14
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/beta_asv_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/beta_asv_%A_%a.err

# Beta diversity at ASV rank: Bray-Curtis, Jaccard, Aitchison, robust Aitchison
# - the four non-phylogenetic distances 05_beta_diversity.sh ran only at
# family/genus/species. prev01 is used so the feature set (1,459 ASVs for
# any_ige) is IDENTICAL to the UniFrac analysis (05b), letting the ASV row
# compare phylogenetic and non-phylogenetic distances on the same ASVs.
#
# Same script, unchanged; it is level-agnostic. Sizing: previous array max task
# 52 min, MaxRSS 0.3 GB at 2 CPU / 16G (job 39121154). PERMANOVA cost scales
# with n (permutations), not with the number of taxa, so the same request holds.
set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"
PARTITIONS=(cluster_A cluster_B cluster_C cluster_D cluster_E cluster_F \
            cluster_G cluster_H cluster_I any_ige single multiple \
            high_score medium_score low_score)
PART=${PARTITIONS[${SLURM_ARRAY_TASK_ID}]}
set +u; source ~/miniconda3/etc/profile.d/conda.sh; conda activate slr-env; set -u
echo "=== $(date) === asv / prev01 / ${PART} ==="
Rscript revision/analysis/05_beta_diversity.R \
  --level asv --prevalence prev01 --partition "${PART}" --nperm 5000
echo "=== $(date) === Done ==="
