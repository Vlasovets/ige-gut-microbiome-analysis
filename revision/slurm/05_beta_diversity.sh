#!/bin/bash
#SBATCH --job-name=rev_beta
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=0-44
#SBATCH --output=revision/slurm/logs/beta_%A_%a.out
#SBATCH --error=revision/slurm/logs/beta_%A_%a.err

# Beta diversity / PERMANOVA / betadisper across the IgE subgroups.
#
# Scope: 3 levels x 15 partitions at prev10 = 45 tasks. prev10 is the primary
# prevalence threshold. (An earlier version of this header claimed any_ige was
# also run at prev05/prev01 by hand; no such outputs exist, so the claim is
# removed.) ASV rank is covered by 05_beta_diversity_asv.sh.
#
# Timing measured on family/prev10/any_ige: 40 min for 4 distances x 7 model-by
# combinations at nperm=5000. Marginal (type III) tests dominate - they refit per
# term, and cost ~3.5x a sequential test. 4 h gives headroom for genus (112 taxa).
#
# qos=cpu_normal, not preemptible: a 40+ min task killed mid-array leaves an
# incomplete grid that reads like a real absence of results.
#
# Usage:  sbatch revision/slurm/05_beta_diversity.sh

set -euo pipefail

REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"

LEVELS=(family genus species)
PARTITIONS=(cluster_A cluster_B cluster_C cluster_D cluster_E cluster_F \
            cluster_G cluster_H cluster_I any_ige single multiple \
            high_score medium_score low_score)
N_PART=${#PARTITIONS[@]}

I=${SLURM_ARRAY_TASK_ID}
PART=${PARTITIONS[$(( I % N_PART ))]}
LEVEL=${LEVELS[$(( I / N_PART ))]}

source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env

echo "=== $(date) === task ${I}: ${LEVEL} / prev10 / ${PART} ==="
Rscript revision/analysis/05_beta_diversity.R \
  --level "${LEVEL}" --prevalence prev10 --partition "${PART}" --nperm 5000
echo "=== $(date) === Done ==="
