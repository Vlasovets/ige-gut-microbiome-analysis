#!/bin/bash
#SBATCH --job-name=rev_da_ext
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --array=0-179
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/daext_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/daext_%A_%a.err

# Extended differential abundance: ZicoSeq + paired Wilcoxon-CLR + limma-voom.
#
# 4 levels x 3 prevalence x 15 partitions = 180 tasks. ASV level IS included
# here (unlike module 02's wrapper) because none of these three methods carries
# ANCOM-BC2's pseudo_sens cost.
#
# MEASURED: species/prev05/any_ige (60 taxa, nperm=4999) = 1m21s wall, ~8 min
# CPU across ~6 cores. ZicoSeq dominates and scales roughly linearly in
# features, so asv/prev01 (1459 taxa) is expected around 35-45 min. 10h is
# deliberate headroom, not an estimate.
#
# PERMUTATION FLOOR - read before interpreting ZicoSeq at ASV rank.
#   A permutation p cannot go below 1/(nperm+1), so under BH over m features the
#   smallest attainable q is m/(nperm+1). At nperm=4999 that is 0.012 for m=60
#   (fine) but 0.29 for m=1459 (structurally incapable of significance).
#   The script prints this floor for every task and warns when it exceeds 0.05.
#   At ASV rank, interpret ZicoSeq via q_native (its own permutation FDR, which
#   is not floor-limited in the same way), NOT via q_within.
#   Raising nperm to clear the floor at m=1459 would need ~30,000 permutations;
#   that is a deliberate non-goal here, and the reason is recorded rather than
#   silently absorbed.
#
# Usage:  sbatch revision/slurm/02d_extended_da.sh

set -euo pipefail

REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"

LEVELS=(family genus species asv)
PREVS=(prev10 prev05 prev01)
PARTITIONS=(cluster_A cluster_B cluster_C cluster_D cluster_E cluster_F \
            cluster_G cluster_H cluster_I any_ige single multiple \
            high_score medium_score low_score)

N_PART=${#PARTITIONS[@]}          # 15
N_PREV=${#PREVS[@]}               # 3

I=${SLURM_ARRAY_TASK_ID}
PART=${PARTITIONS[$(( I % N_PART ))]}
PREV=${PREVS[$(( (I / N_PART) % N_PREV ))]}
LEVEL=${LEVELS[$(( I / (N_PART * N_PREV) ))]}

set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env
set -u

echo "=== $(date) === task ${I}: ${LEVEL} / ${PREV} / ${PART} ==="
Rscript revision/analysis/02d_extended_da.R \
  --level "${LEVEL}" --prevalence "${PREV}" --partition "${PART}" \
  --nperm 4999
echo "=== $(date) === Done ==="
