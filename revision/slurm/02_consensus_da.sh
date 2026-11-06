#!/bin/bash
#SBATCH --job-name=rev_da
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=0-134
#SBATCH --output=revision/slurm/logs/da_%A_%a.out
#SBATCH --error=revision/slurm/logs/da_%A_%a.err

# Consensus differential abundance across the analysis grid.
#
# 3 levels x 3 prevalence thresholds x 15 partitions = 135 tasks.
# ASV level is excluded here: 288-1459 features would make ANCOM-BC2 with
# pseudo_sens the dominant cost. Submit it separately with a longer wall time
# once the aggregated levels are in.
#
# qos=cpu_normal rather than preemptible: ALDEx2 at mc.samples=128 plus
# ANCOM-BC2 with sensitivity analysis is a ~2-20 min task per cell, and losing
# tasks to preemption mid-array leaves an incomplete grid that is easy to
# mistake for a real absence of results.
#
# Usage:  sbatch revision/slurm/02_consensus_da.sh

set -euo pipefail

REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"

LEVELS=(family genus species)
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

source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env

echo "=== $(date) === task ${I}: ${LEVEL} / ${PREV} / ${PART} ==="
Rscript revision/analysis/02_consensus_da.R \
  --level "${LEVEL}" --prevalence "${PREV}" --partition "${PART}" \
  --mc-samples 128
echo "=== $(date) === Done ==="
