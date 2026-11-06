#!/bin/bash
#SBATCH --job-name=rev_da_asv
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-44
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/da_asv_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/da_asv_%A_%a.err

# Consensus DA at ASV rank - the level 02_consensus_da.sh deliberately excluded.
#
# WHY THIS MATTERS
#   The published claim was 61 differentially abundant ASVs. The matched,
#   covariate-adjusted re-analysis reports a null at ASV rank - but from LinDA
#   ALONE. Answering a "use more than one DA method" request with a
#   single-method number, at exactly the rank where the retracted claim lived,
#   is the weakest point in the response letter. This closes it.
#
# COST
#   3 prevalence x 15 partitions = 45 tasks. Features per task: 288 (prev10),
#   489 (prev05), 1459 (prev01) - versus 44-188 at the aggregated levels, where
#   measured cost was 2.5 min/task at 4 CPU/16G (job 39116423, MaxRSS 1.51 GiB).
#   ANCOM-BC2 with pseudo_sens grows faster than linearly in features, so
#   resources are raised to 8 CPU / 64G / 12h. Generous on purpose: an OOM
#   across 45 tasks costs far more than an over-request.
#
#   Gated behind a single smoke test on the heaviest cell (see
#   02c_smoke_asv.sh) via --dependency=afterok, so a bad config cannot burn the
#   whole array.
#
# Usage:
#   SMOKE=$(sbatch --parsable revision/slurm/02c_smoke_asv.sh)
#   sbatch --dependency=afterok:${SMOKE} revision/slurm/02c_consensus_da_asv.sh

set -euo pipefail

REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"

PREVS=(prev10 prev05 prev01)
PARTITIONS=(cluster_A cluster_B cluster_C cluster_D cluster_E cluster_F \
            cluster_G cluster_H cluster_I any_ige single multiple \
            high_score medium_score low_score)

N_PART=${#PARTITIONS[@]}          # 15

I=${SLURM_ARRAY_TASK_ID}
PART=${PARTITIONS[$(( I % N_PART ))]}
PREV=${PREVS[$(( I / N_PART ))]}

set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env
set -u

echo "=== $(date) === task ${I}: asv / ${PREV} / ${PART} ==="
Rscript revision/analysis/02_consensus_da.R \
  --level asv --prevalence "${PREV}" --partition "${PART}" \
  --mc-samples 128
echo "=== $(date) === Done ==="
