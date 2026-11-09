#!/bin/bash
#SBATCH --job-name=rev_ml_perm
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment='agent:claude-code; per-model permutation nulls, family/genus prev10 grid'
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00

#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/mlperm_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/mlperm_%A_%a.err

# Per-model permutation nulls for L1/L2/RF/HistGB on the 30 cells the August run
# covered (family+genus x prev10 x 15 partitions; that run had an L2-only null).
# Same seeds/folds, so AUCs reproduce. --tag perm so the CSVs
# behind already-drafted numbers are never rewritten.
# Est. RF null dominates (200 perms x 5 folds x 1000 trees); 8 h wall is headroom.
# Usage: sbatch --array=0-29 revision/slurm/04e_predictive_ml_perm_backfill.sh
set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"
PY=${REPO}/revision/env/xgb-venv/bin/python

PARTITIONS=(cluster_A cluster_B cluster_C cluster_D cluster_E cluster_F \
            cluster_G cluster_H cluster_I any_ige single multiple \
            high_score medium_score low_score)
LEVELS=(family genus)
I=${SLURM_ARRAY_TASK_ID}
PART=${PARTITIONS[$(( I % 15 ))]}
LEVEL=${LEVELS[$(( I / 15 ))]}

echo "=== $(date) === task ${I}: ${LEVEL} / prev10 / ${PART} ==="
"${PY}" revision/analysis/04_predictive_ml.py \
  --level "${LEVEL}" --prevalence prev10 --partition "${PART}" \
  --models L1-logistic,L2-logistic,RandomForest,HistGradBoost --tag perm
echo "=== $(date) === Done ==="
