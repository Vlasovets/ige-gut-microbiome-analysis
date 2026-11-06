#!/bin/bash
#SBATCH --job-name=rev_ml_xgb
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment='agent:claude-code; XGBoost backfill onto the existing family/genus prev10 grid'
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --array=0-29
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/ml_xgb_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/ml_xgb_%A_%a.err

# Completes the model-comparison table: XGBoost on the 30 cells the August run
# covered (family+genus x prev10 x 15 partitions). --tag xgb so the CSVs
# behind already-drafted numbers are never rewritten; Dummy rides along as the
# in-run floor. Est 5-10 min/task (anchor: 37 min for the whole August array).
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
  --models Dummy,XGBoost --tag xgb
echo "=== $(date) === Done ==="
