#!/bin/bash
#SBATCH --job-name=rev_ml_primary
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment='agent:claude-code; ML benchmark at the primary contrast rank (species/prev05) + ASV'
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --array=0-15
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/ml_primary_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/ml_primary_%A_%a.err

# ML benchmark where the paper's primary contrast lives. The existing 30 runs
# cover family/genus at prev10 only; the pre-specified primary contrast is
# SPECIES at prev05, and the retracted 61-ASV claim lived at ASV rank. Tasks
# 0-14 = species/prev05 x 15 partitions; task 15 = asv/prev05/any_ige.
# Full model set (incl. XGBoost from revision/env/xgb-venv) with a per-model
# within-pair permutation null - a null computed for one model does not
# license a p-value for another.
#
# Sizing: rev_ml (Aug 5) 30 tasks ~37 min total wall at genus/112 features
# with an L2-only null. Species=60 features but ~6 model-nulls (RF's 1000
# trees x 200 perms x 5 folds dominates): est 30-90 min/task; ASV/prev05
# (489 features) up to ~4h. 8 CPU/16G/6h = headroom, not an estimate.
set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"
PY=${REPO}/revision/env/xgb-venv/bin/python

PARTITIONS=(cluster_A cluster_B cluster_C cluster_D cluster_E cluster_F \
            cluster_G cluster_H cluster_I any_ige single multiple \
            high_score medium_score low_score)
I=${SLURM_ARRAY_TASK_ID}
if [ "${I}" -le 14 ]; then
  LEVEL=species; PREV=prev05; PART=${PARTITIONS[$I]}
else
  LEVEL=asv; PREV=prev05; PART=any_ige
fi

echo "=== $(date) === task ${I}: ${LEVEL} / ${PREV} / ${PART} ==="
"${PY}" revision/analysis/04_predictive_ml.py \
  --level "${LEVEL}" --prevalence "${PREV}" --partition "${PART}" \
  --models Dummy,L1-logistic,L2-logistic,RandomForest,HistGradBoost,XGBoost
echo "=== $(date) === Done ==="
