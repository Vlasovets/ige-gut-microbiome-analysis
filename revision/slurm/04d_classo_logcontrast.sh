#!/bin/bash
#SBATCH --job-name=rev_classo
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment="agent:claude-code; c-lasso log-contrast models, same folds as RF/XGB"
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/classo_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/classo_%A_%a.err

# Sparse log-contrast classification (c-lasso C1, the main-paper method) on the
# SAME 46 analyses, matched pairs and pair-grouped folds as the RF/XGBoost
# benchmark (module 04). c-lasso from revision/env/classo-venv (pinned commit
# with the 2026-11-04 classification-path fixes). Writes ml_*_classo_* only.
#
# Tasks: 0-14 family/prev10, 15-29 genus/prev10, 30-44 species/prev05 (15 IgE
# groups each), 45 asv/prev05/any_ige - the existing RF/XGB grid.
# N_PERM_CLASSO: within-pair permutations; set from the smoke-test timing.
#
# Usage:
#   sbatch --array=39 --export=ALL,N_PERM_CLASSO=2 revision/slurm/04d_classo_logcontrast.sh   # smoke: species/any_ige
#   sbatch --array=0-45 revision/slurm/04d_classo_logcontrast.sh
set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"
PY=${REPO}/revision/env/classo-venv/bin/python
N_PERM=${N_PERM_CLASSO:-100}

PARTITIONS=(cluster_A cluster_B cluster_C cluster_D cluster_E cluster_F \
            cluster_G cluster_H cluster_I any_ige single multiple \
            high_score medium_score low_score)
I=${SLURM_ARRAY_TASK_ID}
if   (( I < 15 )); then LEVEL=family;  PREV=prev10; PART=${PARTITIONS[$I]}
elif (( I < 30 )); then LEVEL=genus;   PREV=prev10; PART=${PARTITIONS[$(( I - 15 ))]}
elif (( I < 45 )); then LEVEL=species; PREV=prev05; PART=${PARTITIONS[$(( I - 30 ))]}
else                    LEVEL=asv;     PREV=prev05; PART=any_ige
fi

F=revision/results/ml/ml_${LEVEL}_${PREV}_${PART}_classo_metrics.csv
echo "=== $(date) === task ${I}: ${LEVEL} / ${PREV} / ${PART} (n_perm=${N_PERM}) ==="
[[ -s "$F" && -z "${N_PERM_CLASSO:-}" ]] && { echo "exists, skipping"; exit 0; }
"${PY}" revision/analysis/04_predictive_ml.py \
  --level "${LEVEL}" --prevalence "${PREV}" --partition "${PART}" \
  --classo --n-perm-classo "${N_PERM}"
echo "=== $(date) === Done ==="
