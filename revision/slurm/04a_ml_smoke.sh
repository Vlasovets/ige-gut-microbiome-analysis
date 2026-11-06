#!/bin/bash
#SBATCH --job-name=rev_ml_smoke
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment='agent:claude-code; ML plumbing smoke test gating 04b/04c arrays'
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/ml_smoke_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/ml_smoke_%j.err

# Plumbing check for the patched 04_predictive_ml.py before 46 array tasks:
# XGBoost imports from xgb-venv, per-model nulls compute, --tag suffixes the
# outputs, originals untouched. Tiny parameters on purpose - this validates
# code paths, not statistics. Runs on a compute node because a login-node
# attempt was OOM-killed (exit 137).
set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"
PY=${REPO}/revision/env/xgb-venv/bin/python

echo "=== $(date) === smoke: genus/prev10/any_ige, Dummy+XGBoost, 1 rep, 10 perms ==="
"${PY}" revision/analysis/04_predictive_ml.py \
  --level genus --prevalence prev10 --partition any_ige \
  --models Dummy,XGBoost --tag smoketest --repeats 1 --n-perm 10

echo "=== checks ==="
ls -1 revision/results/ml/ | grep smoketest
# the original (untagged) metrics file must be byte-identical to git HEAD
git diff --quiet -- revision/results/ml/ml_genus_prev10_any_ige_metrics.csv \
  && echo "original metrics untouched: PASS" \
  || { echo "original metrics CHANGED: FAIL"; exit 1; }
echo "=== $(date) === smoke done ==="
