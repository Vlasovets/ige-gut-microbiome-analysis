#!/bin/bash
#SBATCH --job-name=rev_da_full
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment="agent:claude-code; full-cohort covariate-adjusted DA"
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/dafull_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/dafull_%A_%a.err

# Full-cohort, covariate-adjusted DA: the adjustment-only alternative to the
# matched design. Same 5 covariate-adjusted methods, same formula
# (~ ige_status + sex + age + bmi), same feature sets; every sensitized member
# of a partition vs all 233 controls. Paired Wilcoxon has no pairs here and is
# skipped (see 02d).
#
# Two modes:
#   inputs   build per-analysis tables (02f_full_cohort_inputs.py), one task
#   da       one array task = module 02 + module 02d for one cell of
#            4 levels x 3 prevalence x 15 partitions = 180 tasks
#
# MEASURED on the matched grid: <=12 min per script per task, peak RSS 10.2 GB.
#
# Usage:
#   jid=$(sbatch --parsable revision/slurm/02f_full_cohort_da.sh inputs)
#   sbatch --array=9 --dependency=afterok:$jid revision/slurm/02f_full_cohort_da.sh da   # smoke: family/prev10/any_ige
#   sbatch --array=0-179 revision/slurm/02f_full_cohort_da.sh da

set -euo pipefail

REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"
MODE=${1:?mode required: inputs | da}

set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env
set -u

if [[ "${MODE}" == "inputs" ]]; then
  echo "=== $(date) === building full-cohort inputs ==="
  python revision/analysis/02f_full_cohort_inputs.py
  echo "=== $(date) === Done ==="
  exit 0
fi

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

# Retry-safe: a cell whose two outputs both exist is not recomputed.
F1=revision/results/da_full/da_long_${LEVEL}_${PREV}_${PART}.csv
F2=revision/results/da_extended_full/daext_long_${LEVEL}_${PREV}_${PART}.csv

echo "=== $(date) === task ${I}: ${LEVEL} / ${PREV} / ${PART} [full] ==="
if [[ ! -s "${F1}" ]]; then
  Rscript revision/analysis/02_consensus_da.R --design full \
    --level "${LEVEL}" --prevalence "${PREV}" --partition "${PART}" --mc-samples 128
fi
if [[ ! -s "${F2}" ]]; then
  Rscript revision/analysis/02d_extended_da.R --design full \
    --level "${LEVEL}" --prevalence "${PREV}" --partition "${PART}" --nperm 4999
fi
echo "=== $(date) === Done ==="
