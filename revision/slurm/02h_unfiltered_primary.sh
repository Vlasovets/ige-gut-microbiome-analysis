#!/bin/bash
#SBATCH --job-name=rev_da_unf
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment="agent:claude-code; unfiltered full-cohort DA, primary contrast"
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/daunf_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/daunf_%A_%a.err

# Unfiltered ("prev00") feature set of the submitted manuscript, run through the
# same 5 covariate-adjusted methods as the filtered grid, primary contrast only
# (all 275 sensitized vs 233 controls, ~ IgE status + sex + age + BMI).
# Answers "what does the unfiltered table give under the revised methods".
#
#   inputs   build prev00 inputs (02f --unfiltered), asserts 6,958 ASVs
#   da       array 0-3 = family, genus, species, asv
#
# Cost: ASV = 6,958 features vs 1,459 at prev01 (<=12 min/script measured), so
# expect ~1-4 h; ZicoSeq dominates. At m = 6,958 and nperm = 4999 ZicoSeq's
# smallest attainable BH q is 1.39, so its q_within CANNOT reach 0.05 at ASV
# rank - a structural zero, reported as such.
#
# Usage:
#   jid=$(sbatch --parsable revision/slurm/02h_unfiltered_primary.sh inputs)
#   sbatch --array=0-3 --dependency=afterok:$jid revision/slurm/02h_unfiltered_primary.sh da

# Rscript is NOT used: slr-env was relocated and its compiled Rscript still
# points at the old prefix ("Rscript execution error"). The R launcher is a
# script, was prefix-fixed, and gives the same commandArgs (--file=, --args).
set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"
MODE=${1:?mode required: inputs | da}

set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env
set -u

if [[ "${MODE}" == "inputs" ]]; then
  python revision/analysis/02f_full_cohort_inputs.py --unfiltered
  exit 0
fi

LEVELS=(family genus species asv)
LEVEL=${LEVELS[${SLURM_ARRAY_TASK_ID}]}
PREV=prev00; PART=any_ige
F1=revision/results/da_full/da_long_${LEVEL}_${PREV}_${PART}.csv
F2=revision/results/da_extended_full/daext_long_${LEVEL}_${PREV}_${PART}.csv

echo "=== $(date) === ${LEVEL} / ${PREV} / ${PART} [full] ==="
[[ -s "${F1}" ]] || R --no-echo --no-restore --file=revision/analysis/02_consensus_da.R --args --design full \
  --level "${LEVEL}" --prevalence "${PREV}" --partition "${PART}" --mc-samples 128
[[ -s "${F2}" ]] || R --no-echo --no-restore --file=revision/analysis/02d_extended_da.R --args --design full \
  --level "${LEVEL}" --prevalence "${PREV}" --partition "${PART}" --nperm 4999
echo "=== $(date) === Done ==="
