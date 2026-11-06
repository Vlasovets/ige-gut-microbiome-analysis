#!/bin/bash
#SBATCH --job-name=rev_nested
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment='agent:claude-code; M0-M3 nested-model chain explaining 61 unadjusted ASV hits -> matched+adjusted'
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=0-102
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/nested_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/nested_%A_%a.err

# Array map: tasks 0-99 = M0 seeds 42..141 (14 LinDA runs each - any_ige is
# seed-invariant and lives in M1); task 100 = M1 (15 runs); 101 = M2; 102 = M3.
# Sizing: LinDA at ASV rank on 4-7k features x <=500 samples, permutation-free
# linear modelling; est 1-4 min per partition run -> 15-60 min/task. 4h cap is
# headroom (a task nearing it means something is wrong, not slow).
set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env
set -u

I=${SLURM_ARRAY_TASK_ID}
if [ "${I}" -le 99 ]; then
  SEED=$(( 42 + I ))
  echo "=== $(date) === task ${I}: M0 seed ${SEED} ==="
  Rscript revision/analysis/07b_nested_linda.R --model M0 --seed "${SEED}"
else
  MODELS=(M1 M2 M3)
  M=${MODELS[$(( I - 100 ))]}
  echo "=== $(date) === task ${I}: ${M} ==="
  Rscript revision/analysis/07b_nested_linda.R --model "${M}"
fi
echo "=== $(date) === Done ==="
