#!/bin/bash
#SBATCH --job-name=rev_ml
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --array=0-29
#SBATCH --output=revision/slurm/logs/ml_%A_%a.out
#SBATCH --error=revision/slurm/logs/ml_%A_%a.err

# Leakage-free predictive benchmark across the IgE subgroups.
# 2 levels (family, genus) x 15 partitions = 30 tasks.
# Partitions with fewer than 10 complete pairs skip themselves.
#
# Timing: ~10 min at 4 repeats / 60 permutations on genus (112 features);
# 10 repeats / 200 permutations is ~5x that, so 3 h gives headroom.

set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"

LEVELS=(family genus)
PARTITIONS=(cluster_A cluster_B cluster_C cluster_D cluster_E cluster_F \
            cluster_G cluster_H cluster_I any_ige single multiple \
            high_score medium_score low_score)
N=${#PARTITIONS[@]}
I=${SLURM_ARRAY_TASK_ID}
PART=${PARTITIONS[$(( I % N ))]}
LEVEL=${LEVELS[$(( I / N ))]}

source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env
export PYTHONWARNINGS=ignore

echo "=== $(date) === task ${I}: ${LEVEL} / prev10 / ${PART} ==="
python -u revision/analysis/04_predictive_ml.py \
  --level "${LEVEL}" --prevalence prev10 --partition "${PART}" \
  --repeats 10 --n-perm 200
echo "=== $(date) === Done ==="
