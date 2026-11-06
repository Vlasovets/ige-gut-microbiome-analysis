#!/bin/bash
#SBATCH --job-name=rev_da_asv_smoke
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/da_asv_smoke_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/da_asv_smoke_%j.err

# One task, the HEAVIEST cell in the ASV grid: prev01 / any_ige = 1459 features.
# Everything else in the array is cheaper, so if this fits, the array fits.
#
# Purpose is sizing and config validation, not a result - though the result it
# produces is real and the array will simply overwrite it.
#
# Read the true cost afterwards with:
#   sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,State

set -euo pipefail

REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"

set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env
set -u

echo "=== $(date) === SMOKE: asv / prev01 / any_ige (1459 features) ==="
Rscript revision/analysis/02_consensus_da.R \
  --level asv --prevalence prev01 --partition any_ige \
  --mc-samples 128
echo "=== $(date) === Smoke test done ==="
