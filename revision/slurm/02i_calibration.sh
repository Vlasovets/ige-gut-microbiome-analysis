#!/bin/bash
#SBATCH --job-name=rev_calib
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment="agent:claude-code; DA calibration check (permuted labels)"
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/calib_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/calib_%A_%a.err
# Calibration check (module 02i): tasks 0-20 = family, perm 0..20;
# tasks 21-41 = asv, perm 0..20. perm 0 = real labels (reproducibility check).
# Usage: sbatch --array=0-41 revision/slurm/02i_calibration.sh
set -euo pipefail
cd /home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
set +u; source ~/miniconda3/etc/profile.d/conda.sh; conda activate slr-env; set -u
I=${SLURM_ARRAY_TASK_ID}
LEVEL=$([[ $I -lt 21 ]] && echo family || echo asv)
PERM=$(( I % 21 ))
F=revision/results/calibration/calib_${LEVEL}_perm$(printf %02d $PERM).csv
[[ -s "$F" ]] || Rscript revision/analysis/02i_calibration_check.R --level "$LEVEL" --perm "$PERM"
