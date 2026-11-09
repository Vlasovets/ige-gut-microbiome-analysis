#!/bin/bash
#SBATCH --job-name=rev_abl
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment="agent:claude-code; ANCOM-BC2 ablation on permuted labels"
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=2:00:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/abl_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/abl_%A_%a.err
# Module 02j: 4 configs x perms 0..10 = 44 tasks (task = config_index * 11 + perm).
# Usage: sbatch --array=0-43 revision/slurm/02j_ancombc2_ablation.sh
set -euo pipefail
cd /home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
set +u; source ~/miniconda3/etc/profile.d/conda.sh; conda activate slr-env; set -u
CONFIGS=(base pseudo1 nostruc nocov)
I=${SLURM_ARRAY_TASK_ID}; CONFIG=${CONFIGS[$(( I / 11 ))]}; PERM=$(( I % 11 ))
F=revision/results/calibration/ablation_${CONFIG}_perm$(printf %02d $PERM).csv
[[ -s "$F" ]] || Rscript revision/analysis/02j_ancombc2_ablation.R --config "$CONFIG" --perm "$PERM"
