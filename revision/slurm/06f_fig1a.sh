#!/bin/bash
#SBATCH --job-name=rev_fig1a
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment='agent:claude-code; recover published Fig 1A order and rebuild Fig 1A'
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=0:20:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/fig1a_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/fig1a_%j.err
# 06f reads the 113-Mpx published PNG; 06g renders the rebuilt panel.
# STEP=recover|rebuild|all (default all)
set -euo pipefail
cd /home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
set +u; source ~/miniconda3/etc/profile.d/conda.sh; conda activate scont_plots; set -u
STEP=${STEP:-all}
[[ $STEP == recover || $STEP == all ]] && python3 revision/analysis/06f_fig1a_recover_order.py
[[ $STEP == rebuild || $STEP == all ]] && python3 revision/analysis/06g_fig1a_rebuild.py
exit 0
