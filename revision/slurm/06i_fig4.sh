#!/bin/bash
#SBATCH --job-name=rev_fig4
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment='agent:claude-code; render revised Figure 4B-C'
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0:15:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/fig4_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/fig4_%j.err
# Plot rendering; runs here because the login node killed it (exit 137).
set -euo pipefail
cd /home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
set +u; source ~/miniconda3/etc/profile.d/conda.sh; conda activate scont_plots; set -u
python3 revision/analysis/06i_fig4_da.py
