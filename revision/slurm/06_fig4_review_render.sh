#!/bin/bash
#SBATCH --job-name=rev_fig4_render
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment="agent:claude-code; render Figure 4 panels + DA heatmap after co-author review"
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/fig4render_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/fig4render_%j.err
# Renders Figure 4A (06o), 4B (06n), 4C volcanoes (06k, LinDA + ANCOM-BC2, score set)
# and the appendix DA heatmap (06h). Plot rendering is killed on the login node.
set -euo pipefail
cd /home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/analysis
set +u; source ~/miniconda3/etc/profile.d/conda.sh; conda activate scont_plots; set -u
python3 06o_fig4A_alpha.py
python3 06n_fig4B_family_composition.py
python3 06k_fig4_volcanoes_contrasts.py --set score
python3 06k_fig4_volcanoes_contrasts.py --set score --method ANCOMBC2
python3 06h_method_heatmap_nearing.py
