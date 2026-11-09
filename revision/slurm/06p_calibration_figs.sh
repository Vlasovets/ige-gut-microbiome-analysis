#!/bin/bash
#SBATCH --job-name=rev_calib_figs
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment="agent:claude-code; p-value histograms + re-render Fig 4A and DA heatmap"
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/calibfigs_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/calibfigs_%j.err
set -euo pipefail
cd /home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/analysis
set +u; source ~/miniconda3/etc/profile.d/conda.sh; conda activate scont_plots; set -u
python3 06p_pvalue_histograms.py
python3 06o_fig4A_alpha.py > /dev/null
python3 06h_method_heatmap_nearing.py
