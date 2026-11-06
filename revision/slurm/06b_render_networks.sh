#!/bin/bash
#SBATCH --job-name=rev_net_render
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment='agent:claude-code; render network PDFs for offline review (per-sample counts stay on cluster)'
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/net_render_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/net_render_%j.err

# Renders the 3 network views per partition to results/figures/interactive/
# from the cached netprops RDS (seeded, verified edge-identical). Exists
# because the interactive script needs per-sample counts, which must not leave
# the cluster, so the figures are produced here for offline review. Caches
# make each partition a few seconds; 30 min is generous headroom.
set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"
set +u; source ~/miniconda3/etc/profile.d/conda.sh; conda activate slr-env; set -u
for P in any_ige multiple medium_score low_score cluster_I; do
  echo "=== $(date) === ${P} ==="
  NETCOMI_PARTITION=${P} NETCOMI_SAVE_PLOTS=1 Rscript revision/analysis/06b_networks_interactive.R
done
ls -la revision/results/figures/interactive/
