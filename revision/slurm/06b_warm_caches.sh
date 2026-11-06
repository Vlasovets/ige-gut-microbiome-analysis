#!/bin/bash
#SBATCH --job-name=rev_net_warm
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment='agent:claude-code; headless smoke test + cache warm-up of the interactive NetCoMi script'
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/net_warm_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/net_warm_%j.err

# Runs 06b_networks_interactive.R headless once per partition, with plots to a
# scratch PDF device. Purpose: (1) prove the script the user will open in
# RStudio actually runs end to end; (2) pre-build the netprops RDS caches so
# their first interactive session is instant; (3) emit the overlay annotation
# and supplementary edge CSVs. Seeded StARS -> the QC line must print MATCH
# against the stored edge tables for every partition.
set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env
set -u
for P in any_ige multiple medium_score low_score cluster_I; do
  echo "=== $(date) === ${P} ==="
  NETCOMI_PARTITION=${P} Rscript revision/analysis/06b_networks_interactive.R
done
echo "=== $(date) === all caches warm ==="
ls -la revision/results/networks/*/netprops_*.rds
