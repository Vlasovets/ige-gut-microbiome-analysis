#!/bin/bash
#SBATCH --job-name=rev_net
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --array=0-4
#SBATCH --output=revision/slurm/logs/net_%A_%a.out
#SBATCH --error=revision/slurm/logs/net_%A_%a.err

# Genus-level networks with permutation-tested differential edges.
#
# Scope: only partitions where graphical-lasso edge selection is identifiable.
# Edge count is monotone in n at fixed p (measured: 6 edges at n=23 ... 66 at
# n=221), so cluster networks at n=22-24 would measure sample size, not biology.
# Partitions with <15 complete pairs skip themselves.
#
# nPerm = 5000, not 1000. A permutation p-value has a floor of 1/(nPerm+1), so
# with m testable edges BH can only reach q<0.05 if nPerm > m/0.05 - 1. At
# m = 89 testable edges that needs nPerm > 1779; 1000 would leave the test
# structurally unable to return a positive result. 5000 gives headroom.
#
# Timing: ~2.7 s/permutation with stored permuted associations reused between
# netCompare and diffnet, so 5000 perms is ~4 h. 10 h gives margin.

set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"

PARTITIONS=(any_ige multiple medium_score low_score cluster_I)
PART=${PARTITIONS[${SLURM_ARRAY_TASK_ID}]}
# matched (default) or full: sbatch --export=ALL,DESIGN=full ...
DESIGN=${DESIGN:-matched}
NPERM=${NPERM:-5000}

source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env

echo "=== $(date) === networks: ${PART} (${DESIGN}, nperm=${NPERM}) ==="
Rscript revision/analysis/03_genus_networks.R --partition "${PART}" --nperm "${NPERM}" --design "${DESIGN}"
echo "=== $(date) === Done ==="
