#!/bin/bash
#SBATCH --job-name=rev_alpha
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment='agent:claude-code; alpha diversity raw+rarefied, matched+adjusted, single job'
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/alpha_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/alpha_%j.err

# Alpha diversity (module 05c). SINGLE job, deliberately not an array: the
# expensive step (per-sample alpha over 100 rarefactions of a 508 x ~6,958
# matrix) is partition-invariant; the 15 partition tests are arithmetic.
#
# Sizing anchor: module 05b loaded the same feature table PLUS a 15,170-tip
# tree in 64G/8h with UniFrac on top. Here: no tree, no distance matrices.
# Table load ~3-5 min, 100 x (rrarefy + 3 indices) est. 10-40 min, tests
# seconds. 4 CPU / 32G / 2h is comfortable headroom, not a workload estimate
# ceiling. --nice=10000 is the HMGU baseline.
#
# Usage:  sbatch revision/slurm/05c_alpha_diversity.sh

set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"

set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env
set -u

echo "=== $(date) === alpha diversity ==="
Rscript revision/analysis/05c_alpha_diversity.R
echo "=== $(date) === Done ==="
