#!/bin/bash
#SBATCH --job-name=rev_phylo
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --array=0-14
#SBATCH --output=revision/slurm/logs/phylo_%A_%a.out
#SBATCH --error=revision/slurm/logs/phylo_%A_%a.err

# Faith's PD and UniFrac across the IgE subgroups.
# Unblocked by the phylogeny rebuild: the tree reproduces the published ASV set
# exactly (15,170/15,170 tips), so these metrics apply to the analyses the
# manuscript reports.
#
# Memory is set high because GUniFrac accumulates abundance along every branch
# of a 15,170-tip tree pruned per partition.

set -euo pipefail
REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
cd "${REPO}"

PARTITIONS=(cluster_A cluster_B cluster_C cluster_D cluster_E cluster_F \
            cluster_G cluster_H cluster_I any_ige single multiple \
            high_score medium_score low_score)
PART=${PARTITIONS[${SLURM_ARRAY_TASK_ID}]}

set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env
set -u

echo "=== $(date) === phylo diversity: ${PART} ==="
Rscript revision/analysis/05b_phylo_diversity.R --partition "${PART}" --nperm 5000
echo "=== $(date) === Done ==="
