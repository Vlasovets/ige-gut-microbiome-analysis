#!/bin/bash
#SBATCH --job-name=rev_pack_offline
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --comment='agent:claude-code; build offline review zip (aggregate results only)'
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=0:45:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/pack_offline_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/pack_offline_%j.err
# Lustre-bound copy + zip of ~150 MB; took >2 min on the login node, so it
# runs here. set -x in the wrapper call logs which step is slow.
set -euo pipefail
bash -x /home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/tools/pack_offline.sh
