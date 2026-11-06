#!/bin/bash
#SBATCH --job-name=rev_package
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/package_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/package_%j.err

# Runs after the DA arrays finish: merge every method, compute the overlap, and
# build a self-contained archive to work on off-cluster.
#
# WHAT GOES IN THE ARCHIVE
#   All revision code, documentation, planning deliverables and AGGREGATE
#   result tables - everything needed to write the manuscript, the response
#   letter and the figures on a laptop.
#
# WHAT IS DELIBERATELY EXCLUDED, AND WHY
#   Individual-level KORA data (feature tables, per-subject metadata, the
#   matched-pair source tables). These are participant records under a data
#   agreement; a zip file in a downloads folder is the wrong place for them.
#   The aggregate results are sufficient for all downstream writing and
#   plotting. If a specific local task genuinely needs the source tables, add
#   them deliberately rather than shipping them by default.
#
#   Also excluded: the ~195 MB/partition permuted association matrices, which
#   are pure regenerable intermediates.
#
# Usage:
#   sbatch --dependency=afterany:<asv_array>:<ext_array> \
#          revision/slurm/99_finalize_and_package.sh

set -euo pipefail

REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
WORKSPACE=/home/itg/oleg.vlasovets/slr_example
STAMP=$(date '+%Y%m%d')
STAGE="/localscratch/${USER}/rev_package_${SLURM_JOB_ID}"
ARCHIVE_DIR="/ictstr01/scratch/users/${USER}/revision_archive"
NAME="ige_revision_${STAMP}"

cd "${REPO}"
mkdir -p "${STAGE}/${NAME}" "${ARCHIVE_DIR}"

set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env
set -u

echo "=== $(date) === [1/5] merge consensus DA grid ==="
python revision/analysis/02b_consensus_merge.py || echo "  (02b merge reported an issue - continuing)"

echo "=== $(date) === [2/5] six-method overlap ==="
python revision/analysis/02e_method_overlap.py | tee revision/results/tables/da_method_overlap_report.txt

echo "=== $(date) === [3/5] export Markdown briefing ==="
# Regenerated AFTER the overlap step so every number in docs/ reflects the
# final results rather than a partial grid.
python revision/analysis/10_export_markdown.py

echo "=== $(date) === [4/5] stage archive contents ==="
# Code, docs, config, slurm wrappers - everything but generated bulk.
rsync -a --exclude 'env/rlibs' --exclude 'results/networks/*/assoPerm*' \
      --exclude 'slurm/logs' --exclude '__pycache__' \
      revision/ "${STAGE}/${NAME}/revision/"

# Planning deliverables live in the workspace root, not the repo.
cp "${WORKSPACE}/REVISION_PLAN.md" "${WORKSPACE}/REVIEWER_RESPONSE_OUTLINE.md" \
   "${STAGE}/${NAME}/" 2>/dev/null || true

# Repo-level context a reader needs.
cp README.md requirements.txt "${STAGE}/${NAME}/" 2>/dev/null || true
rsync -a docs/ "${STAGE}/${NAME}/docs/" 2>/dev/null || true

# A short orientation note, so the archive is usable without this session.
cat > "${STAGE}/${NAME}/START_HERE.md" <<'NOTE'
# iScience revision — offline working copy

Built automatically at the end of the compute batch. Everything here is code,
documentation and **aggregate** results. No individual-level participant data is
included, by design — see `revision/slurm/99_finalize_and_package.sh`.

## Where things are

- `REVISION_PLAN.md` — the execution plan; §6 is the manuscript-correction inventory
- `REVIEWER_RESPONSE_OUTLINE.md` — response letter, keyed to R1/R2/R3
- `revision/TAXONOMY_COMPARISON.md` — the SILVA 138.2 reclassification analysis
- `revision/UPSTREAM_METHODS.md` — recovered DADA2 / taxonomy / phylogeny provenance
- `revision/results/tables/` — every aggregate result table
- `revision/analysis/` — all analysis modules
- `revision/slurm/` — the SLURM wrappers, each documenting its own cost

## The differential-abundance answer

Six methods were run on identical input with the same formula:

  LinDA, ANCOM-BC2, ALDEx2      (module 02  → `results/da/`)
  ZicoSeq, Wilcoxon-CLR, limma-voom (module 02d → `results/da_extended/`)

Merged and compared in `results/tables/`:

- `da_method_rates.csv` — per-method significance rate
- `da_method_overlap_pairwise.csv` — shared calls and Jaccard, every pair
- `da_method_upset.csv` — how many methods called each taxon (UpSet input)
- `da_method_effect_correlation.csv` — Spearman on effect sizes
- `da_method_overlap_report.txt` — the printed summary
- `da_all_methods_long.csv` — the merged long table behind all of the above

Two cautions carried forward in the code comments:

1. **Effect sizes are never averaged** across methods — the scales differ.
2. **ZicoSeq at ASV rank must be read via `q_native`**, not `q_within`: a
   permutation p has a floor of 1/(nperm+1), so BH over 1,459 features cannot
   reach 0.05 at nperm=4999. The scripts print this floor per analysis.

## What still needs the cluster

Nothing for writing. Regenerating results needs `slr-env` (R 4.4.3 / Bioc 3.20)
and the individual-level KORA tables, neither of which is in this archive.
NOTE

echo "=== $(date) === [5/5] compress ==="
cd "${STAGE}"
zip -qr "${ARCHIVE_DIR}/${NAME}.zip" "${NAME}"
cd /
rm -rf "${STAGE}"

echo
echo "=== $(date) === archive ready ==="
ls -lh "${ARCHIVE_DIR}/${NAME}.zip"
echo
echo "Download with:"
echo "  scp ${USER}@<cluster>:${ARCHIVE_DIR}/${NAME}.zip ."
