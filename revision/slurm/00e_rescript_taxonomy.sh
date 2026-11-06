#!/bin/bash
#SBATCH --job-name=kora_rescript
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_max
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/KORA_data/preprocessing/logs/rescript_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/KORA_data/preprocessing/logs/rescript_%j.err

# A SEPARATE, fully documented taxonomy: RESCRIPt SILVA 138.2, curated and
# extracted to the V3-V4 amplicon, then used to reclassify the KORA ASVs.
#
# WHY
#   taxonomy_clean.csv came from a RESCRIPt classification, but no RESCRIPt call
#   exists in either preprocessing notebook, so its reference construction is
#   undocumented and unreproducible. Rather than guess at the original settings,
#   this builds a reference with the published best-practice parameters and runs
#   it as an independent analysis.
#
#   It does NOT overwrite taxonomy_clean.csv. Outputs go to a separate directory
#   and the two taxonomies are compared, so the original remains the record of
#   what the manuscript reports while this provides a documented reference point.
#
# WHAT "OPTIMAL" MEANS HERE
#   The one substantive choice is REGION EXTRACTION. A classifier trained on
#   full-length 16S loses precision on a short amplicon; extracting the reference
#   to the exact V3-V4 window bounded by 341F/785R is the documented way to
#   improve assignment, most visibly at species level. Since only 12.4% of KORA
#   ASVs currently carry a species label, this is directly testable: if extraction
#   materially raises that, it matters for the paper's species-level claims; if it
#   does not, the sparsity is inherent to the amplicon and the reference, which is
#   itself worth reporting.
#
#   Everything else follows the RESCRIPt SILVA protocol:
#     cull-seqs                    5 degenerate bases, 8-mer homopolymers
#     filter-seqs-length-by-taxon  Archaea 900, Bacteria 1200, Eukaryota 1400
#     dereplicate                  mode 'uniq' - retains distinct taxonomies for
#                                  identical sequences, which is the mode that
#                                  preserves species-level information
#
# VERSION
#   q2-2026.7-slr (the current QIIME 2), NOT the 2021.11 used for the reproduction.
#   This is deliberate and is why the run is framed as a separate analysis rather
#   than a reproduction: it changes two things at once relative to the original
#   taxonomy - the reference construction AND the software - so it cannot be read
#   as "what the original would have given with better settings". It is a
#   best-practice reference point, and the comparison against taxonomy_clean.csv
#   is reported with that caveat attached.
#
#   The 2021.11 rep-seqs artifact was verified to load in 2026.7 (FeatureData
#   [Sequence], UUID 7b9d4dae).
#
#   SILVA 138.2 is used - the newest curated release, and 2026.7's default. The
#   original analysis used SILVA 138.
#
# Usage:  sbatch run_rescript_taxonomy.sh

set -eo pipefail

KORA=/home/itg/oleg.vlasovets/slr_example/KORA_data
REPSEQS="${KORA}/q2_reprocessing/representative_sequences.qza"
OUT="${KORA}/q2_rescript_taxonomy"
SCRATCH="/localscratch/${USER}/kora_rescript_${SLURM_JOB_ID:-manual}"
NCPU="${SLURM_CPUS_PER_TASK:-8}"

# Reference-building intermediates are GB-scale and take ~40 min to regenerate.
# Checkpoint them on persistent scratch so a downstream failure never costs that
# work again - the first attempt died at extract-reads and threw away five
# completed stages.
CKPT="/ictstr01/scratch/users/${USER}/kora_rescript_ckpt"

mkdir -p "${OUT}" "${SCRATCH}" "${CKPT}" "${KORA}/preprocessing/logs"

# Restore any checkpointed artifacts into scratch, so completed stages are
# skipped. Explicit per-stage guards below are clearer than a generic wrapper,
# because several stages emit more than one artifact.
have() { for f in "$@"; do [ -f "${CKPT}/${f}" ] || return 1; done; return 0; }
restore() { for f in "$@"; do cp "${CKPT}/${f}" "${SCRATCH}/${f}"; done; }
keep()    { for f in "$@"; do cp "${SCRATCH}/${f}" "${CKPT}/${f}"; done; }
export TMPDIR="${SCRATCH}/tmp"; mkdir -p "${TMPDIR}"

echo "=== $(date) === RESCRIPt SILVA 138 V3-V4 taxonomy ==="
echo "  rep-seqs : ${REPSEQS}"
echo "  outputs  : ${OUT}"
echo "  scratch  : ${SCRATCH}"
echo "  cpus     : ${NCPU}"

[ -f "${REPSEQS}" ] || { echo "FATAL: rep-seqs missing; run the DADA2 rebuild first" >&2; exit 1; }

# q2-2026.7-slr is not a relocated env, so the Rscript prefix problem that
# affected qiime2-2021.11 does not apply here.
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate q2-2026.7-slr
set -u
qiime --version

# get-silva-data reaches out to the SILVA servers. Fail fast if the node is
# offline rather than after an hour of setup.
echo "=== $(date) === preflight: connectivity ==="
curl -sSf -I --max-time 60 https://www.arb-silva.de >/dev/null \
  && echo "  SILVA reachable" \
  || { echo "FATAL: cannot reach SILVA from this node" >&2; exit 1; }

cd "${SCRATCH}"

# ── 1. Fetch SILVA 138 SSU NR99 with species labels ───────────────────────
echo "=== $(date) === [1/7] get-silva-data ==="
if have silva-rna-seqs.qza silva-tax.qza; then
  echo "  restoring checkpoint: silva-rna-seqs.qza silva-tax.qza"; restore silva-rna-seqs.qza silva-tax.qza
else
  qiime rescript get-silva-data \
    --p-version 138.2 --p-target SSURef_NR99 \
    --p-include-species-labels \
    --o-silva-sequences silva-rna-seqs.qza \
    --o-silva-taxonomy silva-tax.qza \
    --verbose
  cp silva-tax.qza "${OUT}/"
  keep silva-rna-seqs.qza silva-tax.qza
fi

# ── 2. RNA -> DNA ─────────────────────────────────────────────────────────
echo "=== $(date) === [2/7] reverse-transcribe ==="
if have silva-seqs.qza; then
  echo "  restoring checkpoint: silva-seqs.qza"; restore silva-seqs.qza
else
  qiime rescript reverse-transcribe \
    --i-rna-sequences silva-rna-seqs.qza \
    --o-dna-sequences silva-seqs.qza
  keep silva-seqs.qza
fi

# ── 3. Cull low-quality reference sequences ───────────────────────────────
echo "=== $(date) === [3/7] cull-seqs ==="
if have silva-seqs-culled.qza; then
  echo "  restoring checkpoint: silva-seqs-culled.qza"; restore silva-seqs-culled.qza
else
  qiime rescript cull-seqs \
    --i-sequences silva-seqs.qza \
    --p-num-degenerates 5 --p-homopolymer-length 8 \
    --p-n-jobs "${NCPU}" \
    --o-clean-sequences silva-seqs-culled.qza
  keep silva-seqs-culled.qza
fi

# ── 4. Length filter, per taxonomic group ─────────────────────────────────
echo "=== $(date) === [4/7] filter-seqs-length-by-taxon ==="
if have silva-seqs-filt.qza silva-seqs-discard.qza; then
  echo "  restoring checkpoint: silva-seqs-filt.qza silva-seqs-discard.qza"; restore silva-seqs-filt.qza silva-seqs-discard.qza
else
  qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-seqs-culled.qza \
    --i-taxonomy silva-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs silva-seqs-filt.qza \
    --o-discarded-seqs silva-seqs-discard.qza
  keep silva-seqs-filt.qza silva-seqs-discard.qza
fi

# ── 5. Dereplicate the full-length reference ──────────────────────────────
echo "=== $(date) === [5/7] dereplicate (full length) ==="
if have silva-derep-seqs.qza silva-derep-tax.qza; then
  echo "  restoring checkpoint: silva-derep-seqs.qza silva-derep-tax.qza"; restore silva-derep-seqs.qza silva-derep-tax.qza
else
  qiime rescript dereplicate \
    --i-sequences silva-seqs-filt.qza \
    --i-taxa silva-tax.qza \
    --p-mode uniq --p-threads "${NCPU}" \
    --o-dereplicated-sequences silva-derep-seqs.qza \
    --o-dereplicated-taxa silva-derep-tax.qza
  keep silva-derep-seqs.qza silva-derep-tax.qza
fi

# ── 6. Extract the V3-V4 amplicon, then dereplicate again ─────────────────
# This is the step that distinguishes this reference from a full-length one.
echo "=== $(date) === [6/7] extract-reads (341F / 785R) ==="
# NOTE: QIIME 2 2026.7 added a REQUIRED --o-read-extraction-stats output that
# 2021.11 did not have. Omitting it is what killed the first attempt. The stats
# artifact is genuinely useful anyway - per-sequence primer match outcome - so it
# is kept and exported for the Methods.
qiime feature-classifier extract-reads \
  --i-sequences silva-derep-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-n-jobs "${NCPU}" \
  --p-read-orientation forward \
  --o-reads silva-v3v4.qza \
  --o-read-extraction-stats silva-v3v4-extraction-stats.qza \
  --verbose
cp silva-v3v4-extraction-stats.qza "${OUT}/" 2>/dev/null || true

qiime rescript dereplicate \
  --i-sequences silva-v3v4.qza \
  --i-taxa silva-derep-tax.qza \
  --p-mode uniq --p-threads "${NCPU}" \
  --o-dereplicated-sequences silva-v3v4-derep-seqs.qza \
  --o-dereplicated-taxa silva-v3v4-derep-tax.qza
cp silva-v3v4-derep-seqs.qza silva-v3v4-derep-tax.qza "${OUT}/"

# ── 7. Train and classify ─────────────────────────────────────────────────
echo "=== $(date) === [7/7] fit-classifier-naive-bayes ==="
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-v3v4-derep-seqs.qza \
  --i-reference-taxonomy silva-v3v4-derep-tax.qza \
  --o-classifier silva-138.2-v3v4-nb-classifier.qza \
  --verbose
cp silva-138.2-v3v4-nb-classifier.qza "${OUT}/"

echo "=== $(date) === classify-sklearn on the KORA ASVs ==="
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138.2-v3v4-nb-classifier.qza \
  --i-reads "${REPSEQS}" \
  --p-n-jobs "${NCPU}" \
  --o-classification kora-taxonomy-rescript-v3v4.qza \
  --verbose
cp kora-taxonomy-rescript-v3v4.qza "${OUT}/"

qiime tools export --input-path kora-taxonomy-rescript-v3v4.qza \
  --output-path export_tax
cp export_tax/taxonomy.tsv "${OUT}/kora_taxonomy_rescript_v3v4.tsv"

echo "=== $(date) === cleaning scratch ==="
rm -rf "${SCRATCH}"
echo "=== $(date) === Done. Outputs in ${OUT} ==="
ls -la "${OUT}"
