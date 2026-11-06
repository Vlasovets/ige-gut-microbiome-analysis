#!/bin/bash
#SBATCH --job-name=kora_dada2
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_max
#SBATCH --cpus-per-task=32
#SBATCH --mem=480G
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/KORA_data/preprocessing/logs/dada2_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/KORA_data/preprocessing/logs/dada2_%j.err

# Reconstruct the KORA 16S phylogeny (and the DADA2 intermediates) so that
# UniFrac and Faith's PD become computable for the iScience revision.
#
# WHY THIS IS NEEDED
#   The published analyses used data/q2/tree/rooted_tree_15170/tree.nwk, which is
#   not on disk. ASV identifiers are MD5 hashes of the sequences, so the
#   sequences cannot be recovered from feature_table.tsv - the tree has to be
#   rebuilt from the reads.
#
# PARAMETERS
#   Taken verbatim from q2_dada2_KORA.ipynb (cell 8), which is the record of what
#   produced the published feature table:
#     --p-trim-left-f 17   (length of 341F, CCTACGGGNGGCWGCAG)
#     --p-trim-left-r 21   (length of 785R, GACTACHVGGGTATCTAATCC)
#     --p-trunc-len-f 245  --p-trunc-len-r 245   (reads are 250 bp; verified)
#     --p-max-ee-f 2       --p-max-ee-r 2
#     --p-trunc-q 2
#   Chimera removal is left at the default, 'consensus'.
#
# VERSION
#   qiime2-2021.11 deliberately, NOT the newer 2026.7. The manuscript states
#   QIIME 2 v2021.11, and matching it is what gives the re-run a chance of
#   reproducing the same 15,170 ASVs. Because an ASV id is the MD5 of its
#   sequence, identical sequences yield identical ids - so if the ASVs reproduce,
#   the new tree's tip labels line up with the existing feature table and every
#   downstream analysis stays valid. A different QIIME 2 version would risk a
#   tree whose tips do not match the data the manuscript reports.
#
# RESOURCES
#   cpu_max (3-day ceiling, 32-core cap) rather than cpu_normal (24 h, 512
#   cores). DADA2 cannot checkpoint, so a wall-clock kill discards the entire
#   run; the extra cores are not worth that risk. Sized for completion, not
#   speed. cpu_long is not available to this account.

set -euo pipefail

KORA=/home/itg/oleg.vlasovets/slr_example/KORA_data
SRC=/lustre/groups/bds01/datasets/KORA/paired-end-demux.qza
OUT="${KORA}/q2_reprocessing"
SCRATCH="/localscratch/${USER}/kora_dada2_${SLURM_JOB_ID:-manual}"

mkdir -p "${OUT}" "${SCRATCH}" "${KORA}/preprocessing/logs"

# QIIME 2 unpacks artifacts into TMPDIR. Keep that on node-local scratch, never
# on Lustre: the extracted reads are tens of GB and the I/O pattern is hostile
# to a shared filesystem.
export TMPDIR="${SCRATCH}/tmp"
mkdir -p "${TMPDIR}"

echo "=== $(date) === KORA DADA2 + phylogeny ==="
echo "  source   : ${SRC}"
echo "  scratch  : ${SCRATCH}  ($(df -h /localscratch | tail -1 | awk '{print $4}') free)"
echo "  outputs  : ${OUT}"
echo "  cpus     : ${SLURM_CPUS_PER_TASK}   mem: ${SLURM_MEM_PER_NODE:-?}MB"

# The qiime2-2021.11 activation hook references ADDR2LINE unguarded, which trips
# `set -u`. Disable nounset just across activation, then restore it.
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2021.11
set -u
qiime --version

# ── Preflight: R must actually work before we stage 8 GB and burn hours ────
#
# The qiime2-2021.11 env was relocated to /lustre/groups/shared/apps/conda. The
# `R` launcher is a shell script and was prefix-fixed during that move, but
# `Rscript` is a COMPILED BINARY with the original prefix baked in:
#     /home/icb/oleg.vlasovetc/miniconda3/envs/qiime2-2021.11/lib/R
# so DADA2 died with "Rscript execution error: No such file or directory" after
# five minutes. Exporting R_HOME does not help - Rscript fails before reading it.
#
# Fix: the legacy prefix is repopulated with symlinks to the relocated env
# (see LEGACY_PREFIX below), which satisfies the baked-in path. Verify it here
# rather than discovering the breakage again mid-run.
LEGACY_PREFIX=/home/icb/oleg.vlasovetc/miniconda3/envs/qiime2-2021.11
if [ ! -e "${LEGACY_PREFIX}/lib/R/bin/exec/R" ]; then
  echo "FATAL: legacy R prefix missing at ${LEGACY_PREFIX}" >&2
  echo "  Recreate it with symlinks to the relocated env, then resubmit." >&2
  exit 1
fi
echo "=== preflight: R and dada2 ==="
Rscript -e 'suppressMessages(library(dada2));
  cat(sprintf("  %s | dada2 %s\n", R.version.string, packageVersion("dada2")))'   || { echo "FATAL: Rscript/dada2 not usable; aborting before staging" >&2; exit 1; }

# Stage the 8.45 GB artifact to local scratch.
echo "=== $(date) === staging artifact ==="
cp "${SRC}" "${SCRATCH}/paired-end-demux.qza"
echo "  staged $(du -h "${SCRATCH}/paired-end-demux.qza" | cut -f1)"

cd "${SCRATCH}"

# ── 1. DADA2 ──────────────────────────────────────────────────────────────
echo "=== $(date) === DADA2 denoise-paired ==="
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${SCRATCH}/paired-end-demux.qza" \
  --p-trim-left-f 17 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 245 \
  --p-trunc-len-r 245 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-trunc-q 2 \
  --p-chimera-method consensus \
  --p-n-threads "${SLURM_CPUS_PER_TASK}" \
  --output-dir "${SCRATCH}/DADA2_denoising_output" \
  --verbose

echo "=== $(date) === DADA2 done ==="
ls -la "${SCRATCH}/DADA2_denoising_output"

# Copy the DADA2 products back immediately, before attempting the tree, so a
# later failure cannot cost the expensive step.
cp "${SCRATCH}/DADA2_denoising_output"/*.qza "${OUT}/"
echo "  DADA2 artifacts copied to ${OUT}"

# ── 2. Phylogeny ──────────────────────────────────────────────────────────
# Same command as q2_dada2_KORA.ipynb cell 27: MAFFT align -> mask -> FastTree
# -> midpoint root.
echo "=== $(date) === align-to-tree-mafft-fasttree ==="
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${SCRATCH}/DADA2_denoising_output/representative_sequences.qza" \
  --p-n-threads "${SLURM_CPUS_PER_TASK}" \
  --output-dir "${SCRATCH}/tree" \
  --verbose

cp "${SCRATCH}/tree"/*.qza "${OUT}/"

# ── 3. Export to plain formats ────────────────────────────────────────────
echo "=== $(date) === exporting ==="
qiime tools export --input-path "${SCRATCH}/tree/rooted_tree.qza" \
  --output-path "${SCRATCH}/export_tree"
cp "${SCRATCH}/export_tree/tree.nwk" "${OUT}/rooted_tree.nwk"

qiime tools export \
  --input-path "${SCRATCH}/DADA2_denoising_output/representative_sequences.qza" \
  --output-path "${SCRATCH}/export_seqs"
cp "${SCRATCH}/export_seqs/dna-sequences.fasta" "${OUT}/rep_seqs.fasta"

qiime tools export \
  --input-path "${SCRATCH}/DADA2_denoising_output/table.qza" \
  --output-path "${SCRATCH}/export_table"
biom convert -i "${SCRATCH}/export_table/feature-table.biom" \
  -o "${OUT}/feature_table_reprocessed.tsv" --to-tsv

qiime metadata tabulate \
  --m-input-file "${SCRATCH}/DADA2_denoising_output/denoising_stats.qza" \
  --o-visualization "${OUT}/denoising_stats.qzv"

# ── 4. Verify the ASVs reproduce ──────────────────────────────────────────
# The whole point: tip labels must match the ASV ids the manuscript analyses use.
echo "=== $(date) === verifying ASV reproduction ==="
python - <<'PY'
import os
kora = "/home/itg/oleg.vlasovets/slr_example/KORA_data"
out = f"{kora}/q2_reprocessing"

def ids_from_tsv(path, skip):
    s = set()
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i < skip or line.startswith("#"):
                continue
            t = line.split("\t", 1)[0].strip()
            if t:
                s.add(t)
    return s

old = ids_from_tsv(f"{kora}/feature_table.tsv", 1)
new = ids_from_tsv(f"{out}/feature_table_reprocessed.tsv", 1)

tips = set()
with open(f"{out}/rooted_tree.nwk") as fh:
    import re
    tips = set(re.findall(r"[(,]([0-9a-f]{32})", fh.read()))

print(f"  published feature table : {len(old):6d} ASVs")
print(f"  reprocessed             : {len(new):6d} ASVs")
print(f"  tree tips               : {len(tips):6d}")
inter = old & new
print(f"  shared with published   : {len(inter):6d} "
      f"({100*len(inter)/max(1,len(old)):.1f}% of published)")
cov = old & tips
print(f"  published ASVs on tree  : {len(cov):6d} "
      f"({100*len(cov)/max(1,len(old)):.1f}%)")
if len(cov) / max(1, len(old)) > 0.99:
    print("  VERDICT: tree covers the published ASVs. UniFrac/Faith's PD are usable "
          "directly on the existing analyses.")
elif len(cov) / max(1, len(old)) > 0.80:
    print("  VERDICT: partial coverage. Usable, but report the covered fraction and "
          "restrict phylogenetic metrics to ASVs present on the tree.")
else:
    print("  VERDICT: POOR overlap. Do NOT use this tree with the published feature "
          "table - investigate before proceeding.")
PY

echo "=== $(date) === cleaning scratch ==="
rm -rf "${SCRATCH}"
echo "=== $(date) === Done. Outputs in ${OUT} ==="
ls -la "${OUT}"
