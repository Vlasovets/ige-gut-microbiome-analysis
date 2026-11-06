#!/bin/bash
#SBATCH --job-name=rev_taxseaenv
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=revision/slurm/logs/taxseaenv_%j.out
#SBATCH --error=revision/slurm/logs/taxseaenv_%j.err

# Create a dedicated environment for module 07 (TaxSEA re-computation).
#
# Why a separate env: TaxSEA hard-requires R >= 4.5.0. Verified by the failed
# install in job 39114950:
#   "ERROR: this R is version 4.4.3, package 'TaxSEA' requires R >= 4.5.0"
# slr-env is pinned to R 4.4.3 / Bioconductor 3.20 and must stay that way, since
# it produced the published results. bioconda ships a PRECOMPILED r45 build
# (bioconductor-taxsea 1.2.0 r45hdfd78af_0), so nothing is compiled here.
#
# `conda create` from scratch is used, NOT `conda create --clone`, which
# DEADLOCKS on this cluster.
#
# Module 07 only needs to take a ranked taxon vector and run enrichment, so this
# env is deliberately minimal - it is not a full analysis environment.
#
# Usage:  sbatch revision/slurm/00c_create_taxsea_env.sh

set -euo pipefail

export REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
ENV_NAME=taxsea-env
SCRATCH="/localscratch/${USER}/rev_taxseaenv_${SLURM_JOB_ID:-manual}"
mkdir -p "${SCRATCH}"
export TMPDIR="${SCRATCH}"

cd "${REPO}"
source ~/miniconda3/etc/profile.d/conda.sh

echo "=== $(date) === creating ${ENV_NAME} (R 4.5 + TaxSEA) ==="

if conda env list | grep -qE "^${ENV_NAME}\s"; then
  echo "  ${ENV_NAME} already exists; skipping creation"
else
  conda create -n "${ENV_NAME}" -y \
    -c conda-forge -c bioconda \
    bioconductor-taxsea \
    r-yaml
fi

conda activate "${ENV_NAME}"

echo ""
echo "=== $(date) === Verification ==="
Rscript --vanilla - <<'RSCRIPT'
cat(sprintf("  R version: %s\n", getRversion()))
if (requireNamespace("BiocManager", quietly = TRUE))
  cat(sprintf("  Bioconductor: %s\n", as.character(BiocManager::version())))
for (p in c("TaxSEA", "yaml")) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")
  cat(sprintf("  %-10s %s\n", p, v))
}
if (requireNamespace("TaxSEA", quietly = TRUE)) {
  ex <- getNamespaceExports("TaxSEA")
  cat(sprintf("\n  TaxSEA exports (%d): %s\n", length(ex),
              paste(head(sort(ex), 12), collapse = ", ")))
  # The reference sets ship with the package; record what is available so the
  # manuscript can state which database was used. The published analysis never
  # documented this.
  d <- tryCatch(data(package = "TaxSEA")$results[, "Item"],
                error = function(e) character(0))
  if (length(d)) cat(sprintf("  Bundled datasets: %s\n", paste(d, collapse = ", ")))
}
RSCRIPT

# Record the spec so the env is reproducible, per the project convention.
SPECDIR="${REPO}/revision/env"
mkdir -p "${SPECDIR}"
conda env export -n "${ENV_NAME}" --no-builds > "${SPECDIR}/taxsea-env.yml" 2>/dev/null || \
  echo "  (env export failed; record manually)"
python - <<'PY' || true
import json, os, glob
prefix = os.popen("conda env list | awk '$1==\"taxsea-env\"{print $NF}'").read().strip()
rows = []
for f in glob.glob(os.path.join(prefix, "conda-meta", "*.json")):
    try:
        d = json.load(open(f))
        rows.append((d.get("name"), d.get("version"), d.get("build"), d.get("channel", "")))
    except Exception:
        pass
out = os.path.join(os.environ["REPO"], "revision", "env", "taxsea-env.pkglist.tsv")
with open(out, "w") as fh:
    fh.write("name\tversion\tbuild\tchannel\n")
    for r in sorted(rows):
        fh.write("\t".join(str(x) for x in r) + "\n")
print(f"  wrote {out} ({len(rows)} packages)")
PY

rm -rf "${SCRATCH}"
echo "=== $(date) === Done ==="
