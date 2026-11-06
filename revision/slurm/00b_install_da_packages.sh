#!/bin/bash
#SBATCH --job-name=rev_dainstall
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=revision/slurm/logs/dainstall_%j.out
#SBATCH --error=revision/slurm/logs/dainstall_%j.err

# Install ANCOMBC and ALDEx2 for the consensus DA module (03_consensus_da.R).
#
# Why conda rather than BiocManager into a separate R library:
#   The first attempt (00_install_r_extras.sh, job 39114435) compiled from
#   source and FAILED because the conda env lacks build tools -
#     RcppParallel -> "requires cmake (>= 3.5); cmake was not found"
#     highs        -> configuration failed (cmake)
#     energy       -> dependency 'gsl' not available
#     haven        -> fatal error: zlib.h: No such file or directory
#   bioconda ships PRECOMPILED builds against R 4.4 (r44 tag) at exactly the
#   Bioconductor 3.20 versions we want, so nothing is compiled.
#
# Safety: `conda install --dry-run` was inspected before running this. The
# transaction adds 48 new r44 packages and UPDATES ONLY ca-certificates,
# certifi and openssl. No analysis package, no numerical library and no r-base
# change - so LinDA/vegan/phyloseq/SpiecEasi/NetCoMi results are unaffected.
# Rollback point recorded below.
#
# ROLLBACK:  conda install -n slr-env --revision 3
#
# UpSetR already installed into revision/env/rlibs by job 39114435.
# TaxSEA is not on bioconda (it entered Bioconductor at 3.21, which needs
# R >= 4.5); it is installed from GitHub here now that its deps are satisfiable.

set -euo pipefail

export REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
RLIBS="${REPO}/revision/env/rlibs"
ENV_NAME=slr-env
SCRATCH="/localscratch/${USER}/rev_dainstall_${SLURM_JOB_ID:-manual}"

mkdir -p "${RLIBS}" "${SCRATCH}"
export TMPDIR="${SCRATCH}"
cd "${REPO}"

source ~/miniconda3/etc/profile.d/conda.sh

echo "=== $(date) === installing DA packages ==="
echo "  env:      ${ENV_NAME}"
echo "  rollback: conda install -n ${ENV_NAME} --revision 3"

conda list -n "${ENV_NAME}" --revisions | tail -5

conda install -n "${ENV_NAME}" -y \
  -c conda-forge -c bioconda \
  bioconductor-ancombc=2.8.0 \
  bioconductor-aldex2=1.38.0

conda activate "${ENV_NAME}"

echo ""
echo "=== $(date) === TaxSEA from GitHub (absent from Bioc 3.20) ==="
Rscript --vanilla - <<'RSCRIPT'
rlibs <- file.path(Sys.getenv("REPO"), "revision", "env", "rlibs")
.libPaths(c(rlibs, .libPaths()))
options(Ncpus = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
        repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("TaxSEA", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", lib = rlibs)
  tryCatch(
    remotes::install_github("feargalr/TaxSEA", lib = rlibs,
                            upgrade = "never", dependencies = TRUE),
    error = function(e) cat("  TaxSEA install ERROR: ", conditionMessage(e), "\n"))
}
RSCRIPT

echo ""
echo "=== $(date) === Verification ==="
Rscript --vanilla - <<'RSCRIPT'
rlibs <- file.path(Sys.getenv("REPO"), "revision", "env", "rlibs")
.libPaths(c(rlibs, .libPaths()))

pkgs <- c("ANCOMBC", "ALDEx2", "TaxSEA", "UpSetR",
          "LinDA", "MicrobiomeStat", "vegan", "phyloseq",
          "NetCoMi", "SpiecEasi", "permApprox", "GUniFrac", "lmerTest")
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")
  cat(sprintf("  %-16s %s\n", p, v))
}

cat("\n  Required entry points:\n")
for (spec in list(c("ANCOMBC", "ancombc2"), c("ALDEx2", "aldex.clr"),
                  c("ALDEx2", "aldex.glm"), c("LinDA", "linda"))) {
  ok <- tryCatch(spec[2] %in% getNamespaceExports(spec[1]),
                 error = function(e) FALSE)
  cat(sprintf("    %-10s::%-12s %s\n", spec[1], spec[2],
              if (ok) "available" else "MISSING"))
}

# Regression check: the previously published numbers must be unchanged after
# touching the environment. vegan and LinDA are the load-bearing packages.
cat("\n  Regression check (must match the pre-install values):\n")
cat(sprintf("    vegan   %s (expected 2.7-5)\n", as.character(packageVersion("vegan"))))
cat(sprintf("    LinDA   %s (expected 0.2.0)\n", as.character(packageVersion("LinDA"))))
cat(sprintf("    SpiecEasi %s (expected 1.1.1)\n", as.character(packageVersion("SpiecEasi"))))

writeLines(capture.output(sessionInfo()),
           file.path(Sys.getenv("REPO"), "revision", "env", "sessionInfo_da.txt"))
RSCRIPT

rm -rf "${SCRATCH}"
echo "=== $(date) === Done ==="
