#!/bin/bash
#SBATCH --job-name=rev_rinstall
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=3:00:00
#SBATCH --output=revision/slurm/logs/rinstall_%j.out
#SBATCH --error=revision/slurm/logs/rinstall_%j.err

# Install revision-only R packages into a SEPARATE library path.
#
# Why not a cloned conda env: `conda create --clone` DEADLOCKS on this cluster
# (hangs after "Retrieving notices"). A separate R library gives the same
# isolation instantly - slr-env is never modified, so the previously published
# results stay reproducible, and revision scripts see the extra packages by
# prepending this directory to .libPaths().
#
# Bioconductor here is 3.20 (pinned by R 4.4.3):
#   ANCOMBC, ALDEx2 -> available in 3.20 via BiocManager
#   TaxSEA          -> NOT in 3.20 (entered at 3.21), install from GitHub
#   UpSetR          -> CRAN
#
# Usage (from the repo root):
#   sbatch revision/slurm/00_install_r_extras.sh

set -euo pipefail

export REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
RLIBS="${REPO}/revision/env/rlibs"
SCRATCH="/localscratch/${USER}/rev_rinstall_${SLURM_JOB_ID:-manual}"

mkdir -p "${RLIBS}" "${SCRATCH}"
export TMPDIR="${SCRATCH}"

cd "${REPO}"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate slr-env

echo "=== $(date) === R package install into a separate library ==="
echo "  R:        $(which R)"
echo "  lib:      ${RLIBS}"
echo "  TMPDIR:   ${TMPDIR}"
echo "  cpus:     ${SLURM_CPUS_PER_TASK:-1}"

Rscript --vanilla - <<'RSCRIPT'
rlibs <- file.path(Sys.getenv("REPO"), "revision", "env", "rlibs")
if (!dir.exists(rlibs)) dir.create(rlibs, recursive = TRUE)
.libPaths(c(rlibs, .libPaths()))

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
options(Ncpus = ncpus, repos = c(CRAN = "https://cloud.r-project.org"))

cat("Library paths:\n"); print(.libPaths())
cat(sprintf("R %s | Bioconductor %s\n",
            getRversion(), as.character(BiocManager::version())))

want_bioc <- c("ANCOMBC", "ALDEx2")
want_cran <- c("UpSetR")

install_if_missing <- function(pkgs, installer) {
  for (p in pkgs) {
    if (requireNamespace(p, quietly = TRUE)) {
      cat(sprintf("  SKIP  %s already available (%s)\n",
                  p, as.character(packageVersion(p))))
      next
    }
    cat(sprintf("  INSTALL  %s ...\n", p))
    ok <- tryCatch({ installer(p); TRUE },
                   error = function(e) { cat("    ERROR: ", conditionMessage(e), "\n"); FALSE })
    if (ok && requireNamespace(p, quietly = TRUE)) {
      cat(sprintf("    OK  %s %s\n", p, as.character(packageVersion(p))))
    } else {
      cat(sprintf("    FAILED  %s\n", p))
    }
  }
}

install_if_missing(want_cran, function(p)
  install.packages(p, lib = rlibs, dependencies = TRUE))

install_if_missing(want_bioc, function(p)
  BiocManager::install(p, lib = rlibs, ask = FALSE, update = FALSE))

# TaxSEA is not in Bioconductor 3.20. Install from GitHub at a pinned ref so the
# version behind the re-run is recorded - the version behind the ORIGINAL run is
# unknown, which is part of why the p/FDR error cannot be diagnosed from disk.
if (!requireNamespace("TaxSEA", quietly = TRUE)) {
  cat("  INSTALL  TaxSEA from GitHub (absent from Bioc 3.20) ...\n")
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", lib = rlibs)
  tryCatch(
    remotes::install_github("feargalr/TaxSEA", lib = rlibs,
                            upgrade = "never", dependencies = TRUE),
    error = function(e) cat("    ERROR: ", conditionMessage(e), "\n"))
}

cat("\n=== Verification ===\n")
for (p in c("ANCOMBC", "ALDEx2", "TaxSEA", "UpSetR",
            "LinDA", "vegan", "NetCoMi", "SpiecEasi", "GUniFrac", "lmerTest")) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")
  cat(sprintf("  %-12s %s\n", p, v))
}

# ancombc2() is the function the revision needs, not ancombc().
if (requireNamespace("ANCOMBC", quietly = TRUE))
  cat(sprintf("\n  ancombc2 exported: %s\n",
              "ancombc2" %in% getNamespaceExports("ANCOMBC")))

writeLines(capture.output(sessionInfo()),
           file.path(Sys.getenv("REPO"), "revision", "env",
                     "sessionInfo_rlibs.txt"))
RSCRIPT

echo "=== $(date) === cleaning scratch ==="
rm -rf "${SCRATCH}"
echo "=== $(date) === Done ==="
