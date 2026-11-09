#!/usr/bin/env Rscript
# Module 02i - is the LinDA / ANCOM-BC2 disagreement a pipeline bug or a
# property of the methods? Calibration check on the primary contrast.
#
# Primary contrast, full cohort (275 vs 233), ~ ige_status + sex + age + bmi,
# 5% prevalence filter, at one rank. Three methods with IDENTICAL arguments to
# module 02, plus ANCOM-BC (v1), which a co-author expected to agree with LinDA:
#   LinDA       prev.cut = 0, lib.cut = 1, winsor.quan = NULL, pseudo.cnt = 0.5
#   ANCOM-BC2   group = ige_status, struc_zero, neg_lb, prv_cut = lib_cut = 0, BH
#   ANCOM-BC    same group/struc_zero/neg_lb/prv_cut/lib_cut, BH, package defaults
#               otherwise
#
# --perm 0   real labels: must reproduce module 02's stored LinDA and ANCOM-BC2
#            results exactly (reproducibility check), and adds ANCOM-BC v1.
# --perm k   IgE labels permuted across the 508 participants (seed 1000 + k);
#            covariates untouched. No true signal exists, so a calibrated method
#            calls ~0 taxa at q < 0.05 and has uniform p-values.
#
# Note on zeros (read from the ANCOMBC 2.8.0 source): with pseudo = 0 (default)
# ancombc2 takes log(count) and sets log(0) to NA, i.e. each taxon is modelled
# only in the samples where it is present. LinDA adds 0.5 to zeros instead.
#
# Output: results/calibration/calib_<level>_perm<k>.csv
# Usage:  Rscript revision/analysis/02i_calibration_check.R --level asv --perm 3

suppressPackageStartupMessages({
  library(yaml); library(LinDA); library(ANCOMBC); library(TreeSummarizedExperiment)
})
args_raw <- commandArgs(trailingOnly = FALSE)
script <- sub("--file=", "", args_raw[grep("--file=", args_raw)])
REVISION <- normalizePath(file.path(dirname(script), ".."))
.libPaths(c(file.path(REVISION, "env", "rlibs"), .libPaths()))
a <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(a == flag); if (length(i) && i < length(a)) a[i + 1] else default
}
LEVEL <- getarg("--level", "asv"); PERM <- as.integer(getarg("--perm", "0"))
PREV <- "prev05"; PART <- "any_ige"
FORMULA <- "ige_status + sex + age + bmi"
SEED <- yaml::read_yaml(file.path(REVISION, "config", "config.yaml"))$statistics$seed
OUT <- file.path(REVISION, "results", "calibration")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
message(sprintf("=== 02i calibration: %s / %s / %s, perm %d ===", LEVEL, PREV, PART, PERM))

src <- file.path(REVISION, "results", "full_cohort", "inputs", LEVEL, PREV, PART, "source")
oi <- read.csv(file.path(src, "otu_ige.csv"), row.names = 1, check.names = FALSE)
oc <- read.csv(file.path(src, "otu_non_ige.csv"), row.names = 1, check.names = FALSE)
mi <- read.csv(file.path(src, "sample_data_ige.csv"), row.names = 1)
mc <- read.csv(file.path(src, "sample_data_non_ige.csv"), row.names = 1)
taxa <- intersect(colnames(oi), colnames(oc))
otu <- rbind(oi[, taxa, drop = FALSE], oc[, taxa, drop = FALSE])
meta <- rbind(mi, mc)
meta$ige_status <- factor(c(rep("ige", nrow(mi)), rep("non_ige", nrow(mc))),
                          levels = c("non_ige", "ige"))
meta$sex <- factor(meta$u3csex); meta$age <- as.numeric(meta$u3talter)
meta$bmi <- as.numeric(meta$u3tbmi)
keep <- !is.na(meta$age) & !is.na(meta$bmi) & rowSums(otu) > 0
otu <- otu[keep, , drop = FALSE]; meta <- meta[keep, , drop = FALSE]
if (PERM > 0) {
  set.seed(1000 + PERM)
  meta$ige_status <- meta$ige_status[sample(nrow(meta))]
}
counts <- t(as.matrix(otu))
res <- list()

lin <- linda(otu.tab = counts, meta = transform(meta, ige_status = as.integer(ige_status == "ige")),
             formula = paste("~", FORMULA), alpha = 0.05, p.adj.method = "BH",
             prev.cut = 0, lib.cut = 1, winsor.quan = NULL, pseudo.cnt = 0.5,
             corr.cut = 0.1, n.cores = 1)
o <- lin$output[[grep("ige_status", names(lin$output), value = TRUE)[1]]]
res$LinDA <- data.frame(taxon = rownames(o), method = "LinDA", effect_log2 = o$log2FoldChange,
                        p_raw = o$pvalue, q_within = o$padj)

tse <- TreeSummarizedExperiment(assays = list(counts = counts), colData = meta)
set.seed(SEED)
ab <- suppressMessages(ancombc2(
  data = tse, assay_name = "counts", tax_level = NULL, fix_formula = FORMULA,
  p_adj_method = "BH", prv_cut = 0, lib_cut = 0, group = "ige_status",
  struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, n_cl = 1, verbose = FALSE,
  global = FALSE, pairwise = FALSE, pseudo_sens = TRUE))$res
ss <- grep("^passed_ss", names(ab), value = TRUE)
res$ANCOMBC2 <- data.frame(taxon = ab$taxon, method = "ANCOMBC2",
                           effect_log2 = ab$lfc_ige_statusige / log(2),
                           p_raw = ab$p_ige_statusige, q_within = ab$q_ige_statusige,
                           passed_ss = if (length(ss)) ab[[ss[1]]] else NA)

set.seed(SEED)
v1 <- suppressMessages(ancombc(
  data = tse, assay_name = "counts", formula = FORMULA,
  p_adj_method = "BH", prv_cut = 0, lib_cut = 0, group = "ige_status",
  struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = FALSE, n_cl = 1,
  verbose = FALSE))$res
col <- grep("ige_status", names(v1$lfc), value = TRUE)
message("  ANCOM-BC v1 lfc columns: ", paste(names(v1$lfc), collapse = ", "))
stopifnot(length(col) == 1)
tx <- if ("taxon" %in% names(v1$lfc)) v1$lfc$taxon else rownames(v1$lfc)
res$ANCOMBC <- data.frame(taxon = tx, method = "ANCOMBC",
                          # with struc_zero = TRUE ancombc() stores these as 1-column
                          # matrices named "structural_zero (...)": values are right,
                          # the name is not, so flatten
                          effect_log2 = as.numeric(v1$lfc[, col]) / log(2),
                          p_raw = as.numeric(v1$p_val[, col]),
                          q_within = as.numeric(v1$q_val[, col]))

COLS <- c("taxon", "method", "effect_log2", "p_raw", "q_within", "passed_ss")
long <- do.call(rbind, lapply(names(res), function(m) {
  x <- res[[m]]
  if (is.null(x$passed_ss)) x$passed_ss <- NA
  miss <- setdiff(COLS, names(x))
  if (length(miss)) stop(m, " lacks ", paste(miss, collapse = ","), "; has ", paste(names(x), collapse = ","))
  x[, COLS]
}))
long$level <- LEVEL; long$perm <- PERM
f <- file.path(OUT, sprintf("calib_%s_perm%02d.csv", LEVEL, PERM))
write.csv(long, f, row.names = FALSE)
for (m in names(res))
  message(sprintf("  %-9s tested %4d  q<0.05: %4d  p<0.05: %.2f", m, nrow(res[[m]]),
                  sum(res[[m]]$q_within < 0.05, na.rm = TRUE),
                  mean(res[[m]]$p_raw < 0.05, na.rm = TRUE)))
message("wrote ", f)
