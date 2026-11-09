#!/usr/bin/env Rscript
# Module 02j - what makes ANCOM-BC2 anti-conservative on this data? Ablation on
# permuted labels (module 02i showed: median 262 of 489 ASVs at q < 0.05 with
# permuted IgE labels; false-positive rate 74% for ASVs present in <= 10% of
# participants, 8% for ASVs present in > 80%).
#
# Hypothesis: zero handling. ancombc2's default pseudo = 0 takes log(count) and
# sets log(0) = NA, so each taxon is fitted only in its carriers; ancombc (v1),
# which is calibrated here, uses log(count + 1). Each configuration changes ONE
# setting relative to the pipeline (module 02):
#   base     pipeline settings (pseudo = 0, group + struc_zero + neg_lb, covariates)
#   pseudo1  pseudo = 1                     (zeros kept, as in ancombc v1)
#   nostruc  struc_zero = FALSE, neg_lb = FALSE, group = NULL
#   nocov    fix_formula = "ige_status"     (no sex/age/BMI)
# pseudo_sens = FALSE throughout: it only produces the passed_ss flag and does
# not change p-values (checked: base perm 1 vs module 02i perm 1).
# Permutations use the same seeds as 02i (1000 + k); perm 0 = real labels.
#
# Data: primary contrast, ASV, 5% prevalence filter, full cohort.
# Output: results/calibration/ablation_<config>_perm<k>.csv
# Usage:  Rscript revision/analysis/02j_ancombc2_ablation.R --config pseudo1 --perm 3

suppressPackageStartupMessages({ library(yaml); library(ANCOMBC); library(TreeSummarizedExperiment) })
args_raw <- commandArgs(trailingOnly = FALSE)
script <- sub("--file=", "", args_raw[grep("--file=", args_raw)])
REVISION <- normalizePath(file.path(dirname(script), ".."))
.libPaths(c(file.path(REVISION, "env", "rlibs"), .libPaths()))
a <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(a == flag); if (length(i) && i < length(a)) a[i + 1] else default
}
CONFIG <- getarg("--config", "base"); PERM <- as.integer(getarg("--perm", "1"))
stopifnot(CONFIG %in% c("base", "pseudo1", "nostruc", "nocov"))
SEED <- yaml::read_yaml(file.path(REVISION, "config", "config.yaml"))$statistics$seed
OUT <- file.path(REVISION, "results", "calibration")
message(sprintf("=== 02j ANCOM-BC2 ablation: %s, perm %d ===", CONFIG, PERM))

src <- file.path(REVISION, "results", "full_cohort", "inputs", "asv", "prev05", "any_ige", "source")
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
if (PERM > 0) {                       # identical permutation to module 02i
  set.seed(1000 + PERM)
  meta$ige_status <- meta$ige_status[sample(nrow(meta))]
}
tse <- TreeSummarizedExperiment(assays = list(counts = t(as.matrix(otu))), colData = meta)

argl <- list(data = tse, assay_name = "counts", tax_level = NULL,
             fix_formula = "ige_status + sex + age + bmi", p_adj_method = "BH",
             prv_cut = 0, lib_cut = 0, group = "ige_status", struc_zero = TRUE,
             neg_lb = TRUE, alpha = 0.05, n_cl = 1, verbose = FALSE, global = FALSE,
             pairwise = FALSE, pseudo_sens = FALSE, pseudo = 0)
if (CONFIG == "pseudo1") argl$pseudo <- 1
if (CONFIG == "nostruc") { argl$struc_zero <- FALSE; argl$neg_lb <- FALSE; argl$group <- NULL }
if (CONFIG == "nocov")   argl$fix_formula <- "ige_status"

set.seed(SEED)
r <- suppressWarnings(suppressMessages(do.call(ancombc2, argl)))$res
d <- data.frame(taxon = r$taxon, effect_log2 = r$lfc_ige_statusige / log(2),
                se = r$se_ige_statusige / log(2), p_raw = r$p_ige_statusige,
                q_within = r$q_ige_statusige, config = CONFIG, perm = PERM)
f <- file.path(OUT, sprintf("ablation_%s_perm%02d.csv", CONFIG, PERM))
write.csv(d, f, row.names = FALSE)
message(sprintf("  tested %d  q<0.05: %d  p<0.05: %.2f", nrow(d),
                sum(d$q_within < 0.05, na.rm = TRUE), mean(d$p_raw < 0.05, na.rm = TRUE)))
