#!/usr/bin/env Rscript
# Module 07b - one step of the M0 -> M3 nested-model chain, at ASV rank.
#
#   M0  published design: per-partition, unadjusted (~w), controls down-sampled
#       to partition size (except any_ige), one draw per seed  [needs --seed]
#   M1  same, but ALL 233 controls for every partition (removes down-sampling)
#   M2  matched 221-pair design, unadjusted (~w)               (adds matching)
#   M3  matched + sex/age/BMI in the formula                   (adds adjustment)
#
# Everything else is held FIXED across the chain, so each step isolates one
# design change: counts always come from the raw feature table; the feature set
# is always "ASVs non-zero within the selected samples" (the published filter,
# freq_threshold = 0, with zero-mean samples dropped); the test is always
# LinDA with the published arguments (formula, alpha = 0.05, prev.cut = 0,
# lib.cut = 1, package defaults otherwise - exactly what
# utils/helper.calculate_obs_stat_and_pvalues passed); multiplicity is always
# BH per partition at 0.05, hits pooled as the UNION across partitions -
# which is how the published "61 ASVs" was counted.
#
# Usage:
#   Rscript revision/analysis/07b_nested_linda.R --model M0 --seed 42
#   Rscript revision/analysis/07b_nested_linda.R --model M1|M2|M3

suppressPackageStartupMessages({ library(yaml); library(LinDA) })

args_raw <- commandArgs(trailingOnly = FALSE)
script <- sub("--file=", "", args_raw[grep("--file=", args_raw)])
REVISION <- normalizePath(file.path(dirname(script), ".."))
.libPaths(c(file.path(REVISION, "env", "rlibs"), .libPaths()))

a <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(a == flag); if (length(i) && i < length(a)) a[i + 1] else default
}
MODEL <- getarg("--model"); SEED <- getarg("--seed")
stopifnot(MODEL %in% c("M0", "M1", "M2", "M3"))
if (MODEL == "M0") stopifnot(!is.null(SEED))

cfg  <- yaml::read_yaml(file.path(REVISION, "config", "config.yaml"))
KORA <- cfg$paths$kora_data
PIPE <- cfg$paths$pipeline_results
DES  <- file.path(REVISION, "results", "nested_models", "designs")
OUT  <- file.path(REVISION, "results", "nested_models", MODEL)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ALPHA <- 0.05
PARTITIONS <- c("any_ige", "single", "multiple", "high_score", "medium_score",
                "low_score", paste0("cluster_", LETTERS[1:9]))

message(sprintf("=== 07b nested LinDA: %s%s ===", MODEL,
                if (!is.null(SEED)) paste0(" seed ", SEED) else ""))

ft <- read.csv(file.path(KORA, "feature_table.tsv"), sep = "\t",
               row.names = 1, check.names = FALSE)
colnames(ft) <- as.character(colnames(ft))

base <- read.csv(file.path(DES, "base.csv"), colClasses = "character")
ctrl_pool <- base$sample_id[base$partition == "__controls__"]
members <- function(p) base$sample_id[base$partition == p & base$role == "case"]

# One LinDA run under the published pipeline rules.
run_one <- function(part, case_ids, ctrl_ids, covar = NULL) {
  ids <- intersect(c(case_ids, ctrl_ids), colnames(ft))
  M <- as.matrix(ft[, ids, drop = FALSE])
  M <- M[rowSums(M > 0) > 0, , drop = FALSE]          # published filter, thr=0
  keep_s <- colMeans(M) > 0                           # drop zero-mean samples
  M <- M[, keep_s, drop = FALSE]
  meta <- data.frame(row.names = colnames(M),
                     w = as.integer(colnames(M) %in% case_ids))
  fml <- "~w"
  if (!is.null(covar)) {
    meta <- cbind(meta, covar[colnames(M), , drop = FALSE])
    meta$u3csex <- factor(meta$u3csex)
    meta <- meta[stats::complete.cases(meta), , drop = FALSE]
    M <- M[, rownames(meta), drop = FALSE]
    M <- M[rowSums(M > 0) > 0, , drop = FALSE]
    fml <- "~w+u3csex+u3talter+u3tbmi"
  }
  if (length(unique(meta$w)) < 2 || nrow(M) < 20) return(NULL)
  lo <- LinDA::linda(M, meta, formula = fml, alpha = ALPHA,
                     prev.cut = 0, lib.cut = 1)
  o <- lo$output$w
  q <- p.adjust(o$pvalue, method = "BH")              # published BH-per-partition
  data.frame(model = MODEL, seed = ifelse(is.null(SEED), NA, SEED),
             partition = part, asv = rownames(o),
             lfc = o$log2FoldChange, p = o$pvalue, q = q,
             n_case = sum(meta$w == 1), n_ctrl = sum(meta$w == 0),
             n_asv = nrow(M), stringsAsFactors = FALSE)
}

res <- list()
if (MODEL == "M0") {
  draw <- read.csv(file.path(DES, sprintf("seed_%s.csv", SEED)),
                   colClasses = "character")
  # any_ige has no down-sampling, so it is identical in every seed; it is run
  # once in M1 and merged back by 07c. Skipping it here saves 100 LinDA runs.
  for (p in setdiff(PARTITIONS, "any_ige")) {
    res[[p]] <- run_one(p, members(p), draw$sample_id[draw$partition == p])
    message(sprintf("  %-13s hits=%d", p, sum(res[[p]]$q < ALPHA, na.rm = TRUE)))
  }
} else if (MODEL == "M1") {
  for (p in PARTITIONS) {
    res[[p]] <- run_one(p, members(p), ctrl_pool)
    message(sprintf("  %-13s hits=%d", p, sum(res[[p]]$q < ALPHA, na.rm = TRUE)))
  }
} else {                                              # M2 / M3: matched pairs
  for (p in PARTITIONS) {
    src <- file.path(PIPE, "asv", "prev01", p, "source")
    mi <- read.csv(file.path(src, "sample_data_ige.csv"),     row.names = 1)
    mc <- read.csv(file.path(src, "sample_data_non_ige.csv"), row.names = 1)
    cp <- sort(intersect(mi$pair_nb, mc$pair_nb))
    mi <- mi[match(cp, mi$pair_nb), , drop = FALSE]
    mc <- mc[match(cp, mc$pair_nb), , drop = FALSE]
    if (length(cp) < 5) { message(sprintf("  %-13s <5 pairs, skipped", p)); next }
    covar <- NULL
    if (MODEL == "M3") {
      covar <- rbind(mi, mc)[, c("u3csex", "u3talter", "u3tbmi")]
      rownames(covar) <- as.character(c(rownames(mi), rownames(mc)))
    }
    res[[p]] <- run_one(p, as.character(rownames(mi)),
                        as.character(rownames(mc)), covar)
    if (!is.null(res[[p]]))
      message(sprintf("  %-13s pairs=%d hits=%d", p, length(cp),
                      sum(res[[p]]$q < ALPHA, na.rm = TRUE)))
  }
}

long <- do.call(rbind, res); rownames(long) <- NULL
stopifnot(all(long$q >= long$p - 1e-12, na.rm = TRUE))   # BH sanity

hits <- long[!is.na(long$q) & long$q < ALPHA, ]
stem <- if (MODEL == "M0") sprintf("%s_seed%s", MODEL, SEED) else MODEL
write.csv(hits, file.path(OUT, paste0(stem, "_hits.csv")), row.names = FALSE)
summ <- aggregate(cbind(n_asv, n_case, n_ctrl) ~ model + partition, long, max)
summ$n_hits <- vapply(summ$partition,
                      function(p) sum(hits$partition == p), integer(1))
summ$seed <- ifelse(is.null(SEED), NA, SEED)
write.csv(summ, file.path(OUT, paste0(stem, "_summary.csv")), row.names = FALSE)
message(sprintf("\n  union of hit ASVs across partitions in this task: %d",
                length(unique(hits$asv))))
message(sprintf("Wrote %s_{hits,summary}.csv", file.path(OUT, stem)))
