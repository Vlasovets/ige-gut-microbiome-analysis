#!/usr/bin/env Rscript
# Module 02d - three ADDITIONAL differential-abundance methods, so the
# cross-method comparison spans more than one family of assumptions.
#
# WHY
#   Module 02 runs LinDA, ANCOM-BC2 and ALDEx2. All three are log-ratio methods;
#   ANCOM-BC2 turned out to be anti-conservative on this data (18.3% of taxa
#   significant vs 0.07% LinDA and 0.00% ALDEx2), so a reviewer can fairly ask
#   whether the near-null result is a property of that one family. This module
#   adds three methods with genuinely different machinery:
#
#     ZicoSeq      permutation-based, reference-frame, winsorized; makes no
#                  parametric distributional assumption      (GUniFrac 1.9)
#     Wilcoxon-CLR paired non-parametric test on CLR values; the assumption-light
#                  baseline used in the Nearing et al. 2022 benchmark. This is
#                  the ONLY method here that uses the 1:1 matching directly
#                  rather than adjusting for the covariates.
#     limma-voom   mean-variance-weighted linear model on log-CPM; the
#                  RNA-seq-derived approach, included because it is widely used
#                  in the microbiome literature despite the count model being a
#                  poor fit for compositional data          (limma 3.62.2)
#
#   All three are already installed in slr-env, so no environment mutation is
#   needed. MaAsLin2 / metagenomeSeq / edgeR / DESeq2 were considered and NOT
#   installed: bioconda ships them for R 4.3 or R 4.5 while slr-env is R 4.4.3,
#   so installing any of them would force an R version change and invalidate the
#   environment behind the already-computed results. See revision/README.md.
#
# DESIGN RULES (identical to module 02, so the outputs are comparable)
#   * Same input tables, already prevalence-filtered upstream and label-blind.
#   * Same covariates: ~ ige_status + sex + age + bmi, EXCEPT Wilcoxon-CLR,
#     which uses the matched pairs instead (documented in the output).
#   * Effect sizes are NOT averaged across methods - the scales differ. Each is
#     recorded with its own effect_type.
#   * BH within each analysis, matching module 02's q_within.
#
# Usage:
#   Rscript revision/analysis/02d_extended_da.R \
#     --level species --prevalence prev05 --partition any_ige

suppressPackageStartupMessages({
  library(yaml); library(GUniFrac); library(limma)
})

args_raw <- commandArgs(trailingOnly = FALSE)
script <- sub("--file=", "", args_raw[grep("--file=", args_raw)])
REVISION <- normalizePath(file.path(dirname(script), ".."))
.libPaths(c(file.path(REVISION, "env", "rlibs"), .libPaths()))

a <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(a == flag)
  if (length(i) && i < length(a)) a[i + 1] else default
}
LEVEL <- getarg("--level", "species")
PREV  <- getarg("--prevalence", "prev05")
PART  <- getarg("--partition", "any_ige")
NPERM <- as.integer(getarg("--nperm", "4999"))

cfg  <- yaml::read_yaml(file.path(REVISION, "config", "config.yaml"))
PIPE <- cfg$paths$pipeline_results
OUT  <- file.path(REVISION, "results", "da_extended")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED <- cfg$statistics$seed
set.seed(SEED)
FORMULA <- "ige_status + sex + age + bmi"

message(sprintf("=== Extended DA: %s / %s / %s ===", LEVEL, PREV, PART))
message(sprintf("  formula: ~ %s | ZicoSeq nperm=%d | seed=%d",
                FORMULA, NPERM, SEED))

# ── Load exactly what module 02 loads ────────────────────────────────────
src <- file.path(PIPE, LEVEL, PREV, PART, "source")
oi <- read.csv(file.path(src, "otu_ige.csv"),     row.names = 1, check.names = FALSE)
oc <- read.csv(file.path(src, "otu_non_ige.csv"), row.names = 1, check.names = FALSE)
mi <- read.csv(file.path(src, "sample_data_ige.csv"),     row.names = 1)
mc <- read.csv(file.path(src, "sample_data_non_ige.csv"), row.names = 1)

taxa <- intersect(colnames(oi), colnames(oc))
otu  <- rbind(oi[, taxa, drop = FALSE], oc[, taxa, drop = FALSE])
meta <- rbind(mi, mc)
meta$ige_status <- factor(c(rep("ige", nrow(mi)), rep("non_ige", nrow(mc))),
                          levels = c("non_ige", "ige"))
meta$sex <- factor(meta$u3csex)
meta$age <- as.numeric(meta$u3talter)
meta$bmi <- as.numeric(meta$u3tbmi)

keep <- !is.na(meta$age) & !is.na(meta$bmi) & rowSums(otu) > 0
otu <- otu[keep, , drop = FALSE]; meta <- meta[keep, , drop = FALSE]
otu <- otu[rownames(meta), , drop = FALSE]
stopifnot(identical(rownames(otu), rownames(meta)))

n_ige <- sum(meta$ige_status == "ige"); n_ctrl <- sum(meta$ige_status == "non_ige")
message(sprintf("  samples=%d (ige=%d, non_ige=%d)  taxa=%d",
                nrow(otu), n_ige, n_ctrl, ncol(otu)))

counts_pxn <- t(as.matrix(otu))          # taxa x samples
res_list <- list()

# ── 1. ZicoSeq ───────────────────────────────────────────────────────────
# ZicoSeq refuses to run if ANY feature is constant ("identical values (e.g. all
# 0s)"). Prevalence filtering happened upstream on the full 2034-sample cohort,
# so a feature can be non-rare cohort-wide yet all-zero inside a small
# partition. Left unhandled this aborted ZicoSeq for 104 of 180 analyses - and
# because the failure is caught and the script continues, it showed up as
# ZicoSeq "finding nothing" rather than as ZicoSeq never running. That is the
# absence-read-as-evidence trap, so it is fixed and counted, not tolerated.
#
# Constant features are uninformative for EVERY method, so dropping them costs
# no signal. They are recorded as NA for ZicoSeq (never as non-significant),
# and module 02e computes pairwise overlap only over taxa both methods actually
# evaluated.
#
# TWO filters are needed, not one. Removing raw constants is not sufficient:
# ZicoSeq winsorizes the top `outlier.pct` of each feature before testing, and a
# sparse feature that is non-constant ONLY because of its few largest values
# becomes constant afterwards - e.g. (0,0,...,0,5) with one value replaced. That
# raised a second, different error ("After winsorization, feature ... have
# identical values") in 89 of 180 analyses even after the first fix.
#
# The winsorization rule is: values above the (1 - outlier.pct) quantile are
# replaced by that quantile. So with k = ceiling(n * outlier.pct) values
# replaced, the feature is constant afterwards exactly when the remaining
# n - k values are all equal. That is checked directly rather than guessed.
OUTLIER_PCT <- 0.03
n_samp <- ncol(counts_pxn)
k_wins <- ceiling(n_samp * OUTLIER_PCT)
survives <- apply(counts_pxn, 1, function(x) {
  if (length(unique(x)) < 2) return(FALSE)                 # raw constant
  v <- sort(x, decreasing = TRUE)
  length(unique(v[(k_wins + 1):length(v)])) > 1             # constant post-winsor
})
n_drop <- sum(!survives)
if (n_drop)
  message(sprintf("    %d of %d features are constant raw or after winsorization (k=%d of n=%d) - excluded from ZicoSeq",
                  n_drop, nrow(counts_pxn), k_wins, n_samp))
zico_in <- counts_pxn[survives, , drop = FALSE]
zi <- if (nrow(zico_in) < 5) {
  message("    fewer than 5 non-constant features; skipping ZicoSeq"); NULL
} else tryCatch(
  GUniFrac::ZicoSeq(
    meta.dat = meta, feature.dat = zico_in,
    grp.name = "ige_status", adj.name = c("sex", "age", "bmi"),
    feature.dat.type = "count",
    # Filtering is already done upstream; disable ZicoSeq's own so the feature
    # set matches every other method exactly.
    prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0,
    is.winsor = TRUE, outlier.pct = 0.03, winsor.end = "top",
    is.post.sample = TRUE, post.sample.no = 25,
    link.func = list(function(x) x^0.5),
    stats.combine.func = max,
    perm.no = NPERM, strata = NULL,
    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
    verbose = FALSE, return.feature.dat = FALSE),
  error = function(e) { message("    ERROR (winsorized): ", conditionMessage(e)); NULL })

# Last-resort fallback. If the exact pre-filter still leaves ZicoSeq unable to
# winsorize, re-run it WITHOUT winsorization rather than let the method silently
# disappear from this analysis - "did not run" was previously indistinguishable
# from "found nothing", which is the failure mode this whole module guards.
# The deviation is recorded per analysis in zico_winsorized so it can never be
# read as if every analysis used identical settings.
ZICO_WINSORIZED <- TRUE
if (is.null(zi) && nrow(zico_in) >= 5) {
  message("    retrying ZicoSeq with is.winsor = FALSE")
  ZICO_WINSORIZED <- FALSE
  zi <- tryCatch(
    GUniFrac::ZicoSeq(
      meta.dat = meta, feature.dat = zico_in,
      grp.name = "ige_status", adj.name = c("sex", "age", "bmi"),
      feature.dat.type = "count",
      prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0,
      is.winsor = FALSE,
      is.post.sample = TRUE, post.sample.no = 25,
      link.func = list(function(x) x^0.5), stats.combine.func = max,
      perm.no = NPERM, strata = NULL,
      ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
      verbose = FALSE, return.feature.dat = FALSE),
    error = function(e) { message("    ERROR (unwinsorized): ", conditionMessage(e)); NULL })
}
if (!is.null(zi)) {
  # ZicoSeq's R2 is unsigned; recover direction from the mean CLR difference so
  # sign concordance can be assessed against the other methods.
  cl <- t(apply(counts_pxn + 0.5, 2, function(x) log(x) - mean(log(x))))  # samples x taxa
  dir <- colMeans(cl[meta$ige_status == "ige", , drop = FALSE]) -
         colMeans(cl[meta$ige_status == "non_ige", , drop = FALSE])
  tx <- names(zi$p.raw)
  # ZicoSeq's p.adj.fdr is a PERMUTATION-based FDR, not BH: it pools statistics
  # across features, so q < p is legitimate for large p and is NOT an error.
  # Two columns are therefore reported:
  #   q_within  BH(p.raw)  - comparable with every other method and with module
  #                          02, which is what the overlap analysis requires
  #   q_native  ZicoSeq's own FDR - the method as its authors intend it
  # Comparing methods that use DIFFERENT multiplicity procedures would make any
  # "overlap" partly an artifact of the correction rather than of the test.
  res_list$ZicoSeq <- data.frame(
    taxon = tx, method = "ZicoSeq",
    effect_log2 = as.numeric(dir[tx]) / log(2),
    effect_type = "CLR_mean_difference_log2_direction_only",
    se = NA_real_, stat = as.numeric(zi$R2[tx]),
    p_raw = as.numeric(zi$p.raw[tx]),
    q_within = p.adjust(as.numeric(zi$p.raw[tx]), method = "BH"),
    q_native = as.numeric(zi$p.adj.fdr[tx]),
    passed_ss = NA, stringsAsFactors = FALSE)
  message(sprintf("    ZicoSeq evaluated %d of %d features (winsorized=%s)",
                  length(tx), nrow(counts_pxn), ZICO_WINSORIZED))
}

# ── 2. Paired Wilcoxon on CLR ────────────────────────────────────────────
# The assumption-light baseline. Uses the 1:1 matching directly: each sensitized
# individual is compared with the control they were matched to on sex (exact),
# age (+-5) and BMI (+-4), so the covariates are handled by design rather than
# by adjustment. This is a DIFFERENT route to the same confounding control, and
# is reported as such.
message("\n  [2/3] Paired Wilcoxon on CLR ...")
wil <- tryCatch({
  stopifnot("pair_nb" %in% names(meta))
  clr <- t(apply(counts_pxn + 0.5, 2, function(x) log(x) - mean(log(x))))
  ig <- meta$ige_status == "ige"
  pi <- meta$pair_nb[ig];  pc <- meta$pair_nb[!ig]
  common <- intersect(pi, pc)
  message(sprintf("    complete matched pairs: %d", length(common)))
  if (length(common) < 5) stop("too few complete pairs for a paired test")
  A <- clr[ig, , drop = FALSE][match(common, pi), , drop = FALSE]
  B <- clr[!ig, , drop = FALSE][match(common, pc), , drop = FALSE]
  D <- A - B                                    # pairs x taxa, CLR difference
  p <- apply(D, 2, function(d) {
    d <- d[is.finite(d)]
    if (length(d) < 5 || all(abs(d) < 1e-12)) return(NA_real_)
    tryCatch(stats::wilcox.test(d, mu = 0, exact = FALSE)$p.value,
             error = function(e) NA_real_)
  })
  data.frame(
    taxon = colnames(D), method = "Wilcoxon_CLR_paired",
    effect_log2 = colMeans(D, na.rm = TRUE) / log(2),
    effect_type = "paired_CLR_difference_log2",
    se = apply(D, 2, function(d) sd(d, na.rm = TRUE) / sqrt(sum(is.finite(d)))),
    stat = NA_real_, p_raw = p,
    q_within = p.adjust(p, method = "BH"), q_native = NA_real_,
    passed_ss = NA, stringsAsFactors = FALSE)
}, error = function(e) { message("    ERROR: ", conditionMessage(e)); NULL })
if (!is.null(wil)) res_list$Wilcoxon_CLR_paired <- wil

# ── 3. limma-voom ────────────────────────────────────────────────────────
# Included deliberately as a method whose count model is a POOR fit for
# compositional data. If it behaves like ANCOM-BC2 (many hits) while the
# log-ratio methods do not, that is informative about method calibration rather
# than about biology - which is the point being made to the reviewer.
message("\n  [3/3] limma-voom ...")
lv <- tryCatch({
  des <- model.matrix(as.formula(paste("~", FORMULA)), data = meta)
  stopifnot(nrow(des) == nrow(meta))
  cn <- grep("ige_status", colnames(des), value = TRUE)[1]
  # voom needs library sizes; log-CPM with quantile normalisation between
  # samples. No edgeR, so no TMM - stated in effect_type.
  v <- limma::voom(counts_pxn, design = des, normalize.method = "quantile",
                   plot = FALSE)
  fit <- limma::eBayes(limma::lmFit(v, des))
  tt <- limma::topTable(fit, coef = cn, number = Inf, sort.by = "none",
                        adjust.method = "BH")
  data.frame(
    taxon = rownames(tt), method = "limma_voom",
    effect_log2 = tt$logFC, effect_type = "log2FC_logCPM_quantile_norm",
    se = NA_real_, stat = tt$t, p_raw = tt$P.Value,
    q_within = tt$adj.P.Val, q_native = NA_real_,
    passed_ss = NA, stringsAsFactors = FALSE)
}, error = function(e) { message("    ERROR: ", conditionMessage(e)); NULL })
if (!is.null(lv)) res_list$limma_voom <- lv

if (!length(res_list)) { message("No method produced output. Aborting."); quit(status = 1) }
long <- do.call(rbind, res_list); rownames(long) <- NULL

# ── QC 1: BH sanity ──────────────────────────────────────────────────────
# Applies to q_within only. q_native is a permutation FDR for which q < p is
# legitimate, so including it here would fail the check for a non-reason.
bad <- long[!is.na(long$p_raw) & !is.na(long$q_within) &
              long$q_within < long$p_raw - 1e-12, ]
message(sprintf("\n  QC: rows with q_within < p: %d  %s", nrow(bad),
                if (nrow(bad)) "FAIL" else "PASS"))
if (nrow(bad)) quit(status = 1)

# ── QC 2: permutation floor ──────────────────────────────────────────────
# A permutation p has a hard floor of 1/(nperm+1). Under BH over m features the
# smallest attainable q is therefore m/(nperm+1). If that exceeds alpha, the
# ZicoSeq BH column is STRUCTURALLY incapable of significance and reporting
# "no hits" from it would be meaningless rather than evidential. This exact
# trap already invalidated an earlier network permutation test in this project
# at nPerm=200, so it is checked explicitly rather than assumed.
# m is the number of features ZicoSeq actually tested, not the number present:
# BH is applied over the tested set, so using the full count would overstate
# the floor.
m <- if (!is.null(zi)) length(zi$p.raw) else ncol(otu)
floor_q <- m / (NPERM + 1)
message(sprintf("  QC: permutation floor - min attainable BH q for ZicoSeq = %.4g (m=%d, nperm=%d)",
                floor_q, m, NPERM))
if (floor_q > 0.05) {
  message(sprintf("    WARNING: floor %.3g > 0.05. The ZicoSeq q_within column CANNOT", floor_q))
  message( "    reach significance at this nperm. Interpret ZicoSeq via q_native")
  message(sprintf("    (permutation FDR, not floor-limited this way), or raise nperm above %d.",
                  ceiling(m / 0.05) - 1))
}
long$zicoseq_bh_floor <- floor_q

# QC: every method must have returned the same feature set, else "overlap" is
# an artifact of differing denominators - the exact trap module 02 guards.
nt <- tapply(long$taxon, long$method, function(x) length(unique(x)))
message("  QC: taxa per method: ", paste(sprintf("%s=%d", names(nt), nt), collapse = " "))
if (length(unique(nt)) > 1)
  message("  NOTE: feature sets differ (ZicoSeq drops constant features). ",
          "Module 02e computes overlap only over taxa BOTH methods evaluated, ",
          "so this does not bias the comparison - but a bare count of ",
          "significant hits per method is not comparable across methods.")
if (length(res_list) < 3)
  message(sprintf("  NOTE: only %d of 3 methods produced output (%s). A method that ",
                  length(res_list), paste(names(res_list), collapse = ", ")),
          "did not run is NOT a method that found nothing.")

long$level <- LEVEL; long$prevalence <- PREV; long$partition <- PART
long$n_ige <- n_ige; long$n_ctrl <- n_ctrl; long$n_taxa <- ncol(otu)
long$formula <- paste("~", FORMULA); long$nperm <- NPERM; long$seed <- SEED
long$zico_winsorized <- if (exists("ZICO_WINSORIZED")) ZICO_WINSORIZED else NA

f <- file.path(OUT, sprintf("daext_long_%s_%s_%s.csv", LEVEL, PREV, PART))
write.csv(long, f, row.names = FALSE)
message(sprintf("\nWrote %s (%d rows)", f, nrow(long)))

message("\n  Significant at q_within < 0.05, by method:")
for (m in unique(long$method)) {
  s <- long[long$method == m & !is.na(long$q_within) & long$q_within < 0.05, ]
  detail <- if (nrow(s) == 0) "" else paste0(": ", paste(
    sprintf("%s(%+.2f)", s$taxon, s$effect_log2), collapse = ", "))
  if (nchar(detail) > 300) detail <- paste0(substr(detail, 1, 300), " ...")
  message(sprintf("    %-20s %d%s", m, nrow(s), detail))
}
message("Done.")
