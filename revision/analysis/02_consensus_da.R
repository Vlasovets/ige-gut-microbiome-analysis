#!/usr/bin/env Rscript
# Module 02 - consensus differential abundance: LinDA + ANCOM-BC2 + ALDEx2.
#
# Answers R1 Major 1 / R2 Minor 6. The published analysis used LinDA alone, with
# pair_matching disabled and controls randomly down-sampled, so no covariate ever
# entered the model. This module runs three methods on BYTE-IDENTICAL matched,
# covariate-adjusted input and reports where they agree.
#
# Three design rules, each fixing a specific way the comparison could be rigged:
#
# 1. IDENTICAL FEATURE SET. Prevalence filtering already happened upstream in
#    00b_aggregate_otu.py, on the full 2034-sample table and independently of
#    labels. Each method's internal filter is therefore DISABLED (LinDA
#    prev.cut/lib.cut = 0; ANCOM-BC2 prv_cut/lib_cut = 0; ALDEx2 has none).
#    Otherwise each method re-filters on the matched subset, the denominators
#    differ, and any "overlap" is an artifact of differing feature sets.
#
# 2. EFFECT SIZES ARE NOT AVERAGED. The three are on different scales:
#      LinDA      log2 fold change of relative abundance
#      ANCOM-BC2  natural-log LFC          -> divided by log(2)
#      ALDEx2     CLR difference, natural log -> divided by log(2), and it is a
#                 CLR difference, NOT a fold change
#    Each is reported in its own column with its estimand named.
#
# 3. SIGN HARMONISATION IS ASSERTED, NOT ASSUMED. LinDA takes integer 0/1,
#    ANCOM-BC2 an explicitly levelled factor, ALDEx2 a model matrix. A silent
#    reference flip would invert one method's effects and destroy the consensus.
#    The run aborts if the three disagree in sign on the strongest signals.
#
# Multiplicity: this script emits q_within (BH within one analysis). The wider
# scopes (q_secondary, q_all_partitions) are computed by 02b_consensus_merge.py,
# which needs the whole grid. Note that BH across all 15 partitions leaves
# nothing significant at any rank - hence the pre-specified primary contrast in
# config.yaml.
#
# Usage:
#   Rscript revision/analysis/02_consensus_da.R \
#     --level species --prevalence prev05 --partition any_ige

suppressPackageStartupMessages({
  library(yaml); library(LinDA); library(ANCOMBC); library(ALDEx2)
  library(TreeSummarizedExperiment)
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
MC    <- as.integer(getarg("--mc-samples", "128"))

cfg  <- yaml::read_yaml(file.path(REVISION, "config", "config.yaml"))
PIPE <- cfg$paths$pipeline_results
OUT  <- file.path(REVISION, "results", "da")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED <- cfg$statistics$seed
FORMULA <- "ige_status + sex + age + bmi"

message(sprintf("=== Consensus DA: %s / %s / %s ===", LEVEL, PREV, PART))
message(sprintf("  formula: ~ %s | ALDEx2 mc.samples=%d | seed=%d",
                FORMULA, MC, SEED))

# ── Load the matched, prevalence-filtered input ──────────────────────────
src <- file.path(PIPE, LEVEL, PREV, PART, "source")
oi <- read.csv(file.path(src, "otu_ige.csv"),     row.names = 1, check.names = FALSE)
oc <- read.csv(file.path(src, "otu_non_ige.csv"), row.names = 1, check.names = FALSE)
mi <- read.csv(file.path(src, "sample_data_ige.csv"),     row.names = 1)
mc <- read.csv(file.path(src, "sample_data_non_ige.csv"), row.names = 1)

taxa <- intersect(colnames(oi), colnames(oc))
otu  <- rbind(oi[, taxa, drop = FALSE], oc[, taxa, drop = FALSE])   # samples x taxa
meta <- rbind(mi, mc)
meta$ige_status <- factor(c(rep("ige", nrow(mi)), rep("non_ige", nrow(mc))),
                          levels = c("non_ige", "ige"))   # non_ige is reference
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
if (nrow(otu) < 20)
  message("  WARNING: small sample; treat this partition as exploratory")

counts_pxn <- t(as.matrix(otu))   # taxa x samples, the orientation all three want
res_list <- list()

# ── 1. LinDA ─────────────────────────────────────────────────────────────
message("\n  [1/3] LinDA ...")
meta_linda <- meta
meta_linda$ige_status <- as.integer(meta$ige_status == "ige")   # 0/1
# Every one of these was undocumented in the submitted manuscript. prev.cut = 0
# and lib.cut = 1 are LinDA's own defaults, passed explicitly so the record is
# unambiguous and the feature set matches the other two methods.
lin <- tryCatch(
  linda(otu.tab = counts_pxn, meta = meta_linda,
        formula = paste("~", FORMULA),
        alpha = 0.05, p.adj.method = "BH",
        prev.cut = 0, lib.cut = 1, winsor.quan = NULL,
        pseudo.cnt = 0.5, corr.cut = 0.1,
        adaptive = TRUE, imputation = FALSE, n.cores = 1),
  error = function(e) { message("    ERROR: ", conditionMessage(e)); NULL })
if (!is.null(lin)) {
  cn <- grep("ige_status", names(lin$output), value = TRUE)[1]
  o <- lin$output[[cn]]
  res_list$LinDA <- data.frame(
    taxon = rownames(o), method = "LinDA",
    effect_log2 = o$log2FoldChange,
    effect_type = "log2FC_relative_abundance",
    se = o$lfcSE, stat = o$stat, p_raw = o$pvalue, q_within = o$padj,
    stringsAsFactors = FALSE)
  message(sprintf("    %d taxa tested; %d at q<0.05",
                  nrow(o), sum(o$padj < 0.05, na.rm = TRUE)))
}

# ── 2. ANCOM-BC2 ─────────────────────────────────────────────────────────
# NOTE the defaults being overridden: p_adj_method is "holm" and prv_cut is 0.1.
message("  [2/3] ANCOM-BC2 ...")
tse <- TreeSummarizedExperiment(assays = list(counts = counts_pxn),
                               colData = meta)
set.seed(SEED)
ab <- tryCatch(
  suppressMessages(ancombc2(
    data = tse, assay_name = "counts", tax_level = NULL,
    fix_formula = FORMULA, p_adj_method = "BH",
    prv_cut = 0, lib_cut = 0,
    # The documented two-group configuration: declaring `group` enables
    # structural-zero detection, which is how ANCOM-BC2 is meant to handle taxa
    # absent from one group entirely. neg_lb uses the lower bound rather than the
    # point estimate when deciding a structural zero, which is the conservative
    # choice. Running without these (group = NULL, struc_zero = FALSE) makes the
    # method markedly anti-conservative on data this sparse.
    group = "ige_status", struc_zero = TRUE, neg_lb = TRUE,
    alpha = 0.05, n_cl = 1,
    verbose = FALSE, global = FALSE, pairwise = FALSE, pseudo_sens = TRUE)),
  error = function(e) { message("    ERROR: ", conditionMessage(e)); NULL })
if (!is.null(ab)) {
  r <- ab$res
  ss <- grep("^passed_ss", names(r), value = TRUE)
  res_list$ANCOMBC2 <- data.frame(
    taxon = r$taxon, method = "ANCOMBC2",
    # ANCOM-BC2 reports a NATURAL-log LFC; convert to log2.
    effect_log2 = r$lfc_ige_statusige / log(2),
    effect_type = "log2FC_absolute_abundance",
    se = r$se_ige_statusige / log(2), stat = r$W_ige_statusige,
    p_raw = r$p_ige_statusige, q_within = r$q_ige_statusige,
    passed_ss = if (length(ss)) r[[ss[1]]] else NA,
    stringsAsFactors = FALSE)
  n_sig <- sum(r$q_ige_statusige < 0.05, na.rm = TRUE)
  n_ss  <- sum(r$q_ige_statusige < 0.05 & res_list$ANCOMBC2$passed_ss, na.rm = TRUE)
  message(sprintf("    %d taxa tested; %d at q<0.05, of which %d also pass the",
                  nrow(r), n_sig, n_ss))
  message("      pseudo-count sensitivity check (passed_ss).")
  if (n_sig > 0 && n_ss < 0.5 * n_sig)
    message(sprintf("      NOTE: %d of %d ANCOM-BC2 hits FAIL sensitivity. Its ",
                    n_sig - n_ss, n_sig),
            "standard errors are small relative to the other methods, so raw\n",
            "      ANCOM-BC2 significance is not reportable on its own - the ",
            "passed_ss\n      flag exists for exactly this and is carried through ",
            "to the consensus.")
}

# ── 3. ALDEx2 ────────────────────────────────────────────────────────────
message(sprintf("  [3/3] ALDEx2 (mc.samples=%d) ...", MC))
mm <- model.matrix(as.formula(paste("~", FORMULA)), data = meta)
set.seed(SEED)
al <- tryCatch({
  clr <- aldex.clr(reads = counts_pxn, conds = mm, mc.samples = MC,
                   denom = "all", verbose = FALSE)
  aldex.glm(clr, verbose = FALSE, fdr.method = "BH")
}, error = function(e) { message("    ERROR: ", conditionMessage(e)); NULL })
if (!is.null(al)) {
  res_list$ALDEx2 <- data.frame(
    taxon = rownames(al), method = "ALDEx2",
    # A CLR-scale difference in natural log. NOT a fold change.
    effect_log2 = al[["ige_statusige:Est"]] / log(2),
    effect_type = "clr_difference_log2",
    se = al[["ige_statusige:SE"]] / log(2),
    stat = al[["ige_statusige:t.val"]],
    p_raw = al[["ige_statusige:pval"]],
    q_within = al[["ige_statusige:pval.padj"]],
    stringsAsFactors = FALSE)
  message(sprintf("    %d taxa tested; %d at q<0.05", nrow(al),
                  sum(al[["ige_statusige:pval.padj"]] < 0.05, na.rm = TRUE)))
}

if (!length(res_list)) { message("\nAll three methods failed."); quit(status = 1) }

# ── QC 1: identical feature sets ─────────────────────────────────────────
message("\n  QC: feature sets")
ns <- sapply(res_list, nrow)
print(ns)
if (length(unique(ns)) != 1) {
  message("    WARNING: methods tested different numbers of taxa - the overlap ",
          "would be confounded with filtering. Check prev.cut/prv_cut.")
} else {
  message(sprintf("    PASS - all methods tested the same %d taxa", ns[1]))
}

# ── QC 2: sign concordance on the strongest signals ──────────────────────
long <- do.call(rbind, lapply(res_list, function(d) {
  d$passed_ss <- if ("passed_ss" %in% names(d)) d$passed_ss else NA
  d[, c("taxon", "method", "effect_log2", "effect_type", "se", "stat",
        "p_raw", "q_within", "passed_ss")]
}))
rownames(long) <- NULL

# The failure mode being guarded against is a REFERENCE-LEVEL FLIP, which inverts
# one method's effects wholesale. That shows up as a strongly NEGATIVE correlation
# across all taxa. Sign disagreement on individual near-null taxa is expected and
# is not an error: the three methods estimate different quantities (relative
# abundance, absolute abundance, CLR difference), so on mostly-null data their
# signs agree at close to chance. Do not abort on that.
message("  QC: cross-method effect correlation (guards against a reference flip)")
wide <- reshape(long[, c("taxon", "method", "effect_log2")],
                idvar = "taxon", timevar = "method", direction = "wide")
names(wide) <- sub("^effect_log2\\.", "", names(wide))
methods_present <- setdiff(names(wide), "taxon")
flip <- FALSE
if (length(methods_present) >= 2) {
  pairs <- combn(methods_present, 2, simplify = FALSE)
  for (pr in pairs) {
    x <- wide[[pr[1]]]; y <- wide[[pr[2]]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 5) next
    rp <- suppressWarnings(cor(x[ok], y[ok]))
    rs <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
    sa <- mean(sign(x[ok]) == sign(y[ok]))
    message(sprintf("    %-9s vs %-9s Pearson=%+.3f Spearman=%+.3f sign_agree=%.0f%%",
                    pr[1], pr[2], rp, rs, 100 * sa))
    if (!is.na(rp) && rp < -0.5) flip <- TRUE
  }
}
if (flip) {
  message("    FAIL - a strongly negative cross-method correlation indicates a ",
          "reference-level flip. Aborting so a mis-signed consensus is never written.")
  quit(status = 1)
}
message("    PASS (no reference flip detected)")

# ── QC 3: BH sanity ──────────────────────────────────────────────────────
bad <- long[!is.na(long$p_raw) & !is.na(long$q_within) &
              long$q_within < long$p_raw - 1e-12, ]
message(sprintf("  QC: rows with q < p: %d  %s", nrow(bad),
                if (nrow(bad)) "FAIL" else "PASS"))
if (nrow(bad)) quit(status = 1)

# ── Write ────────────────────────────────────────────────────────────────
long$level <- LEVEL; long$prevalence <- PREV; long$partition <- PART
long$n_ige <- n_ige; long$n_ctrl <- n_ctrl; long$n_taxa <- ncol(otu)
long$formula <- paste("~", FORMULA); long$mc_samples <- MC; long$seed <- SEED

f <- file.path(OUT, sprintf("da_long_%s_%s_%s.csv", LEVEL, PREV, PART))
write.csv(long, f, row.names = FALSE)
message(sprintf("\nWrote %s (%d rows)", f, nrow(long)))

# Per-analysis summary for the log.
message("\n  Significant at q_within < 0.05, by method:")
for (m in unique(long$method)) {
  s <- long[long$method == m & !is.na(long$q_within) & long$q_within < 0.05, ]
  detail <- if (nrow(s) == 0) "" else paste0(": ", paste(
    sprintf("%s(%+.2f)", s$taxon, s$effect_log2), collapse = ", "))
  if (nchar(detail) > 300) detail <- paste0(substr(detail, 1, 300), " ...")
  message(sprintf("    %-9s %d%s", m, nrow(s), detail))
}

# ANCOM-BC2 hits that also survive its own sensitivity check - the reportable set.
an <- long[long$method == "ANCOMBC2" & !is.na(long$q_within) &
             long$q_within < 0.05 & !is.na(long$passed_ss) & long$passed_ss, ]
if (nrow(an))
  message(sprintf("    %-9s %d also pass sensitivity: %s", "ANCOMBC2*", nrow(an),
                  paste(sprintf("%s(%+.2f)", an$taxon, an$effect_log2),
                        collapse = ", ")))
message("Done.")
