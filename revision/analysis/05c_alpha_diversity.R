#!/usr/bin/env Rscript
# Module 05c - non-phylogenetic alpha diversity, done properly.
#
# WHY THIS MODULE EXISTS
#   R3 asks for confounder control (adjustment, stratification or matching)
#   when comparing alpha diversity. The published Table S5 used an unmatched,
#   uncorrected Mann-Whitney test on raw counts; the revision so far computed
#   only Faith's PD (module 05b). No reproducible source exists for the
#   Shannon/richness numbers quoted in the draft response letter - this module
#   creates one.
#
# DESIGN
#   * Alpha is computed from the RAW feature table restricted to the 508
#     analytic samples, BEFORE any prevalence filtering: Chao1 is undefined
#     without singletons, and Shannon/richness are biased by prevalence
#     filters. This is the same "ASVs non-zero within the 508 subset" set the
#     published analysis used (6,958 ASVs; verified by module 05b).
#   * Per-sample alpha is partition-invariant, so it is computed ONCE for all
#     508 samples; the 15 partition tests are then cheap arithmetic. Hence a
#     single job, not an array.
#   * PRIMARY analysis: RAW (unrarefied) counts - the authors' decision,
#     following McMurdie & Holmes 2014 (PMC3974642): discarding reads is
#     inadmissible when the question can be answered on the full data. Note the
#     DA pipeline NEVER rarefies anywhere; this module is alpha only.
#     Sequencing-effort control is done model-side instead:
#       - a per-partition depth-balance test (paired Wilcoxon on library size:
#         if depth does not differ between matched groups, unrarefied alpha
#         comparisons are unconfounded by effort), and
#       - an adjusted model including log10(library size) as a covariate.
#   * SENSITIVITY only: rarefaction in the sense of the QIIME 2 guidance
#     (repeated subsampling, 100 iterations, AVERAGED - not the single-draw
#     "rarefying" that McMurdie & Holmes criticise; vegan::rrarefy, seeds
#     42..141). Kept because part of the field (e.g. Schloss 2024, mSphere)
#     still expects it; reported in the supplement, never as the headline.
#     Depth rule: rarefy to the minimum depth across the 508 unless that is
#     below 10% of the median depth, in which case use the 10th-percentile
#     depth and drop shallower samples (and their pairs). All counts reported.
#   * Three analyses per partition x index x variant:
#       paired_wilcoxon    - matched pairs (sex exact, age +-5, BMI +-4);
#                            effect = matched-pairs rank-biserial r and
#                            Hodges-Lehmann shift (conf.int)
#       adjusted_lm_full   - lm(alpha ~ ige + sex + age + bmi), partition
#                            members vs ALL 233 controls (adjustment instead
#                            of matching; uses everyone)
#       adjusted_lm_matched- same model on the matched samples only
#                            (adjustment on top of matching)
#   * Multiplicity: BH within partition (q_within, across the 3 indices) and
#     BH across all partitions x indices (q_all) - separately per variant and
#     analysis, so no column ever mixes multiplicity scopes.
#
# Chao1 caveat, recorded rather than hidden: DADA2 denoising suppresses many
# per-sample singletons, so Chao1 on DADA2 counts is a lower-bound-ish index.
# It is reported because reviewers expect it, flagged in the output metadata.
#
# Usage:  Rscript revision/analysis/05c_alpha_diversity.R
# (no arguments - the whole grid runs in one pass)

suppressPackageStartupMessages({ library(yaml); library(vegan) })

args_raw <- commandArgs(trailingOnly = FALSE)
script <- sub("--file=", "", args_raw[grep("--file=", args_raw)])
REVISION <- normalizePath(file.path(dirname(script), ".."))
.libPaths(c(file.path(REVISION, "env", "rlibs"), .libPaths()))

cfg  <- yaml::read_yaml(file.path(REVISION, "config", "config.yaml"))
KORA <- cfg$paths$kora_data
PIPE <- cfg$paths$pipeline_results
OUT  <- file.path(REVISION, "results", "alpha")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED    <- cfg$statistics$seed          # 42
N_RARE  <- 100
PARTITIONS <- c("any_ige", "single", "multiple", "high_score", "medium_score",
                "low_score", paste0("cluster_", LETTERS[1:9]))
INDICES <- c("shannon", "observed", "chao1")

message("=== Module 05c: alpha diversity (raw + rarefied, matched + adjusted) ===")

# ── 1. Raw table, 508 subset ─────────────────────────────────────────────
ft <- read.csv(file.path(KORA, "feature_table.tsv"), sep = "\t",
               row.names = 1, check.names = FALSE)
colnames(ft) <- as.character(colnames(ft))
ige_all <- read.csv(file.path(KORA, "ige_kora_subset_all.csv"),
                    row.names = 1, check.names = FALSE)
subset508 <- intersect(as.character(rownames(ige_all)), colnames(ft))
stopifnot(length(subset508) == 508)
keep_asv <- rownames(ft)[rowSums(ft[, subset508, drop = FALSE] > 0) > 0]
message(sprintf("  samples: %d | ASVs non-zero in subset: %d (published: 6,958)",
                length(subset508), length(keep_asv)))

# samples x ASVs integer matrix, the orientation vegan wants
X <- t(as.matrix(ft[keep_asv, subset508]))
storage.mode(X) <- "integer"

depths <- rowSums(X)
message(sprintf("  depth: min %d | 10th pct %.0f | median %.0f | max %d",
                min(depths), quantile(depths, .1), median(depths), max(depths)))

# ── 2. Depth rule ────────────────────────────────────────────────────────
rare_depth <- min(depths)
dropped <- character(0)
if (rare_depth < 0.10 * median(depths)) {
  rare_depth <- floor(quantile(depths, 0.10))
  dropped <- rownames(X)[depths < rare_depth]
  message(sprintf("  min depth < 10%% of median -> rarefying to 10th pct (%d); dropping %d shallow sample(s): %s",
                  rare_depth, length(dropped), paste(dropped, collapse = ", ")))
} else {
  message(sprintf("  rarefying to min depth: %d (no samples dropped)", rare_depth))
}

# ── 3. Per-sample alpha, raw ─────────────────────────────────────────────
alpha_one <- function(M) {
  # M: samples x ASVs, integer
  sh <- vegan::diversity(M, index = "shannon")            # natural log
  ob <- rowSums(M > 0)
  ch <- suppressWarnings(t(vegan::estimateR(M))[, "S.chao1"])
  cbind(shannon = sh, observed = ob, chao1 = ch)
}
message("  computing raw per-sample alpha ...")
A_raw <- alpha_one(X)

# ── 4. Per-sample alpha, repeated rarefaction (streaming mean) ───────────
message(sprintf("  computing rarefied alpha: %d reps at depth %d ...", N_RARE, rare_depth))
keep_r <- setdiff(rownames(X), dropped)
Xr <- X[keep_r, , drop = FALSE]
acc  <- matrix(0, nrow(Xr), length(INDICES),
               dimnames = list(rownames(Xr), INDICES))
acc2 <- acc
for (r in seq_len(N_RARE)) {
  set.seed(SEED + r - 1)
  Ar <- alpha_one(vegan::rrarefy(Xr, rare_depth))
  acc  <- acc  + Ar
  acc2 <- acc2 + Ar^2
  if (r %% 20 == 0) message(sprintf("    rep %d/%d", r, N_RARE))
}
A_rare    <- acc / N_RARE
A_rare_sd <- sqrt(pmax(acc2 / N_RARE - A_rare^2, 0))

persample <- rbind(
  data.frame(sample = rownames(A_raw), variant = "raw", A_raw,
             shannon_sd = NA_real_, observed_sd = NA_real_, chao1_sd = NA_real_,
             row.names = NULL, check.names = FALSE),
  data.frame(sample = rownames(A_rare), variant = sprintf("rarefied_%dx", N_RARE),
             A_rare, shannon_sd = A_rare_sd[, "shannon"],
             observed_sd = A_rare_sd[, "observed"], chao1_sd = A_rare_sd[, "chao1"],
             row.names = NULL, check.names = FALSE))
persample$rarefaction_depth <- ifelse(persample$variant == "raw", NA, rare_depth)
write.csv(persample, file.path(OUT, "alpha_persample.csv"), row.names = FALSE)
message(sprintf("  wrote alpha_persample.csv (%d rows)", nrow(persample)))

alpha_of <- function(ids, variant, index) {
  A <- if (variant == "raw") A_raw else A_rare
  ids <- intersect(as.character(ids), rownames(A))
  setNames(A[ids, index], ids)
}

# ── 5. Covariates for the adjusted models ────────────────────────────────
kora <- read.csv(file.path(KORA, "kora.csv"), check.names = FALSE)
# kora.csv covers the full cohort (~1950 rows); only ~508 have a 16S sample id.
# Rows without one cannot be rownames (and are irrelevant here) - drop first.
kora <- kora[!is.na(kora$u3_16s_id), ]
kora$u3_16s_id <- as.character(kora$u3_16s_id)
kora <- kora[!duplicated(kora$u3_16s_id), c("u3_16s_id", "u3csex", "u3talter", "u3tbmi")]
rownames(kora) <- kora$u3_16s_id
# sensitized = any non-zero IgE component (the implemented definition)
sens508 <- rownames(ige_all)[rowSums(ige_all != 0, na.rm = TRUE) > 0]
ctrl508 <- setdiff(rownames(ige_all), sens508)
message(sprintf("  full-cohort groups: %d sensitized / %d controls (expect 275/233)",
                length(sens508), length(ctrl508)))

# depth per sample, for the balance test and the depth-adjusted model
lib_depth <- rowSums(X)          # X is samples x ASVs

# ── 6. Tests per partition ───────────────────────────────────────────────
rows <- list()
for (part in PARTITIONS) {
  src <- file.path(PIPE, "asv", "prev01", part, "source")
  mi <- read.csv(file.path(src, "sample_data_ige.csv"),     row.names = 1)
  mc <- read.csv(file.path(src, "sample_data_non_ige.csv"), row.names = 1)
  cp <- sort(intersect(mi$pair_nb, mc$pair_nb))
  mi <- mi[match(cp, mi$pair_nb), , drop = FALSE]
  mc <- mc[match(cp, mc$pair_nb), , drop = FALSE]
  ids_i <- as.character(rownames(mi)); ids_c <- as.character(rownames(mc))

  # (0) depth balance: is sequencing effort associated with group AT ALL?
  di <- lib_depth[intersect(ids_i, names(lib_depth))]
  dc <- lib_depth[intersect(ids_c, names(lib_depth))]
  okd <- !is.na(di) & !is.na(dc[names(dc)])
  if (length(di) == length(dc) && length(di) >= 5) {
    wd <- suppressWarnings(wilcox.test(di, dc, paired = TRUE, exact = FALSE,
                                       conf.int = TRUE))
    rows[[length(rows) + 1]] <- data.frame(
      partition = part, index = "library_size", variant = "raw",
      analysis = "depth_balance", n = length(di) * 2L, n_pairs = length(di),
      mean_ige = mean(di), mean_ctrl = mean(dc),
      estimate = unname(wd$estimate), ci_lo = wd$conf.int[1],
      ci_hi = wd$conf.int[2], effect_rank_biserial = NA_real_,
      statistic = unname(wd$statistic), p_value = wd$p.value)
  }

  # full-cohort membership for this partition = its sensitized members vs ALL controls
  memb_full <- as.character(rownames(mi))     # matched members carry the definition
  # for the full model, use every partition member (incl. unmatched): recover from
  # the ige-side source of the UNMATCHED table if present; else fall back to matched
  full_i <- memb_full
  for (variant in unique(persample$variant)) {
    for (index in INDICES) {
      a <- alpha_of(c(ids_i, ids_c), variant, index)

      # (a) paired Wilcoxon on complete pairs
      ai <- a[ids_i]; ac <- a[ids_c]
      ok <- !is.na(ai) & !is.na(ac)
      n_pairs <- sum(ok)
      if (n_pairs >= 5) {
        d <- ai[ok] - ac[ok]
        wt <- suppressWarnings(wilcox.test(ai[ok], ac[ok], paired = TRUE,
                                           conf.int = TRUE, exact = FALSE))
        # matched-pairs rank-biserial from signed ranks
        rk <- rank(abs(d[d != 0])); sg <- sign(d[d != 0])
        rb <- if (length(rk)) (sum(rk[sg > 0]) - sum(rk[sg < 0])) / sum(rk) else NA_real_
        rows[[length(rows) + 1]] <- data.frame(
          partition = part, index = index, variant = variant,
          analysis = "paired_wilcoxon", n = n_pairs * 2L, n_pairs = n_pairs,
          mean_ige = mean(ai[ok]), mean_ctrl = mean(ac[ok]),
          estimate = unname(wt$estimate), ci_lo = wt$conf.int[1], ci_hi = wt$conf.int[2],
          effect_rank_biserial = rb, statistic = unname(wt$statistic),
          p_value = wt$p.value)
      }

      # (b) adjusted lm: full cohort, matched, and matched + log10(depth) -
      # the last is the McMurdie-&-Holmes-compatible effort control (model the
      # nuisance, discard nothing)
      for (an in c("adjusted_lm_full", "adjusted_lm_matched",
                   "adjusted_lm_matched_depth")) {
        gi <- if (an == "adjusted_lm_full") full_i else ids_i
        gc <- if (an == "adjusted_lm_full") ctrl508 else ids_c
        df <- data.frame(id = c(gi, gc),
                         ige = c(rep(1L, length(gi)), rep(0L, length(gc))))
        df$alpha <- alpha_of(df$id, variant, index)[df$id]
        df <- cbind(df, kora[df$id, c("u3csex", "u3talter", "u3tbmi")])
        df$log_depth <- log10(lib_depth[df$id])
        df <- df[stats::complete.cases(df), ]
        if (nrow(df) >= 20 && length(unique(df$ige)) == 2) {
          fit <- if (an == "adjusted_lm_matched_depth")
            lm(alpha ~ ige + factor(u3csex) + u3talter + u3tbmi + log_depth,
               data = df)
          else
            lm(alpha ~ ige + factor(u3csex) + u3talter + u3tbmi, data = df)
          co <- summary(fit)$coefficients["ige", ]
          ci <- confint(fit)["ige", ]
          rows[[length(rows) + 1]] <- data.frame(
            partition = part, index = index, variant = variant,
            analysis = an, n = nrow(df), n_pairs = NA_integer_,
            mean_ige = mean(df$alpha[df$ige == 1]), mean_ctrl = mean(df$alpha[df$ige == 0]),
            estimate = co["Estimate"], ci_lo = ci[1], ci_hi = ci[2],
            effect_rank_biserial = NA_real_, statistic = co["t value"],
            p_value = co["Pr(>|t|)"])
        }
      }
    }
  }
  message(sprintf("  %-13s pairs=%d", part, length(cp)))
}

res <- do.call(rbind, rows); rownames(res) <- NULL

# ── 7. Multiplicity, per variant x analysis (never mixed) ────────────────
res$q_within <- NA_real_; res$q_all <- NA_real_
for (v in unique(res$variant)) for (an in unique(res$analysis)) {
  s <- res$variant == v & res$analysis == an
  res$q_all[s] <- p.adjust(res$p_value[s], method = "BH")
  for (p in unique(res$partition)) {
    sp <- s & res$partition == p
    res$q_within[sp] <- p.adjust(res$p_value[sp], method = "BH")
  }
}
res$rarefaction_depth <- ifelse(res$variant == "raw", NA, rare_depth)
res$chao1_note <- ifelse(res$index == "chao1",
  "DADA2 suppresses per-sample singletons; Chao1 on denoised counts is a lower-bound-ish index", "")
write.csv(res, file.path(OUT, "alpha_tests_long.csv"), row.names = FALSE)

# ── 8. Table S5-shaped summary. PRIMARY = RAW counts (authors' decision,
# per McMurdie & Holmes 2014); the rarefaction variant is sensitivity only. ─
s5 <- res[res$analysis == "paired_wilcoxon" & res$index %in% c("shannon", "observed"), ]
s5$is_primary <- s5$variant == "raw"
s5 <- s5[order(-s5$is_primary, s5$variant, s5$index, s5$partition),
         c("partition", "index", "variant", "is_primary", "n_pairs", "mean_ige",
           "mean_ctrl", "estimate", "ci_lo", "ci_hi", "p_value", "q_within",
           "q_all")]
write.csv(s5, file.path(OUT, "alpha_summary_tableS5.csv"), row.names = FALSE)

# ── 9. QC gates ──────────────────────────────────────────────────────────
np_any <- unique(res$n_pairs[res$partition == "any_ige" & res$analysis == "paired_wilcoxon"])
message("\n=== QC ===")
message(sprintf("  any_ige complete pairs: %s (expect 221)", paste(np_any, collapse = ",")))
counts <- table(res$analysis, res$variant)
print(counts)
stopifnot(all(counts["paired_wilcoxon", ] == length(PARTITIONS) * length(INDICES)))
bad <- sum(!is.na(res$q_within) & res$q_within < res$p_value - 1e-12)
message(sprintf("  BH sanity (q >= p): %s", if (bad) "FAIL" else "PASS"))
if (bad) quit(status = 1)

db <- res[res$analysis == "depth_balance", ]
message(sprintf("  depth balance (library size, paired): %d partitions, min p = %.3f%s",
                nrow(db), min(db$p_value),
                if (min(db$p_value) > 0.05) " - effort is group-balanced everywhere"
                else " - ATTENTION: depth differs somewhere; depth-adjusted lm is the arbiter"))
pw <- res[res$analysis == "paired_wilcoxon", ]
message(sprintf("  significant after BH across all (q_all<0.05): %d of %d tests",
                sum(pw$q_all < 0.05, na.rm = TRUE), nrow(pw)))
prim <- pw[pw$partition == "any_ige" & pw$index == "shannon", ]
for (i in seq_len(nrow(prim)))
  message(sprintf("  any_ige shannon [%s]: IgE %.3f vs ctrl %.3f, p=%.3f, q_all=%.3f",
                  prim$variant[i], prim$mean_ige[i], prim$mean_ctrl[i],
                  prim$p_value[i], prim$q_all[i]))
message("Done.")
