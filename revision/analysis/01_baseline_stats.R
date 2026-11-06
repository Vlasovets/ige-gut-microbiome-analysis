#!/usr/bin/env Rscript
# Tables 1 and 2 with formal statistical comparisons.
#
# The submitted manuscript reports Tables 1 and 2 with NO p-values, no test
# statistic and no named test, yet draws comparative conclusions from them
# ("These differences confirm that IgE sensitization is associated with
# increased burden of allergic disease", L257-258). This module supplies the
# tests, effect sizes with confidence intervals, and BH correction within each
# table.
#
# Test choice:
#   categorical  -> chi-squared, with Fisher's exact when any expected cell < 5
#                   (mandatory here: asthma is 6/233 in controls, PPI is ~16/233)
#   continuous   -> Wilcoxon rank-sum for 2 groups; Kruskal-Wallis across the
#                   9 clusters, with Dunn post-hoc when significant
#   effect sizes -> Cliff's delta (continuous), Cramer's V and risk difference
#                   with Wilson CI (categorical)
#
# Two versions of Table 1 are produced:
#   (i)  full 508 cohort  - must reproduce the published counts exactly
#   (ii) matched 442 subset - a QC assertion that matching worked, NOT a result;
#        age/sex/BMI should be non-significant here by construction
#
# Usage:
#   Rscript revision/analysis/01_baseline_stats.R

suppressPackageStartupMessages({
  library(yaml)
})

REVISION <- file.path(dirname(sub("--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])), "..")
REVISION <- normalizePath(REVISION)
.libPaths(c(file.path(REVISION, "env", "rlibs"), .libPaths()))

cfg <- yaml::read_yaml(file.path(REVISION, "config", "config.yaml"))
KORA <- cfg$paths$kora_data
PIPE <- cfg$paths$pipeline_results
OUT  <- file.path(REVISION, "results", "tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

set.seed(cfg$statistics$seed)

# ── Variable dictionary ───────────────────────────────────────────────────
# Verified present in kora.csv with zero missingness in the matched set.
VARS <- list(
  age        = list(col = "u3talter",      label = "Age (years)",        type = "cont"),
  bmi        = list(col = "u3tbmi",        label = "BMI (kg/m2)",        type = "cont"),
  sex        = list(col = "u3csex",        label = "Sex",                type = "cat"),
  education  = list(col = "u3tschul",      label = "Education",          type = "cat"),
  smoking    = list(col = "u3tcigreg_sf",  label = "Smoking",            type = "cat"),
  packyears  = list(col = "u3tpyrs",       label = "Pack-years",         type = "cont"),
  alcohol    = list(col = "u3talkkon",     label = "Alcohol (g/day)",    type = "cont"),
  activity   = list(col = "u3tphact",      label = "Physical activity",  type = "cat"),
  ppi        = list(col = "u3tmppi",       label = "PPI use",            type = "cat"),
  asthma     = list(col = "u3da13aa",      label = "Asthma",             type = "cat"),
  dermatitis = list(col = "u3da4aa",       label = "Dermatitis",         type = "cat"),
  hayfever   = list(col = "u3da8aa",       label = "Hay fever",          type = "cat")
)
# Colliders: downstream of sensitization. Tabulated as outcomes, never used as
# adjustment covariates elsewhere in the revision.
COLLIDERS <- c("asthma", "dermatitis", "hayfever")

# ── Effect sizes ──────────────────────────────────────────────────────────
cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  # P(X>Y) - P(X<Y), computed from the Mann-Whitney U
  u <- sum(outer(x, y, ">")) + 0.5 * sum(outer(x, y, "=="))
  2 * u / (length(x) * length(y)) - 1
}

cramers_v <- function(tab) {
  if (any(dim(tab) < 2)) return(NA_real_)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  sqrt(as.numeric(chi) / (sum(tab) * (min(dim(tab)) - 1)))
}

wilson_ci <- function(k, n, conf = 0.95) {
  if (n == 0) return(c(NA, NA))
  z <- qnorm(1 - (1 - conf) / 2); p <- k / n
  d <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d)
}

# ── Two-group comparison ──────────────────────────────────────────────────
compare_two <- function(v, grp, spec) {
  a <- v[grp]; b <- v[!grp]
  if (spec$type == "cont") {
    a <- as.numeric(a); b <- as.numeric(b)
    ok <- sum(!is.na(a)) > 2 && sum(!is.na(b)) > 2
    if (!ok) return(NULL)
    # Shapiro on the pooled sample decides how we *describe*, not how we test:
    # we test with Wilcoxon throughout for consistency with the paired design.
    sw <- tryCatch(shapiro.test(c(a, b)[seq_len(min(5000, length(c(a, b))))])$p.value,
                   error = function(e) NA_real_)
    tt <- suppressWarnings(wilcox.test(a, b, exact = FALSE, conf.int = TRUE))
    list(
      summary_1 = sprintf("%.2f (%.2f)", mean(a, na.rm = TRUE), sd(a, na.rm = TRUE)),
      summary_0 = sprintf("%.2f (%.2f)", mean(b, na.rm = TRUE), sd(b, na.rm = TRUE)),
      test = "Wilcoxon rank-sum", statistic = as.numeric(tt$statistic),
      p = tt$p.value, effect_name = "Cliff's delta",
      effect = cliffs_delta(a, b),
      ci_low = NA_real_, ci_high = NA_real_,
      note = sprintf("Shapiro p=%.3g", sw)
    )
  } else {
    tab <- table(factor(v[!is.na(v)]), grp[!is.na(v)])
    if (any(dim(tab) < 2)) return(NULL)
    ex <- suppressWarnings(chisq.test(tab, correct = FALSE)$expected)
    use_fisher <- any(ex < 5)
    tt <- if (use_fisher) {
      tryCatch(fisher.test(tab, simulate.p.value = nrow(tab) > 2, B = 1e5),
               error = function(e) list(p.value = NA_real_))
    } else suppressWarnings(chisq.test(tab, correct = FALSE))
    # Risk difference for the positive level when binary.
    rd <- ci <- NA_real_; rdci <- c(NA, NA)
    if (nrow(tab) == 2) {
      k1 <- tab[2, "TRUE"]; n1 <- sum(tab[, "TRUE"])
      k0 <- tab[2, "FALSE"]; n0 <- sum(tab[, "FALSE"])
      rd <- k1 / n1 - k0 / n0
      # Newcombe (1998) hybrid-score interval, method 10: combine the two
      # Wilson intervals in quadrature. The previous version subtracted the
      # endpoints (lower1 - upper0, upper1 - lower0), which is not an interval
      # for a difference - it is ~sqrt(2) too wide and was mislabelled as a CI
      # for Cramer's V in the output columns.
      p1 <- k1 / n1; p0 <- k0 / n0
      c1 <- wilson_ci(k1, n1); c0 <- wilson_ci(k0, n0)
      rdci <- c(rd - sqrt((p1 - c1[1])^2 + (c0[2] - p0)^2),
                rd + sqrt((c1[2] - p1)^2 + (p0 - c0[1])^2))
    }
    fmt <- function(col) paste(sprintf("%d (%.1f%%)", tab[, col],
                                       100 * tab[, col] / sum(tab[, col])),
                               collapse = " / ")
    list(
      summary_1 = fmt("TRUE"), summary_0 = fmt("FALSE"),
      test = if (use_fisher) "Fisher's exact" else "Chi-squared",
      statistic = if (use_fisher) NA_real_ else as.numeric(tt$statistic),
      p = tt$p.value, effect_name = "Cramer's V", effect = cramers_v(tab),
      ci_low = NA_real_, ci_high = NA_real_,       # no CI is computed for V
      rd = rd, rd_ci_low = rdci[1], rd_ci_high = rdci[2],
      rd_level = if (nrow(tab) == 2) rownames(tab)[2] else NA_character_,
      note = sprintf("levels: %s%s", paste(rownames(tab), collapse = ","),
                     if (!is.na(rd)) sprintf("; risk diff=%.3f", rd) else "")
    )
  }
}

# ── Load and harmonise metadata ───────────────────────────────────────────
message("=== Loading data ===")
ige <- read.csv(file.path(KORA, "ige_kora_subset_all.csv"),
                row.names = 1, check.names = FALSE)
sens_all <- !apply(ige == 0, 1, all)
message(sprintf("  IgE panel: %d individuals x %d components (NOT 112)",
                nrow(ige), ncol(ige)))
message(sprintf("  Sensitized (any non-zero): %d / non-sensitized: %d",
                sum(sens_all), sum(!sens_all)))

meta <- read.csv(file.path(KORA, "kora.csv"), check.names = FALSE)
meta <- meta[!is.na(meta$u3_16s_id), ]
rownames(meta) <- as.character(as.integer(meta$u3_16s_id))
# -999 is a skip-logic sentinel throughout the KORA survey blocks.
meta[meta == -999] <- NA

# Recode to readable levels. Physical activity 1&2 collapse to "regular",
# which is what reproduces the published Table 1 counts.
recode <- function(m) {
  m$u3csex       <- factor(m$u3csex, levels = c(1, 2), labels = c("Male", "Female"))
  m$u3tschul     <- factor(m$u3tschul)
  m$u3tcigreg_sf <- factor(m$u3tcigreg_sf, levels = c(1, 3, 4),
                           labels = c("Regular", "Ex-smoker", "Never"))
  m$u3tphact     <- factor(ifelse(m$u3tphact %in% c(1, 2), "Regular",
                            ifelse(m$u3tphact == 3, "Irregular",
                              ifelse(m$u3tphact == 4, "None", NA))),
                           levels = c("Regular", "Irregular", "None"))
  for (cc in c("u3tmppi", "u3da13aa", "u3da4aa", "u3da8aa"))
    if (cc %in% names(m))
      m[[cc]] <- factor(ifelse(m[[cc]] == 1, "Yes", "No"), levels = c("No", "Yes"))
  m
}
meta <- recode(meta)

# ── Build a table for one cohort definition ──────────────────────────────
build_table <- function(ids, grp_logical, tag) {
  m <- meta[ids, , drop = FALSE]
  g <- grp_logical
  rows <- list()
  for (nm in names(VARS)) {
    spec <- VARS[[nm]]
    if (!spec$col %in% names(m)) next
    r <- compare_two(m[[spec$col]], g, spec)
    if (is.null(r)) next
    rows[[nm]] <- data.frame(
      table = tag, variable = spec$label, column = spec$col,
      type = spec$type, is_collider = nm %in% COLLIDERS,
      n_sensitized = sum(g), n_non_sensitized = sum(!g),
      sensitized = r$summary_1, non_sensitized = r$summary_0,
      test = r$test, statistic = r$statistic, p_value = r$p,
      effect_name = r$effect_name, effect_size = r$effect,
      effect_ci_low = r$ci_low, effect_ci_high = r$ci_high,
      risk_diff = if (is.null(r$rd)) NA_real_ else r$rd,
      risk_diff_ci_low = if (is.null(r$rd_ci_low)) NA_real_ else r$rd_ci_low,
      risk_diff_ci_high = if (is.null(r$rd_ci_high)) NA_real_ else r$rd_ci_high,
      risk_diff_level = if (is.null(r$rd_level)) NA_character_ else r$rd_level,
      note = r$note, stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  # BH WITHIN this table, across the variables tested.
  out$q_BH <- p.adjust(out$p_value, method = "BH")
  rownames(out) <- NULL
  out
}

# ── Table 1 (i): full 508 cohort — must reproduce published counts ────────
message("\n=== Table 1 (i): full cohort ===")
ids508 <- intersect(rownames(meta), names(sens_all)[!is.na(sens_all)])
g508 <- sens_all[ids508]
t1_full <- build_table(ids508, g508, "Table1_full508")
print(t1_full[, c("variable", "sensitized", "non_sensitized", "test",
                  "p_value", "q_BH")], row.names = FALSE, digits = 3)

# Validation against the published table.
message("\n  Validation against published counts:")
checks <- list(c("Asthma", "39", "6"), c("Dermatitis", "29", "12"),
               c("Hay fever", "129", "8"))
n_fail <- 0
for (ck in checks) {
  row <- t1_full[t1_full$variable == ck[1], ]
  if (!nrow(row)) next
  got_1 <- sub(".*/ *(\\d+) .*", "\\1", row$sensitized)
  got_0 <- sub(".*/ *(\\d+) .*", "\\1", row$non_sensitized)
  ok <- got_1 == ck[2] && got_0 == ck[3]
  if (!ok) n_fail <- n_fail + 1
  message(sprintf("    %-12s published %s/%s  recomputed %s/%s  %s",
                  ck[1], ck[2], ck[3], got_1, got_0, if (ok) "PASS" else "MISMATCH"))
}

# ── Table 1 (ii): matched 442 subset — a QC assertion ────────────────────
message("\n=== Table 1 (ii): matched subset (QC, not a result) ===")
ige_ids  <- readLines(file.path(PIPE, "matched", "any_ige", "ige_ids.txt"))
ctrl_ids <- readLines(file.path(PIPE, "matched", "any_ige", "ctrl_ids.txt"))
ids442 <- c(ige_ids, ctrl_ids)
ids442 <- ids442[ids442 %in% rownames(meta)]
g442 <- c(rep(TRUE, sum(ige_ids %in% rownames(meta))),
          rep(FALSE, sum(ctrl_ids %in% rownames(meta))))
t1_match <- build_table(ids442, g442, "Table1_matched442")
print(t1_match[, c("variable", "sensitized", "non_sensitized", "test",
                   "p_value", "q_BH")], row.names = FALSE, digits = 3)

message("\n  QC: matching variables should now be NON-significant by construction:")
for (v in c("Age (years)", "BMI (kg/m2)", "Sex")) {
  r <- t1_match[t1_match$variable == v, ]
  if (nrow(r))
    message(sprintf("    %-14s p=%.3f  %s", v, r$p_value,
                    if (r$p_value > 0.05) "PASS" else "WARNING - matching may be imperfect"))
}

# ── Table 2: the 9 IgE clusters ──────────────────────────────────────────
message("\n=== Table 2: IgE clusters A-I ===")
cl <- read.csv(file.path(KORA, "cluster_df_latent.csv"), row.names = 1)
rownames(cl) <- as.character(rownames(cl))
ids_cl <- intersect(rownames(cl), rownames(meta))
mc <- meta[ids_cl, , drop = FALSE]
clus <- factor(cl[ids_cl, "cluster"])
message(sprintf("  %d individuals across %d clusters (sizes: %s)",
                length(ids_cl), nlevels(clus),
                paste(sprintf("%s=%d", levels(clus), table(clus)), collapse = " ")))

t2_rows <- list()
for (nm in names(VARS)) {
  spec <- VARS[[nm]]
  if (!spec$col %in% names(mc)) next
  v <- mc[[spec$col]]
  if (spec$type == "cont") {
    v <- as.numeric(v)
    if (sum(!is.na(v)) < 10) next
    kt <- kruskal.test(v ~ clus)
    # eta-squared for Kruskal-Wallis
    eta <- (as.numeric(kt$statistic) - nlevels(clus) + 1) /
           (length(v) - nlevels(clus))
    t2_rows[[nm]] <- data.frame(
      variable = spec$label, column = spec$col, type = spec$type,
      test = "Kruskal-Wallis", statistic = as.numeric(kt$statistic),
      df = as.numeric(kt$parameter), p_value = kt$p.value,
      effect_name = "eta-squared", effect_size = eta, stringsAsFactors = FALSE)
  } else {
    tab <- table(factor(v), clus)
    if (any(dim(tab) < 2)) next
    ex <- suppressWarnings(chisq.test(tab, correct = FALSE)$expected)
    use_f <- any(ex < 5)
    tt <- if (use_f)
      tryCatch(fisher.test(tab, simulate.p.value = TRUE, B = 1e5),
               error = function(e) list(p.value = NA_real_))
      else suppressWarnings(chisq.test(tab, correct = FALSE))
    t2_rows[[nm]] <- data.frame(
      variable = spec$label, column = spec$col, type = spec$type,
      test = if (use_f) "Fisher's exact (simulated)" else "Chi-squared",
      statistic = if (use_f) NA_real_ else as.numeric(tt$statistic),
      df = NA_real_, p_value = tt$p.value,
      effect_name = "Cramer's V", effect_size = cramers_v(tab),
      stringsAsFactors = FALSE)
  }
}
t2 <- do.call(rbind, t2_rows)
t2$q_BH <- p.adjust(t2$p_value, method = "BH")
t2$underpowered_note <- "clusters F (n=10) and A (n=12) are underpowered"
rownames(t2) <- NULL
print(t2[, c("variable", "test", "statistic", "p_value", "q_BH",
             "effect_name", "effect_size")], row.names = FALSE, digits = 3)

# Per-cluster descriptives for the table body.
desc <- do.call(rbind, lapply(levels(clus), function(k) {
  s <- mc[clus == k, , drop = FALSE]
  data.frame(cluster = k, n = nrow(s),
             age_mean = mean(as.numeric(s$u3talter), na.rm = TRUE),
             age_sd = sd(as.numeric(s$u3talter), na.rm = TRUE),
             bmi_mean = mean(as.numeric(s$u3tbmi), na.rm = TRUE),
             bmi_sd = sd(as.numeric(s$u3tbmi), na.rm = TRUE),
             pct_female = 100 * mean(s$u3csex == "Female", na.rm = TRUE),
             stringsAsFactors = FALSE)
}))

# ── Table S1: IgE-burden partitions ──────────────────────────────────────
# Added AFTER the cluster block on purpose: the simulated Fisher tests draw from
# the seeded RNG, so appending (not interleaving) keeps every Table 2 p-value
# byte-identical.
#
# The six partitions OVERLAP (any = everyone; single + multiple = any;
# high + medium + low = any), so one test "across partitions" would be
# invalid. The two disjoint families are tested instead:
#   single vs multiple          two-group, same machinery as Table 1
#   high vs medium vs low       k-group, same machinery as Table 2
message("\n=== Table S1: IgE-burden partitions ===")
burden <- rowSums(ige != 0)
sens_ids <- names(burden)[burden >= 1]
parts <- list(
  any_ige      = sens_ids,
  single       = names(burden)[burden == 1],
  multiple     = names(burden)[burden > 1],
  high_score   = sens_ids[burden[sens_ids] > 15],
  medium_score = sens_ids[burden[sens_ids] <= 15 & burden[sens_ids] > 5],
  low_score    = sens_ids[burden[sens_ids] <= 5])
parts <- lapply(parts, function(x) intersect(as.character(x), rownames(meta)))

# (a) single vs multiple
ids_sm <- c(parts$single, parts$multiple)
g_sm <- c(rep(TRUE, length(parts$single)), rep(FALSE, length(parts$multiple)))
tS1_sm <- build_table(ids_sm, g_sm, "TableS1_single_vs_multiple")
tS1_sm$family <- "single_vs_multiple"

# (b) high vs medium vs low - identical code path to the Table 2 cluster tests
ids_b <- c(parts$high_score, parts$medium_score, parts$low_score)
bgrp <- factor(c(rep("high", length(parts$high_score)),
                 rep("medium", length(parts$medium_score)),
                 rep("low", length(parts$low_score))),
               levels = c("high", "medium", "low"))
mb <- meta[ids_b, , drop = FALSE]
b_rows <- list()
for (nm in names(VARS)) {
  spec <- VARS[[nm]]
  if (!spec$col %in% names(mb)) next
  v <- mb[[spec$col]]
  if (spec$type == "cont") {
    v <- as.numeric(v)
    kt <- kruskal.test(v ~ bgrp)
    eta <- (as.numeric(kt$statistic) - nlevels(bgrp) + 1) /
           (length(v) - nlevels(bgrp))
    b_rows[[nm]] <- data.frame(
      variable = spec$label, column = spec$col, type = spec$type,
      test = "Kruskal-Wallis", statistic = as.numeric(kt$statistic),
      df = as.numeric(kt$parameter), p_value = kt$p.value,
      effect_name = "eta-squared", effect_size = eta, stringsAsFactors = FALSE)
  } else {
    tab <- table(factor(v), bgrp)
    if (any(dim(tab) < 2)) next
    ex <- suppressWarnings(chisq.test(tab, correct = FALSE)$expected)
    use_f <- any(ex < 5)
    tt <- if (use_f)
      tryCatch(fisher.test(tab, simulate.p.value = TRUE, B = 1e5),
               error = function(e) list(p.value = NA_real_))
      else suppressWarnings(chisq.test(tab, correct = FALSE))
    b_rows[[nm]] <- data.frame(
      variable = spec$label, column = spec$col, type = spec$type,
      test = if (use_f) "Fisher's exact (simulated)" else "Chi-squared",
      statistic = if (use_f) NA_real_ else as.numeric(tt$statistic),
      df = if (use_f) NA_real_ else as.numeric(tt$parameter),
      p_value = tt$p.value,
      effect_name = "Cramer's V", effect_size = cramers_v(tab),
      stringsAsFactors = FALSE)
  }
}
tS1_b <- do.call(rbind, b_rows)
tS1_b$q_BH <- p.adjust(tS1_b$p_value, method = "BH")
tS1_b$family <- "high_vs_medium_vs_low"
rownames(tS1_b) <- NULL
print(tS1_sm[, c("variable", "test", "p_value", "q_BH")], row.names = FALSE, digits = 3)
print(tS1_b[, c("variable", "test", "p_value", "q_BH")], row.names = FALSE, digits = 3)

# Descriptives + the matched-pair count each microbiome analysis actually used.
n_pairs <- function(part) {
  src <- file.path(PIPE, "asv", "prev01", part, "source")
  a <- read.csv(file.path(src, "sample_data_ige.csv"))
  b <- read.csv(file.path(src, "sample_data_non_ige.csv"))
  length(intersect(a$pair_nb, b$pair_nb))
}
edu_share <- function(s) {
  e <- factor(s$u3tschul, levels = levels(meta$u3tschul))
  100 * as.numeric(table(e)) / sum(!is.na(e))
}
# Top components: the k most prevalent within the group (count of individuals
# with detectable IgE), with EVERY component tied with the k-th listed - so no
# arbitrary tie-break is ever needed. Bromelain excluded (a cross-reactive
# carbohydrate-determinant marker, not an allergen source). Written as data,
# not chosen by eye: two submitted Table S1 rows did not match any ranking
# (single listed Der p 10, detected in 0 of 47; low listed Phl p 5, 8%, over
# Phl p 1, 18%).
top_components <- function(ids, k = 3) {
  x <- ige[ids, setdiff(colnames(ige), "bromelain"), drop = FALSE]
  cnt <- colSums(x != 0)
  o <- order(-cnt)
  cut <- cnt[o][k]
  keep <- o[cnt[o] >= cut & cnt[o] > 0]
  paste(sprintf("%s:%d/%d", colnames(x)[keep], cnt[keep], length(ids)),
        collapse = ";")
}
all_groups <- c(parts, setNames(lapply(levels(clus), function(k)
  intersect(as.character(ids_cl[clus == k]), rownames(meta))),
  paste0("cluster_", levels(clus))))
descS1 <- do.call(rbind, lapply(names(all_groups), function(g) {
  ids <- all_groups[[g]]; s <- meta[ids, , drop = FALSE]
  ed <- edu_share(s)
  data.frame(group = g, n = length(ids), n_pairs = n_pairs(g),
             age_mean = mean(as.numeric(s$u3talter), na.rm = TRUE),
             age_sd = sd(as.numeric(s$u3talter), na.rm = TRUE),
             bmi_mean = mean(as.numeric(s$u3tbmi), na.rm = TRUE),
             bmi_sd = sd(as.numeric(s$u3tbmi), na.rm = TRUE),
             pct_male = 100 * mean(s$u3csex == "Male", na.rm = TRUE),
             edu_s = ed[1], edu_m = ed[2], edu_h = ed[3],
             top3_components = top_components(ids),
             stringsAsFactors = FALSE)
}))
print(descS1[, c("group", "n", "n_pairs", "top3_components")], row.names = FALSE)

# ── Write outputs ────────────────────────────────────────────────────────
write.csv(t1_full,  file.path(OUT, "table1_full508_tested.csv"), row.names = FALSE)
write.csv(t1_match, file.path(OUT, "table1_matched442_qc.csv"), row.names = FALSE)
write.csv(t2,       file.path(OUT, "table2_clusters_tested.csv"), row.names = FALSE)
write.csv(desc,     file.path(OUT, "table2_clusters_descriptives.csv"), row.names = FALSE)
tS1_sm$df <- NA_real_
cols_S1 <- c("family", "variable", "column", "type", "test", "statistic", "df",
             "p_value", "effect_name", "effect_size", "q_BH")
write.csv(rbind(tS1_sm[, cols_S1], tS1_b[, cols_S1]),
          file.path(OUT, "tableS1_partitions_tested.csv"), row.names = FALSE)
write.csv(descS1, file.path(OUT, "group_descriptives_pairs_topIgE.csv"), row.names = FALSE)

message(sprintf("\nWrote 4 tables to %s", OUT))
message(sprintf("Validation failures: %d", n_fail))
if (n_fail > 0) quit(status = 1)
message("Done.")
