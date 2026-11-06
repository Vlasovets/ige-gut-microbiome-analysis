#!/usr/bin/env Rscript
# Module 05 - beta diversity, PERMANOVA and dispersion.
#
# Answers R2 Major 2 and R3 Major 2/5. The submitted manuscript contains NO
# beta-diversity analysis at all, despite framing itself against a prior
# beta-diversity finding (Goldberg et al.). This module supplies it.
#
# Four decisions that differ from the existing pipeline, each for a stated reason:
#
# 1. by = "margin" IS PRIMARY. The existing script reports only by = "terms" with
#    ige_status entered FIRST, which is the most favourable ordering for the term
#    of interest. Type III (margin) tests each term against the full model. We
#    also report by = "terms" with ige_status entered first AND last, so the
#    ordering sensitivity is visible rather than hidden.
#
# 2. betadisper IS MANDATORY, not optional. A significant adonis2 with unequal
#    dispersion is a dispersion result, not a location result. Reported alongside
#    every PERMANOVA row.
#
# 3. FOUR distances, including robust Aitchison. The pipeline used Bray-Curtis and
#    a pseudo-count CLR ("Aitchison"). rclr avoids the pseudo-count entirely.
#    Jaccard adds a presence/absence view. UniFrac is NOT computed: no phylogeny
#    or representative sequences exist for this dataset (see REVISION_PLAN B2).
#
# 4. COLLIDERS ARE EXCLUDED FROM ADJUSTMENT. Asthma, dermatitis and hay fever are
#    downstream of sensitization; conditioning on them induces collider bias.
#    They are tabulated as outcomes in module 01, never adjusted for here.
#
# Permutations are restricted WITHIN matched pairs throughout, which is the
# exchangeability the 1:1 design implies.
#
# Usage:
#   Rscript revision/analysis/05_beta_diversity.R \
#     --level family --prevalence prev10 --partition any_ige

suppressPackageStartupMessages({
  library(yaml); library(vegan); library(permute)
})

args_raw <- commandArgs(trailingOnly = FALSE)
script <- sub("--file=", "", args_raw[grep("--file=", args_raw)])
REVISION <- normalizePath(file.path(dirname(script), ".."))
.libPaths(c(file.path(REVISION, "env", "rlibs"), .libPaths()))

a <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- which(a == f); if (length(i) && i < length(a)) a[i+1] else d }
LEVEL <- getarg("--level", "family")
PREV  <- getarg("--prevalence", "prev10")
PART  <- getarg("--partition", "any_ige")
NPERM <- as.integer(getarg("--nperm", "5000"))

cfg  <- yaml::read_yaml(file.path(REVISION, "config", "config.yaml"))
PIPE <- cfg$paths$pipeline_results
OUT  <- file.path(REVISION, "results", "beta")
FIGS <- file.path(REVISION, "results", "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)
SEED <- cfg$statistics$seed

message(sprintf("=== Beta diversity: %s / %s / %s (nperm=%d) ===",
                LEVEL, PREV, PART, NPERM))

# ── Load matched data ────────────────────────────────────────────────────
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
meta$sex        <- factor(meta$u3csex)
meta$age        <- as.numeric(meta$u3talter)
meta$bmi        <- as.numeric(meta$u3tbmi)
meta$education  <- factor(meta$u3tschul)
meta$smoking    <- factor(meta$u3tcigreg_sf)
meta$packyears  <- as.numeric(meta$u3tpyrs)
meta$alcohol    <- as.numeric(meta$u3talkkon)
meta$activity   <- factor(ifelse(meta$u3tphact %in% c(1, 2), "regular",
                          ifelse(meta$u3tphact == 3, "irregular", "none")))
meta$ppi        <- factor(meta$u3tmppi)
meta$pair_nb    <- factor(meta$pair_nb)

# Drop samples missing anything the widest model needs, so all three models are
# fitted on an IDENTICAL sample set and their R2 values are comparable.
need <- c("age","bmi","sex","education","smoking","packyears","alcohol","activity","ppi")
ok <- stats::complete.cases(meta[, need]) & rowSums(otu) > 0
message(sprintf("  samples: %d of %d retained (complete covariates, non-zero)",
                sum(ok), nrow(meta)))
otu <- otu[ok, , drop = FALSE]; meta <- meta[ok, , drop = FALSE]
otu <- otu[rownames(meta), , drop = FALSE]

# Matched pairs must remain intact: drop any pair that lost a member, otherwise
# within-pair permutation is ill-defined.
tab <- table(meta$pair_nb)
full_pairs <- names(tab)[tab == 2]
drop_n <- sum(!meta$pair_nb %in% full_pairs)
if (drop_n) {
  message(sprintf("  dropping %d sample(s) from broken pairs", drop_n))
  keep2 <- meta$pair_nb %in% full_pairs
  otu <- otu[keep2, , drop = FALSE]; meta <- meta[keep2, , drop = FALSE]
}
meta$pair_nb <- droplevels(meta$pair_nb)
message(sprintf("  final: n=%d (%d complete pairs), taxa=%d",
                nrow(meta), nlevels(meta$pair_nb), ncol(otu)))
if (nlevels(meta$pair_nb) < 5) {
  message("  too few pairs for a permutation test; writing a stub and exiting")
  write.csv(data.frame(level=LEVEL, prevalence=PREV, partition=PART,
                       note="skipped: <5 complete pairs"),
            file.path(OUT, sprintf("beta_%s_%s_%s.csv", LEVEL, PREV, PART)),
            row.names = FALSE)
  quit(status = 0)
}

# ── Distances ────────────────────────────────────────────────────────────
clr_pseudo <- function(x) { l <- log(x + 1); l - rowMeans(l) }
dists <- list(
  bray             = vegdist(otu, method = "bray"),
  jaccard          = vegdist(otu, method = "jaccard", binary = TRUE),
  aitchison_pseudo = dist(clr_pseudo(otu)),
  aitchison_robust = dist(decostand(otu, method = "rclr"))
)
for (nm in names(dists)) {
  d <- dists[[nm]]
  if (any(!is.finite(as.vector(d)))) {
    message(sprintf("  WARNING: %s contains non-finite values; dropped", nm))
    dists[[nm]] <- NULL
  }
}

# ── Models ───────────────────────────────────────────────────────────────
MODELS <- list(
  M1 = "ige_status + sex + age + bmi",
  M2 = "ige_status + sex + age + bmi + education + smoking + packyears + alcohol + activity",
  M3 = "ige_status + sex + age + bmi + education + smoking + packyears + alcohol + activity + ppi",
  # ige_status entered LAST, to expose the term-order sensitivity that reporting
  # only by="terms" with it first would conceal.
  M1_last = "sex + age + bmi + ige_status"
)

perm_ctrl <- how(nperm = NPERM, blocks = meta$pair_nb)
rows <- list()

for (dn in names(dists)) {
  d <- dists[[dn]]
  for (mn in names(MODELS)) {
    for (by in c("margin", "terms")) {
      # by="margin" on M1_last is identical to M1; skip the duplicate.
      if (mn == "M1_last" && by == "margin") next
      f <- as.formula(paste("d ~", MODELS[[mn]]))
      set.seed(SEED)
      ad <- tryCatch(
        adonis2(f, data = meta, by = by, permutations = perm_ctrl),
        error = function(e) { message(sprintf("    adonis2 %s/%s/%s: %s",
                                              dn, mn, by, conditionMessage(e))); NULL })
      if (is.null(ad)) next
      terms_i <- setdiff(rownames(ad), c("Residual", "Total"))
      for (tm in terms_i) {
        rows[[length(rows) + 1]] <- data.frame(
          level = LEVEL, prevalence = PREV, partition = PART,
          distance = dn, model = mn, by = by, term = tm,
          Df = ad[tm, "Df"], R2 = ad[tm, "R2"],
          F_stat = ad[tm, "F"], p_value = ad[tm, "Pr(>F)"],
          n = nrow(meta), n_pairs = nlevels(meta$pair_nb), n_taxa = ncol(otu),
          nperm = NPERM, stringsAsFactors = FALSE)
      }
    }
  }
}
res <- do.call(rbind, rows)

# BH across the terms tested within each distance x model x by.
res$q_BH <- NA_real_
for (k in unique(paste(res$distance, res$model, res$by))) {
  i <- which(paste(res$distance, res$model, res$by) == k)
  res$q_BH[i] <- p.adjust(res$p_value[i], method = "BH")
}

message("\n=== PERMANOVA: ige_status term ===")
ig <- res[res$term == "ige_status", ]
print(ig[order(ig$distance, ig$model, ig$by),
         c("distance","model","by","R2","F_stat","p_value","q_BH")],
      row.names = FALSE, digits = 4)

# The comparison that carries the paper: is the IgE effect larger than the sex effect?
message("\n=== Effect size in context (M1, by=margin) ===")
ctx <- res[res$model == "M1" & res$by == "margin", ]
for (dn in unique(ctx$distance)) {
  s <- ctx[ctx$distance == dn, ]
  s <- s[order(-s$R2), ]
  message(sprintf("  %-18s %s", dn,
                  paste(sprintf("%s R2=%.4f(p=%.3f)", s$term, s$R2, s$p_value),
                        collapse = "  ")))
}

# ── Dispersion — mandatory companion ─────────────────────────────────────
message("\n=== betadisper (a location result requires equal dispersion) ===")
disp_rows <- list()
for (dn in names(dists)) {
  bd <- tryCatch(betadisper(dists[[dn]], meta$ige_status, type = "median",
                            bias.adjust = TRUE),
                 error = function(e) NULL)
  if (is.null(bd)) next
  set.seed(SEED)
  pt <- permutest(bd, permutations = how(nperm = NPERM, blocks = meta$pair_nb))
  p_disp <- pt$tab[1, "Pr(>F)"]
  mdist <- tapply(bd$distances, meta$ige_status, mean)
  disp_rows[[dn]] <- data.frame(
    level = LEVEL, prevalence = PREV, partition = PART, distance = dn,
    F_stat = pt$tab[1, "F"], p_value = p_disp,
    mean_dist_non_ige = mdist[["non_ige"]], mean_dist_ige = mdist[["ige"]],
    interpretation = if (!is.na(p_disp) && p_disp < 0.05)
      "UNEQUAL dispersion - a significant PERMANOVA may reflect spread, not location"
      else "equal dispersion - PERMANOVA interpretable as a location shift",
    stringsAsFactors = FALSE)
  message(sprintf("  %-18s F=%.3f p=%.4f  mean dist: non_ige=%.3f ige=%.3f  %s",
                  dn, pt$tab[1, "F"], p_disp, mdist[["non_ige"]], mdist[["ige"]],
                  if (!is.na(p_disp) && p_disp < 0.05) "UNEQUAL" else "equal"))
}
disp <- do.call(rbind, disp_rows)

# ── PCoA ─────────────────────────────────────────────────────────────────
if (PART == "any_ige") {
  cols <- c(non_ige = "#2980B9", ige = "#C0392B")
  pdf(file.path(FIGS, sprintf("pcoa_%s_%s_%s.pdf", LEVEL, PREV, PART)),
      width = 9, height = 4.6, pointsize = 9)
  par(mfrow = c(1, min(2, length(dists))), mar = c(4, 4, 3, 1))
  for (dn in head(names(dists), 2)) {
    pc <- cmdscale(dists[[dn]], k = 2, eig = TRUE)
    ve <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
    r <- ig[ig$distance == dn & ig$model == "M1" & ig$by == "margin", ]
    pd <- disp[disp$distance == dn, ]
    plot(pc$points, col = cols[as.character(meta$ige_status)], pch = 19, cex = 0.55,
         xlab = sprintf("PCo1 (%.1f%%)", ve[1]),
         ylab = sprintf("PCo2 (%.1f%%)", ve[2]),
         main = sprintf("%s\nR2=%.4f p=%.3f | disp p=%.3f",
                        dn, if (nrow(r)) r$R2 else NA,
                        if (nrow(r)) r$p_value else NA,
                        if (nrow(pd)) pd$p_value else NA),
         cex.main = 0.95)
    ordiellipse(pc$points, meta$ige_status, kind = "sd", conf = 0.95,
                col = cols, lwd = 1.6)
    legend("topright", legend = names(cols), col = cols, pch = 19,
           bty = "n", cex = 0.9)
  }
  dev.off()
  message(sprintf("\n  Wrote %s",
                  file.path(FIGS, sprintf("pcoa_%s_%s_%s.pdf", LEVEL, PREV, PART))))
}

# ── Write ────────────────────────────────────────────────────────────────
write.csv(res,  file.path(OUT, sprintf("beta_%s_%s_%s.csv", LEVEL, PREV, PART)),
          row.names = FALSE)
if (!is.null(disp))
  write.csv(disp, file.path(OUT, sprintf("betadisper_%s_%s_%s.csv", LEVEL, PREV, PART)),
            row.names = FALSE)
message(sprintf("\nWrote %d PERMANOVA rows to %s", nrow(res), OUT))
message("Done.")
