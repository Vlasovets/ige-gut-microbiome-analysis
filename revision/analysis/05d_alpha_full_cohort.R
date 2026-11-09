#!/usr/bin/env Rscript
# Module 05d - alpha diversity in the FULL-COHORT, covariate-adjusted design.
#
# Same design as the revision's DA: every IgE group vs all 233 IgE-negative,
# adjusted for sex, age and BMI (no matching). Indices as in the submitted
# Figure 4A: Shannon (raw counts, from module 05c) and Faith's PD (computed here
# for all 508; module 05b computed it only for matched samples, 489 of 508).
#
# Faith's PD: same tree, same ASV set (non-zero in the 508 subset, 6,958 ASVs)
# and same implementation as 05b, so per-sample values are identical where they
# overlap (asserted).
#
# Test per group x index: lm(index ~ group + sex + age + bmi); p for the group
# term; BH across the 15 groups within each index.
#
# Output: results/alpha/alpha_full_persample.csv   (per participant; gitignored)
#         results/alpha/alpha_full_tests.csv
# Usage:  Rscript revision/analysis/05d_alpha_full_cohort.R

suppressPackageStartupMessages({ library(yaml); library(ape) })
args_raw <- commandArgs(trailingOnly = FALSE)
script <- sub("--file=", "", args_raw[grep("--file=", args_raw)])
REVISION <- normalizePath(file.path(dirname(script), ".."))
cfg  <- yaml::read_yaml(file.path(REVISION, "config", "config.yaml"))
KORA <- cfg$paths$kora_data
RES  <- file.path(REVISION, "results")

faith_pd <- function(tree, present) {          # identical to 05b
  present <- intersect(present, tree$tip.label)
  if (length(present) < 2) return(NA_real_)
  sum(ape::keep.tip(tree, present)$edge.length, na.rm = TRUE)
}
tree <- ape::read.tree(file.path(KORA, "q2_reprocessing", "rooted_tree.nwk"))

base <- read.csv(file.path(RES, "nested_models", "designs", "base.csv"), colClasses = "character")
ctrl <- base$sample_id[base$partition == "__controls__"]
ids  <- unique(base$sample_id)
stopifnot(length(ctrl) == 233, length(ids) == 508)

ft <- read.csv(file.path(KORA, "feature_table.tsv"), sep = "\t", row.names = 1, check.names = FALSE)
ft <- ft[, ids]
ft <- ft[rowSums(ft > 0) > 0, ]
stopifnot(nrow(ft) == 6958)
pd <- vapply(ids, function(s) faith_pd(tree, rownames(ft)[ft[[s]] > 0]), numeric(1))

# consistency with 05b where samples overlap
old <- do.call(rbind, lapply(Sys.glob(file.path(RES, "phylo", "faithpd_persample_*.csv")),
                             function(f) read.csv(f, colClasses = c(sample = "character"))))
old <- old[!duplicated(old$sample), ]
dif <- max(abs(pd[old$sample] - old$pd))
message(sprintf("Faith PD: %d samples; max |diff| vs 05b on %d overlapping = %.2e",
                length(pd), nrow(old), dif))
stopifnot(dif < 1e-8)

al <- read.csv(file.path(RES, "alpha", "alpha_persample.csv"), colClasses = c(sample = "character"))
al <- al[al$variant == "raw", ]
kora <- read.csv(file.path(KORA, "kora.csv"), check.names = FALSE)
kora$id <- as.character(as.integer(kora$u3_16s_id))
kora <- kora[!duplicated(kora$id), ]
ps <- data.frame(sample = ids, shannon = al$shannon[match(ids, al$sample)], pd = pd[ids],
                 sex = factor(kora$u3csex[match(ids, kora$id)]),
                 age = as.numeric(kora$u3talter[match(ids, kora$id)]),
                 bmi = as.numeric(kora$u3tbmi[match(ids, kora$id)]))
stopifnot(!anyNA(ps))
write.csv(ps, file.path(RES, "alpha", "alpha_full_persample.csv"), row.names = FALSE)

parts <- c("any_ige", "single", "multiple", "high_score", "medium_score", "low_score",
           paste0("cluster_", LETTERS[1:9]))
rows <- list()
for (p in parts) {
  cases <- base$sample_id[base$partition == p & base$role == "case"]
  d <- ps[ps$sample %in% c(cases, ctrl), ]
  d$group <- factor(ifelse(d$sample %in% cases, "ige", "non_ige"), levels = c("non_ige", "ige"))
  for (ix in c("shannon", "pd")) {
    co <- summary(lm(as.formula(paste(ix, "~ group + sex + age + bmi")), data = d))$coefficients
    rows[[length(rows) + 1]] <- data.frame(partition = p, index = ix, n_ige = length(cases),
      n_ctrl = length(ctrl), estimate = co["groupige", "Estimate"],
      se = co["groupige", "Std. Error"], p_value = co["groupige", "Pr(>|t|)"])
  }
}
t <- do.call(rbind, rows)
t$q_value <- ave(t$p_value, t$index, FUN = function(p) p.adjust(p, "BH"))
write.csv(t, file.path(RES, "alpha", "alpha_full_tests.csv"), row.names = FALSE)
print(aggregate(cbind(min_p = p_value, min_q = q_value) ~ index, t, min))
