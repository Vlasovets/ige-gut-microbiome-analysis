#!/usr/bin/env Rscript
# Module 07b - species-level input of the primary TaxSEA analysis, with set
# membership of the two headline sets, for plotting (module 06t).
#
# Primary analysis = module 07 --design full: species ranked by the signed LinDA
# statistic (all 275 sensitized vs 233 non-sensitized, ~ IgE + sex + age + BMI),
# 5% prevalence filter. Name handling is identical to module 07 (SILVA brackets
# stripped, duplicates collapsed to the stronger signal), and membership is
# resolved the way TaxSEA resolves it: species name -> NCBI id (NCBI_ids) ->
# id contained in the set (TaxSEA_db).
#
# Output (species-level summary statistics only, no participant data):
#   results/tables/taxsea_species_membership_full.csv
#
# Environment: taxsea-env (R >= 4.5, TaxSEA 1.2.0)
# Usage: Rscript revision/analysis/07b_taxsea_members.R

suppressPackageStartupMessages(library(TaxSEA))

args_raw <- commandArgs(trailingOnly = FALSE)
script <- sub("--file=", "", args_raw[grep("--file=", args_raw)])
REVISION <- normalizePath(file.path(dirname(script), ".."))

data("NCBI_ids", package = "TaxSEA")
data("TaxSEA_db", package = "TaxSEA")
HEADLINE <- c(folic_acid = "MiMeDB_producers_of_Folic_acid",
              vitamin_a  = "MiMeDB_producers_of_Vitamin_A")

d <- read.csv(file.path(REVISION, "results", "da_full",
                        "da_long_species_prev05_any_ige.csv"), stringsAsFactors = FALSE)
d <- d[d$method == "LinDA" & !is.na(d$stat) & nzchar(d$taxon), ]
d$taxon_taxsea <- gsub("\\[|\\]", "", d$taxon)
d <- d[order(-abs(d$stat)), ]
d <- d[!duplicated(d$taxon_taxsea), ]

ids <- NCBI_ids[d$taxon_taxsea]
d$ncbi_id <- vapply(ids, function(x) if (is.null(x)) NA_character_ else as.character(x), "")
d$mapped <- !is.na(d$ncbi_id)

# TaxSEA's background (TaxSEA 1.2.0 source): mapped taxa that belong to at least
# one set whose size, counted among the input taxa, is within [5, 100].
in_db <- d$ncbi_id[d$mapped & d$ncbi_id %in% unique(unlist(TaxSEA_db))]
sizes <- vapply(TaxSEA_db, function(s) length(unique(intersect(s, in_db))), numeric(1))
kept <- TaxSEA_db[sizes >= 5 & sizes <= 100]
d$in_background <- d$mapped & d$ncbi_id %in% unique(unlist(kept))
for (k in names(HEADLINE))
  d[[paste0("in_", k)]] <- d$in_background & d$ncbi_id %in% unlist(TaxSEA_db[[HEADLINE[[k]]]])

out <- d[, c("taxon", "taxon_taxsea", "ncbi_id", "mapped", "in_background", "stat",
             "effect_log2", "p_raw", "q_within", paste0("in_", names(HEADLINE)))]
out <- out[order(out$stat), ]
f <- file.path(REVISION, "results", "tables", "taxsea_species_membership_full.csv")
write.csv(out, f, row.names = FALSE)

message(sprintf("species %d, mapped %d, TaxSEA background %d, folic acid members %d, vitamin A members %d",
                nrow(out), sum(out$mapped), sum(out$in_background),
                sum(out$in_folic_acid), sum(out$in_vitamin_a)))
# Cross-check against module 07 (taxsea_recomputed_full.csv, prev05 / stat)
ref <- read.csv(file.path(REVISION, "results", "tables", "taxsea_recomputed_full.csv"))
ref <- ref[ref$prevalence == "prev05" & ref$rank_by == "stat" &
           ref$sub_database == "Metabolite_producers", ]
for (k in names(HEADLINE)) {
  r <- out$stat[out$in_background]; s <- out$stat[out[[paste0("in_", k)]]]
  p <- suppressWarnings(ks.test(s, r)$p.value)
  p_ref <- ref$PValue[ref$taxonSetName == HEADLINE[[k]]]
  message(sprintf("  %-10s n=%d  median=%.3f  KS p=%.4f  (module 07: %.4f)", k, length(s), median(s), p, p_ref))
  stopifnot(abs(p - p_ref) < 1e-6)
}
message("wrote ", f)
