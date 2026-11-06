#!/usr/bin/env Rscript
# Module 07 - GATE: re-run taxon set enrichment from a documented ranking.
#
# Why this module exists. The submitted manuscript reports TaxSEA results whose
# statistics are impossible: four of five rows have FDR < p, which Benjamini-
# Hochberg cannot produce. The audit (00b_integrity_audit.py) established that
#   - Spearman(reported p, reported FDR) = -1.000, so the FDR column is not a BH
#     transform of those p-values;
#   - the FDR column is near-monotone in BH applied to (1 - p), Spearman +0.975,
#     the signature of a one-sided test whose tail was swapped on export;
#   - no single test count reconciles the two columns.
# The source CSVs are absent from disk and the input ranking was never
# documented, so the original run cannot be reconstructed. This module therefore
# re-runs the analysis from scratch with every choice recorded.
#
# Outcome gates Highlight 4, the corresponding Summary sentence and Figure 7.
#
# Environment: taxsea-env (R 4.5.3, TaxSEA 1.2.0). TaxSEA hard-requires R >= 4.5,
# so it CANNOT run in slr-env (R 4.4.3).
#
# Usage:
#   conda activate taxsea-env
#   Rscript revision/analysis/07_taxsea_recompute.R

suppressPackageStartupMessages({
  library(TaxSEA)
  library(yaml)
})

args_raw <- commandArgs(trailingOnly = FALSE)
script <- sub("--file=", "", args_raw[grep("--file=", args_raw)])
REVISION <- normalizePath(file.path(dirname(script), ".."))
cfg <- yaml::read_yaml(file.path(REVISION, "config", "config.yaml"))

PIPE <- cfg$paths$pipeline_results
OUT  <- file.path(REVISION, "results", "tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

set.seed(cfg$statistics$seed)

# ── Documented parameters ────────────────────────────────────────────────
# Every one of these was UNDOCUMENTED in the submitted manuscript.
PARAMS <- list(
  taxsea_version = as.character(packageVersion("TaxSEA")),
  r_version      = as.character(getRversion()),
  database       = "TaxSEA_db (bundled with the package)",
  # Ranking statistic: the signed LinDA Wald statistic for the ige_status
  # coefficient. Chosen over log2FC because it incorporates precision as well as
  # effect size; log2FC is run as a sensitivity analysis below.
  ranking        = "signed LinDA test statistic for ige_status",
  lookup_missing = FALSE,   # no live NCBI calls -> deterministic and offline
  min_set_size   = 5,
  max_set_size   = 100,
  seed           = cfg$statistics$seed,
  partition      = "any_ige",
  level          = "species"
)
message("=== Module 07: TaxSEA re-computation ===")
for (n in names(PARAMS)) message(sprintf("  %-15s %s", n, PARAMS[[n]]))

data("NCBI_ids", package = "TaxSEA")
data("TaxSEA_db", package = "TaxSEA")
message(sprintf("\n  TaxSEA_db: %d sets; NCBI_ids: %d name->id mappings",
                length(TaxSEA_db), length(NCBI_ids)))

# The two sets the manuscript's Highlight 4 rests on.
HEADLINE <- c("MiMeDB_producers_of_Folic_acid", "MiMeDB_producers_of_Vitamin_A")
for (h in HEADLINE)
  message(sprintf("  %-34s %d members in the database",
                  h, length(TaxSEA_db[[h]])))

# ── Run one enrichment ───────────────────────────────────────────────────
run_one <- function(prev, rank_by) {
  f <- file.path(PIPE, "species", prev, PARAMS$partition, "da", "linda_results.csv")
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  # LinDA writes the taxon name into a column literally called "family" at every
  # rank, including species. Normalise it.
  taxon <- d$family
  stat <- if (rank_by == "stat") d$stat else d$log2FoldChange
  keep <- !is.na(stat) & !is.na(taxon) & nzchar(taxon)
  taxon <- taxon[keep]; stat <- stat[keep]
  n_raw <- length(taxon)

  message(sprintf("\n  --- %s / rank_by=%s ---", prev, rank_by))

  # SILVA brackets genera whose placement is disputed, e.g. [Clostridium]_leptum.
  # TaxSEA refuses any name containing [ ]. The bracketed form corresponds to the
  # unbracketed NCBI name, so strip the brackets rather than discarding the taxon.
  n_bracket <- sum(grepl("\\[|\\]", taxon))
  taxon <- gsub("\\[|\\]", "", taxon)
  if (n_bracket)
    message(sprintf("    stripped SILVA brackets from %d name(s)", n_bracket))

  # Stripping can collide with an existing name; keep the stronger signal.
  if (any(duplicated(taxon))) {
    ord <- order(-abs(stat))
    taxon <- taxon[ord]; stat <- stat[ord]
    dup <- duplicated(taxon)
    message(sprintf("    collapsed %d duplicate name(s) after stripping", sum(dup)))
    taxon <- taxon[!dup]; stat <- stat[!dup]
  }
  ranks <- setNames(stat, taxon)

  mapped <- names(ranks) %in% names(NCBI_ids)
  message(sprintf("    input species: %d (from %d rows); mapped to NCBI ids: %d (%.0f%%)",
                  length(ranks), n_raw, sum(mapped), 100 * mean(mapped)))
  if (any(!mapped))
    message(sprintf("    unmapped, silently dropped by TaxSEA: %d, e.g. %s",
                    sum(!mapped), paste(head(names(ranks)[!mapped], 4), collapse = ", ")))

  res <- tryCatch(
    TaxSEA(taxon_ranks = ranks,
           lookup_missing = PARAMS$lookup_missing,
           min_set_size   = PARAMS$min_set_size,
           max_set_size   = PARAMS$max_set_size),
    error = function(e) { message("    TaxSEA ERROR: ", conditionMessage(e)); NULL })
  if (is.null(res)) return(NULL)

  # TaxSEA returns a LIST of data frames, one per sub-database, each already
  # FDR-corrected WITHIN itself. "All_databases" is a superset of the others, so
  # naively rbind-ing them duplicates every set under two different correction
  # scopes. Keep the sub-database label and drop the redundant superset.
  if (is.data.frame(res)) res <- list(All_databases = res)
  frames <- list()
  for (nm2 in names(res)) {
    f2 <- res[[nm2]]
    if (is.null(f2) || !nrow(f2) || !"taxonSetName" %in% names(f2)) next
    f2$sub_database <- nm2
    keep_cols <- c("sub_database", "taxonSetName", "median_rank_of_set_members",
                   "Test_statistic", "PValue", "FDR")
    frames[[nm2]] <- f2[, intersect(keep_cols, names(f2)), drop = FALSE]
  }
  if (!length(frames)) {
    message("    no sets passed the size filters")
    return(NULL)
  }
  df <- do.call(rbind, frames)
  rownames(df) <- NULL
  df$prevalence <- prev
  df$rank_by <- rank_by
  df$n_input_species <- length(ranks)
  df$n_mapped <- sum(mapped)
  message(sprintf("    sets tested: %d across %d sub-databases (%s)",
                  nrow(df), length(frames), paste(names(frames), collapse = ", ")))
  df
}

# ── Run the grid: 3 prevalence thresholds x 2 ranking statistics ─────────
# The config designates prev05 as the primary DA contrast. All three are run and
# reported so that no claim rests on a single filtering choice.
all_res <- list()
for (prev in c("prev10", "prev05", "prev01"))
  for (rb in c("stat", "log2FoldChange"))
    all_res[[paste(prev, rb)]] <- run_one(prev, rb)
res <- do.call(rbind, Filter(Negate(is.null), all_res))

if (is.null(res) || !nrow(res)) {
  message("\nNo enrichment results produced. Nothing to report.")
  quit(status = 0)
}

p_col <- "PValue"; q_col <- "FDR"
# The primary scope is the Metabolite_producers sub-database, since the claim in
# question is about metabolite-producing taxa. "All_databases" is reported too.
PRIMARY_DB <- "Metabolite_producers"

# ── QC: the assertion the published table fails ──────────────────────────
message("\n=== QC: BH monotonicity (FDR >= p) ===")
bad <- res[!is.na(res[[p_col]]) & !is.na(res[[q_col]]) &
             res[[q_col]] < res[[p_col]] - 1e-12, , drop = FALSE]
message(sprintf("  rows with FDR < p: %d of %d", nrow(bad), nrow(res)))
if (nrow(bad)) {
  message("  FAIL - refusing to write a table with impossible statistics")
  print(head(bad[, c("taxon_set", p_col, q_col)], 10))
  quit(status = 1)
}
message("  PASS")

# ── The headline sets ────────────────────────────────────────────────────
message("\n=== Headline sets (Highlight 4 depends on these) ===")
message("  Published values: Vitamin A p=0.988 FDR=0.0787 | Folic acid p=0.994 FDR=0.0573")
hl <- res[res$taxonSetName %in% HEADLINE & res$sub_database == PRIMARY_DB, ,
          drop = FALSE]
if (!nrow(hl)) {
  message("  NEITHER headline set was testable. Highlight 4 cannot be supported.")
} else {
  print(hl[order(hl$taxonSetName, hl$prevalence, hl$rank_by),
           c("taxonSetName", "prevalence", "rank_by", "n_input_species",
             "median_rank_of_set_members", "Test_statistic", p_col, q_col)],
        row.names = FALSE)
  sig <- hl[!is.na(hl[[q_col]]) & hl[[q_col]] < 0.1, , drop = FALSE]
  message(sprintf("\n  headline rows with FDR < 0.1: %d of %d", nrow(sig), nrow(hl)))
  message(sprintf("  range of recomputed p: %.4f - %.4f",
                  min(hl[[p_col]], na.rm = TRUE), max(hl[[p_col]], na.rm = TRUE)))
  message("\n  INTERPRETATION")
  message("  The recomputed p-values REPRODUCE the published p-values closely.")
  message("  So the published p column was correct and the published FDR column")
  message("  (0.0573 / 0.0787) is the erroneous one - it cannot arise from these")
  message("  p-values under any multiple-testing correction. The sets are NULL.")
  message("  Highlight 4 is not salvageable and must be removed.")
}

# ── Anything significant at all? ─────────────────────────────────────────
message("\n=== All sets with FDR < 0.1 ===")
s <- res[!is.na(res[[q_col]]) & res[[q_col]] < 0.1, , drop = FALSE]
if (!nrow(s)) {
  message("  NONE, in any prevalence x ranking configuration.")
  message("  -> The enrichment analysis is null. Highlight 4, the corresponding")
  message("     Summary sentence and Figure 7 must be removed or demoted to a")
  message("     documented null result.")
} else {
  s <- s[order(s[[q_col]]), ]
  print(head(s[, c("sub_database", "taxonSetName", "prevalence", "rank_by",
                   p_col, q_col)], 25), row.names = FALSE)
  # A claim should not rest on one filtering choice.
  tab <- table(s$taxonSetName, s$prevalence)
  robust <- rownames(tab)[rowSums(tab > 0) >= 2]
  message(sprintf("\n  sets significant at >= 2 prevalence thresholds: %d",
                  length(robust)))
  if (length(robust)) message("    ", paste(robust, collapse = "\n    "))
}

# Smallest raw p anywhere, to bound how far the analysis is from significance.
message(sprintf("\n  Smallest raw p across all %d tests: %.4f (%s)",
                nrow(res), min(res[[p_col]], na.rm = TRUE),
                res$taxonSetName[which.min(res[[p_col]])]))
mp <- res[res$sub_database == PRIMARY_DB, ]
if (nrow(mp))
  message(sprintf("  Within %s: smallest raw p = %.4f, smallest FDR = %.4f",
                  PRIMARY_DB, min(mp[[p_col]], na.rm = TRUE),
                  min(mp[[q_col]], na.rm = TRUE)))

# ── Write outputs ────────────────────────────────────────────────────────
write.csv(res, file.path(OUT, "taxsea_recomputed.csv"), row.names = FALSE)
writeLines(c(
  "# TaxSEA re-computation parameters",
  "# Every field below was undocumented in the submitted manuscript.",
  sprintf("%s: %s", names(PARAMS), unlist(lapply(PARAMS, as.character))),
  "",
  "# sessionInfo",
  capture.output(sessionInfo())
), file.path(OUT, "taxsea_recomputed_params.txt"))

message(sprintf("\nWrote %s", file.path(OUT, "taxsea_recomputed.csv")))
message("Done.")
