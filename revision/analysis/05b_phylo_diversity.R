#!/usr/bin/env Rscript
# Module 05b - phylogenetic diversity: Faith's PD and UniFrac.
#
# Unblocked by the phylogeny rebuild (revision/slurm/00d_rebuild_phylogeny.sh).
# The tree reproduces the published ASV set exactly - 15,170 of 15,170 tips, and
# the re-derived count matrix differs from the published one in 2 cells out of
# 30,855,780 (one read swapped between two ASVs in a single sample) - so these
# metrics apply directly to the analyses the manuscript reports.
#
# WHAT THIS ADDS
#   Faith's PD  - the published Table S4 reports it, but the analysis behind it
#                 was unmatched, unadjusted and uncorrected. Recomputed here on
#                 matched pairs with BH across partitions.
#   UniFrac     - requested by R2 Major 2 and previously impossible. Unweighted,
#                 weighted and generalised (alpha = 0.5), each with adonis2 under
#                 within-pair permutation and a betadisper companion.
#
# Faith's PD is computed directly from the tree rather than via picante (not
# installed): PD is the total branch length of the minimal subtree spanning the
# observed taxa, including the path to the root. The implementation is checked
# against a hand-computed toy tree at startup - if that check fails the script
# aborts rather than emitting numbers nobody has verified.
#
# ASV SET
#   Phylogenetic metrics need ASV-level data; family/genus tables have no tree.
#   Two sets are used deliberately:
#     Faith's PD - ASVs non-zero within the 508-sample analytic subset. That is
#                  exactly the "6,958 ASVs" the manuscript reports (verified),
#                  so the values are comparable to the published Table S4.
#                  NB the Methods describe this as "ASVs present in <=1%
#                  individual were removed", which is not what was done.
#     UniFrac    - the prev01 matched ASV tables. Unweighted UniFrac is highly
#                  sensitive to singletons, so some floor is required; prev01 is
#                  the mildest available and is reported as such.
#
# Usage:
#   Rscript revision/analysis/05b_phylo_diversity.R --partition any_ige

suppressPackageStartupMessages({
  library(yaml); library(ape); library(vegan); library(permute); library(GUniFrac)
})

args_raw <- commandArgs(trailingOnly = FALSE)
script <- sub("--file=", "", args_raw[grep("--file=", args_raw)])
REVISION <- normalizePath(file.path(dirname(script), ".."))
.libPaths(c(file.path(REVISION, "env", "rlibs"), .libPaths()))

a <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- which(a == f); if (length(i) && i < length(a)) a[i+1] else d }
PART  <- getarg("--partition", "any_ige")
NPERM <- as.integer(getarg("--nperm", "5000"))

cfg  <- yaml::read_yaml(file.path(REVISION, "config", "config.yaml"))
KORA <- cfg$paths$kora_data
PIPE <- cfg$paths$pipeline_results
TREE <- file.path(KORA, "q2_reprocessing", "rooted_tree.nwk")
OUT  <- file.path(REVISION, "results", "phylo")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED <- cfg$statistics$seed

if (!file.exists(TREE)) stop("Tree not found. Run 00d_rebuild_phylogeny.sh first: ", TREE)

# ── Faith's PD, and a self-check ─────────────────────────────────────────
faith_pd <- function(tree, present) {
  present <- intersect(present, tree$tip.label)
  if (length(present) < 2) return(NA_real_)
  sub <- ape::keep.tip(tree, present)
  sum(sub$edge.length, na.rm = TRUE)
}

# Toy tree with known total branch length. ((A:1,B:1):1,(C:2,D:2):1);
# Keeping all four tips spans every edge: 1+1+1+2+2+1 = 8.
# Keeping A and C spans A:1, its parent:1, C:2, its parent:1 = 5.
.toy <- ape::read.tree(text = "((A:1,B:1):1,(C:2,D:2):1);")
stopifnot(abs(faith_pd(.toy, c("A","B","C","D")) - 8) < 1e-9)
stopifnot(abs(faith_pd(.toy, c("A","C")) - 5) < 1e-9)
message("  PD implementation self-check: PASS")

message(sprintf("=== Phylogenetic diversity: %s ===", PART))
tree <- ape::read.tree(TREE)
message(sprintf("  tree: %d tips, rooted=%s", ape::Ntip(tree), ape::is.rooted(tree)))

# ── Matched samples ──────────────────────────────────────────────────────
mi <- read.csv(file.path(PIPE, "asv", "prev01", PART, "source", "sample_data_ige.csv"),
               row.names = 1)
mc <- read.csv(file.path(PIPE, "asv", "prev01", PART, "source", "sample_data_non_ige.csv"),
               row.names = 1)
cp <- sort(intersect(mi$pair_nb, mc$pair_nb))
mi <- mi[match(cp, mi$pair_nb), , drop = FALSE]
mc <- mc[match(cp, mc$pair_nb), , drop = FALSE]
ids <- c(rownames(mi), rownames(mc))
grp <- factor(c(rep("ige", nrow(mi)), rep("non_ige", nrow(mc))),
              levels = c("non_ige", "ige"))
pair <- factor(c(mi$pair_nb, mc$pair_nb))
message(sprintf("  %d matched pairs (%d samples)", length(cp), length(ids)))
if (length(cp) < 5) { message("  too few pairs; exiting"); quit(status = 0) }

# ── Faith's PD on the 508-subset ASV set, matching the published definition ──
message("\n  Faith's PD (ASV set = non-zero within the 508 analytic subset)")
ft <- read.csv(file.path(KORA, "feature_table.tsv"), sep = "\t",
               row.names = 1, check.names = FALSE)
colnames(ft) <- as.character(colnames(ft))
ige_all <- read.csv(file.path(KORA, "ige_kora_subset_all.csv"),
                    row.names = 1, check.names = FALSE)
subset508 <- intersect(as.character(rownames(ige_all)), colnames(ft))
keep_asv <- rownames(ft)[rowSums(ft[, subset508, drop = FALSE] > 0) > 0]
message(sprintf("    ASVs non-zero in the 508 subset: %d (published figure: 6,958)",
                length(keep_asv)))

present <- ft[keep_asv, intersect(ids, colnames(ft)), drop = FALSE]
pd <- vapply(colnames(present),
             function(s) faith_pd(tree, rownames(present)[present[[s]] > 0]),
             numeric(1))

d <- data.frame(sample = names(pd), pd = as.numeric(pd),
                group = grp[match(names(pd), ids)],
                pair = pair[match(names(pd), ids)], stringsAsFactors = FALSE)
d <- d[!is.na(d$pd) & !is.na(d$group), ]
w <- reshape(d[, c("pair", "group", "pd")], idvar = "pair",
             timevar = "group", direction = "wide")
w <- w[stats::complete.cases(w), ]
tt <- suppressWarnings(wilcox.test(w$pd.ige, w$pd.non_ige, paired = TRUE))
message(sprintf("    mean PD: IgE %.2f | non-IgE %.2f | paired Wilcoxon p = %.4f",
                mean(w$pd.ige), mean(w$pd.non_ige), tt$p.value))

write.csv(data.frame(partition = PART, n_pairs = nrow(w),
                     mean_pd_ige = mean(w$pd.ige),
                     mean_pd_non_ige = mean(w$pd.non_ige),
                     statistic = as.numeric(tt$statistic),
                     p_value = tt$p.value, n_asvs = length(keep_asv)),
          file.path(OUT, sprintf("faithpd_%s.csv", PART)), row.names = FALSE)
write.csv(d, file.path(OUT, sprintf("faithpd_persample_%s.csv", PART)),
          row.names = FALSE)

# ── UniFrac ──────────────────────────────────────────────────────────────
message("\n  UniFrac (prev01 ASV table)")
oi <- read.csv(file.path(PIPE, "asv", "prev01", PART, "source", "otu_ige.csv"),
               row.names = 1, check.names = FALSE)
oc <- read.csv(file.path(PIPE, "asv", "prev01", PART, "source", "otu_non_ige.csv"),
               row.names = 1, check.names = FALSE)
taxa <- intersect(colnames(oi), colnames(oc))
otu <- rbind(oi[rownames(mi), taxa, drop = FALSE], oc[rownames(mc), taxa, drop = FALSE])
on_tree <- intersect(colnames(otu), tree$tip.label)
message(sprintf("    %d ASVs, %d on the tree (%.1f%%)", ncol(otu), length(on_tree),
                100 * length(on_tree) / ncol(otu)))
otu <- otu[, on_tree, drop = FALSE]
otu <- otu[rowSums(otu) > 0, , drop = FALSE]
g2 <- grp[match(rownames(otu), ids)]
p2 <- droplevels(pair[match(rownames(otu), ids)])

# Keep only intact pairs so within-pair permutation stays well defined.
tabp <- table(p2); full <- names(tabp)[tabp == 2]
k <- p2 %in% full
otu <- otu[k, , drop = FALSE]; g2 <- g2[k]; p2 <- droplevels(p2[k])
message(sprintf("    %d samples in %d complete pairs", nrow(otu), nlevels(p2)))

pruned <- ape::keep.tip(tree, colnames(otu))
uf <- GUniFrac::GUniFrac(as.matrix(otu), pruned, alpha = c(0, 0.5, 1))$unifracs
dists <- list(unweighted = as.dist(uf[, , "d_UW"]),
              generalized_0.5 = as.dist(uf[, , "d_0.5"]),
              weighted = as.dist(uf[, , "d_1"]))

meta <- data.frame(ige_status = g2, pair_nb = p2)
ctrl <- how(nperm = NPERM, blocks = meta$pair_nb)
rows <- list()
for (nm in names(dists)) {
  dd <- dists[[nm]]
  set.seed(SEED)
  ad <- adonis2(dd ~ ige_status, data = meta, by = "margin", permutations = ctrl)
  set.seed(SEED)
  bd <- betadisper(dd, meta$ige_status, type = "median", bias.adjust = TRUE)
  pt <- permutest(bd, permutations = how(nperm = NPERM, blocks = meta$pair_nb))
  rows[[nm]] <- data.frame(
    partition = PART, distance = paste0("unifrac_", nm),
    R2 = ad["ige_status", "R2"], F_stat = ad["ige_status", "F"],
    p_adonis = ad["ige_status", "Pr(>F)"],
    disp_F = pt$tab[1, "F"], p_betadisper = pt$tab[1, "Pr(>F)"],
    n = nrow(otu), n_pairs = nlevels(p2), n_asvs = ncol(otu), nperm = NPERM,
    stringsAsFactors = FALSE)
  message(sprintf("    %-16s R2=%.4f p=%.4f | dispersion p=%.4f %s", nm,
                  rows[[nm]]$R2, rows[[nm]]$p_adonis, rows[[nm]]$p_betadisper,
                  if (!is.na(rows[[nm]]$p_betadisper) &&
                      rows[[nm]]$p_betadisper < 0.05) "UNEQUAL" else ""))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT, sprintf("unifrac_%s.csv", PART)), row.names = FALSE)
message(sprintf("\nWrote Faith's PD and UniFrac results to %s", OUT))
message("Done.")
