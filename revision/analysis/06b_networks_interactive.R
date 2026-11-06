# Module 06b - interactive NetCoMi network exploration (run me in RStudio).
#
# WHAT THIS IS
#   A self-contained, section-by-section script for exploring the genus-rank
#   SPIEC-EASI networks and deciding how Figure 6 should look. Nothing here
#   recomputes statistics: the permutation edge tests (nPerm = 5000) were run
#   once by module 03 and are LOADED from the stored CSVs. Only the network
#   object itself is (re)built when no cache exists - construction is seeded
#   and takes ~1-3 min per partition; every later re-plot is instant.
#
# HOW TO USE (RStudio)
#   1. Open this file; use the outline (Ctrl+Shift+O) - every "# ---- " line
#      is a section you can run independently with Ctrl+Enter.
#   2. Set PARTITION in SETTINGS, run top to bottom once, then iterate on the
#      plot sections, adjusting SETTINGS as you go.
#   3. Figures you like: set SAVE_PLOTS <- TRUE and rerun a plot section; files
#      land in revision/results/figures/interactive/.
#
# DESIGN (2026-11-06)
#   "full" (default, primary): all sensitized vs all 233 controls, the same
#   samples as the paper's DA (results/networks_full). "matched": the earlier
#   1:1 matched-pair networks (results/networks), kept as the sensitivity
#   analysis. In the full design medium_score / low_score / cluster_I are
#   unbalanced against 233 controls, so their edge COUNTS reflect n; only the
#   per-edge permutation test is interpretable there.
#
# THE OVERLAY = PREVALENCE CORE (R1: "high-confidence, core taxa")
#   Core genera (present in >= 80% of the cohort) are drawn as squares. The
#   core is defined and quantified by module 03b
#   (results/networks_full/network_node_annotation.csv, core_network_stats.csv);
#   this script only reads it. It is cohort-level, so it is the same for both
#   designs. The earlier overlay (genus containing the matched-design consensus
#   DA taxon, networks/network_overlay_annotation.csv) is SUPERSEDED: that
#   consensus is not what the paper reports, and no genus is DA with a
#   calibrated method in the paper's primary analysis.

# ---- SETTINGS (edit freely) ----
PARTITION   <- "any_ige"   # any_ige | multiple | medium_score | low_score | cluster_I
# (a batch cache-warming job can override the choice without editing this file)
if (nzchar(Sys.getenv("NETCOMI_PARTITION"))) PARTITION <- Sys.getenv("NETCOMI_PARTITION")
DESIGN      <- "full"      # full (primary) | matched (sensitivity)
if (nzchar(Sys.getenv("NETCOMI_DESIGN"))) DESIGN <- Sys.getenv("NETCOMI_DESIGN")
RECOMPUTE   <- FALSE       # TRUE: rebuild the network even if a cache exists
EDGE_FILTER <- 0.00        # hide |association| below this in the plots
LABEL_CEX   <- 1.0         # node label size
SHOW_OVERLAY<- TRUE        # draw prevalence-core genera (>= 80%) as squares
SAVE_PLOTS  <- FALSE       # TRUE: also write PDFs to results/figures/interactive/
if (nzchar(Sys.getenv("NETCOMI_SAVE_PLOTS"))) SAVE_PLOTS <- TRUE   # batch rendering override
ALPHA_EDGE  <- 0.05        # q threshold when listing differential edges

# ---- 0. Setup (run once) ----
suppressPackageStartupMessages({
  library(yaml); library(NetCoMi); library(igraph)
})
# Resolve the repo root whether this runs via source(), Rscript, or RStudio.
.find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("--file=", "", a[grep("--file=", a)])
  if (length(f)) return(normalizePath(file.path(dirname(f), "..")))
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)
    if (!inherits(p, "try-error") && nzchar(p))
      return(normalizePath(file.path(dirname(p), "..")))
  }
  # fallback: assume we sit somewhere inside the repo
  normalizePath("/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision")
}
REVISION <- .find_root()
.libPaths(c(file.path(REVISION, "env", "rlibs"), .libPaths()))
cfg  <- yaml::read_yaml(file.path(REVISION, "config", "config.yaml"))
PIPE <- cfg$paths$pipeline_results
KORA <- cfg$paths$kora_data
stopifnot(DESIGN %in% c("full", "matched"))
NETD <- file.path(REVISION, "results", if (DESIGN == "full") "networks_full" else "networks")
TAXA <- file.path(REVISION, "results", "taxa_sets")
FIGD <- file.path(REVISION, "results", "figures", "interactive")
dir.create(FIGD, recursive = TRUE, showWarnings = FALSE)
SEED <- cfg$statistics$seed
PARTS_AVAILABLE <- c("any_ige", "multiple", "medium_score", "low_score", "cluster_I")
stopifnot(PARTITION %in% PARTS_AVAILABLE)
message(sprintf("repo: %s | partition: %s | design: %s", REVISION, PARTITION, DESIGN))

# ---- 1. Load the STORED edge statistics (never recomputed here) ----
edges <- read.csv(file.path(NETD, PARTITION, sprintf("edges_%s.csv", PARTITION)))
nodes <- read.csv(file.path(NETD, PARTITION, sprintf("nodes_%s.csv", PARTITION)))
message(sprintf("edges table: %d node pairs | testable: %d | q<%.2f: %d",
                nrow(edges), sum(edges$testable),
                ALPHA_EDGE, sum(edges$q_adaptBH < ALPHA_EDGE, na.rm = TRUE)))
sig_edges <- edges[!is.na(edges$q_adaptBH) & edges$q_adaptBH < ALPHA_EDGE, ]
if (nrow(sig_edges)) {
  message("differential edges surviving FDR (from the stored 5000-permutation test):")
  print(sig_edges[, c("taxon_i", "taxon_j", "asso_ige", "asso_non_ige",
                      "p_perm", "q_adaptBH")], row.names = FALSE)
} else message("no differential edge survives FDR in this partition")

# ---- 2. Prevalence-core overlay (READ from module 03b; nothing written here) ----
top50 <- read.csv(file.path(TAXA, "genus_top50_cohort.csv"), row.names = 1)
ann_f <- file.path(REVISION, "results", "networks_full", "network_node_annotation.csv")
if (!file.exists(ann_f)) stop("missing ", ann_f, " - run 03b_network_overlay.py first")
ann <- read.csv(ann_f)
ann <- ann[ann$partition == PARTITION, ]
# pandas writes booleans as "True"/"False", which read.csv leaves as character
ann$core80 <- as.logical(ann$core80)
stopifnot(!anyNA(ann$core80), nrow(ann) == 50)
core_genera <- ann$genus[ann$core80]
message(sprintf("prevalence core (>= 80%% of cohort): %d of %d genera",
                length(core_genera), nrow(ann)))
message("calibrated DA genera (LinDA/ALDEx2, q < 0.05, paper primary): ",
        if (any(ann$n_sig_calibrated > 0)) paste(ann$genus[ann$n_sig_calibrated > 0], collapse = ", ")
        else "NONE")

# ---- 3. Build or load the network object (seeded; cache = RDS) ----
cache <- file.path(NETD, PARTITION, sprintf("netprops_%s.rds", PARTITION))
if (file.exists(cache) && !RECOMPUTE) {
  props <- readRDS(cache)
  message("loaded cached netAnalyze object: ", cache)
} else {
  message("constructing network (seeded StARS; ~1-3 min) ...")
  frozen_taxa <- rownames(top50)
  src <- if (DESIGN == "full") {
    file.path(REVISION, "results", "full_cohort", "inputs", "genus", "prev10", PARTITION, "source")
  } else file.path(PIPE, "genus", "prev10", PARTITION, "source")
  oi <- read.csv(file.path(src, "otu_ige.csv"),     row.names = 1, check.names = FALSE)
  oc <- read.csv(file.path(src, "otu_non_ige.csv"), row.names = 1, check.names = FALSE)
  mi <- read.csv(file.path(src, "sample_data_ige.csv"),     row.names = 1)
  mc <- read.csv(file.path(src, "sample_data_non_ige.csv"), row.names = 1)
  if (DESIGN == "matched") {     # same pair alignment as module 03
    common_pairs <- sort(intersect(mi$pair_nb, mc$pair_nb))
    mi <- mi[mi$pair_nb %in% common_pairs, , drop = FALSE]
    mc <- mc[mc$pair_nb %in% common_pairs, , drop = FALSE]
    mi <- mi[order(match(mi$pair_nb, common_pairs)), , drop = FALSE]
    mc <- mc[order(match(mc$pair_nb, common_pairs)), , drop = FALSE]
  }
  avail <- intersect(frozen_taxa, intersect(colnames(oi), colnames(oc)))
  otu <- as.matrix(rbind(oi[rownames(mi), avail, drop = FALSE],
                         oc[rownames(mc), avail, drop = FALSE]))
  grp <- c(rep("ige", nrow(mi)), rep("non_ige", nrow(mc)))
  # IDENTICAL construction arguments to module 03 (same seed -> same network);
  # see 03_genus_networks.R for the reasoning behind every value.
  net <- netConstruct(data = otu, group = grp,
    measure = "spieceasi",
    measurePar = list(method = "glasso", nlambda = 30, lambda.min.ratio = 0.01,
                      sel.criterion = "stars",
                      pulsar.params = list(rep.num = 20, thresh = 0.05,
                                           subsample.ratio = 0.8, ncores = 1)),
    filtTax = "none", filtTaxPar = NULL, filtSamp = "none", filtSampPar = NULL,
    zeroMethod = "none", zeroPar = NULL, normMethod = "none", normPar = NULL,
    sparsMethod = "none", dissFunc = "signed", weighted = TRUE,
    jointPrepro = FALSE,
    matchDesign = if (DESIGN == "matched") c(1, 1) else NULL, seed = SEED,
    logFile = NULL, verbose = 0)
  props <- netAnalyze(net,
    centrLCC = FALSE, weightDeg = FALSE, normDeg = TRUE, normBetw = TRUE,
    normClose = TRUE, normEigen = TRUE,
    clustMethod = "cluster_fast_greedy", hubPar = "eigenvector", hubQuant = 0.95,
    lnormFit = FALSE, graphlet = FALSE, avDissIgnoreInf = FALSE,
    sPathAlgo = "dijkstra", sPathNorm = TRUE, normNatConnect = TRUE,
    weightClustCoef = TRUE, verbose = 0)
  saveRDS(props, cache)
  message("cached: ", cache)
  # QC: the reconstruction must match the stored edge tables (same seed).
  n_e1 <- sum(net$adjaMat1[upper.tri(net$adjaMat1)] != 0)
  n_e2 <- sum(net$adjaMat2[upper.tri(net$adjaMat2)] != 0)
  stored1 <- sum(edges$present_ige); stored2 <- sum(edges$present_non_ige)
  message(sprintf("QC edges rebuilt vs stored: IgE %d vs %d | non-IgE %d vs %d %s",
                  n_e1, stored1, n_e2, stored2,
                  if (n_e1 == stored1 && n_e2 == stored2) "(MATCH)"
                  else "(MISMATCH - investigate before trusting plots)"))
}

# ---- 4. PLOT: the two association networks side by side ----
node_names <- rownames(props$input$adjaMat1)
# core genera as squares via NetCoMi's shape-by-feature mechanism
shape_vec <- if (SHOW_OVERLAY) factor(ifelse(node_names %in% core_genera, "core", "non-core"),
                                      levels = c("core", "non-core")) else NULL
if (!is.null(shape_vec)) names(shape_vec) <- node_names
if (SAVE_PLOTS) pdf(file.path(FIGD, sprintf("net_%s_%s.pdf", DESIGN, PARTITION)),
                    width = 13, height = 7)
plot(props,
     sameLayout = TRUE, layoutGroup = "union",
     nodeColor = "cluster", nodeSize = "eigenvector",
     labelScale = FALSE, cexLabels = LABEL_CEX,
     nodeShape = if (SHOW_OVERLAY) c("square", "circle") else "circle",
     featVecShape = shape_vec,
     edgeFilter = if (EDGE_FILTER > 0) "threshold" else "none",
     edgeFilterPar = EDGE_FILTER,
     highlightHubs = TRUE,
     nodeTransp = 20, hubTransp = 0,
     groupNames = c(sprintf("IgE-sensitized (%s, %s)", PARTITION, DESIGN), "non-sensitized"),
     showTitle = TRUE, cexTitle = 1.1)
if (SHOW_OVERLAY)
  legend("bottomleft", bty = "n", cex = 0.85, pch = c(15, 19),
         legend = c(sprintf("prevalence core (>= 80%%; %d genera)", length(core_genera)),
                    "non-core"))
if (SAVE_PLOTS) dev.off()

# ---- 5. PLOT: differential edges (from the STORED permutation test) ----
de <- edges[edges$testable & !is.na(edges$p_perm), ]
de$sig <- !is.na(de$q_adaptBH) & de$q_adaptBH < ALPHA_EDGE
g <- graph_from_data_frame(
  de[de$present_ige | de$present_non_ige,
     c("taxon_i", "taxon_j", "asso_ige", "asso_non_ige", "q_adaptBH", "sig")],
  directed = FALSE)
E(g)$weight <- abs(E(g)$asso_ige - E(g)$asso_non_ige)
if (SAVE_PLOTS) pdf(file.path(FIGD, sprintf("diffnet_%s_%s.pdf", DESIGN, PARTITION)),
                    width = 9, height = 9)
set.seed(SEED)
plot(g,
     edge.width = 1 + 6 * E(g)$weight / max(E(g)$weight),
     edge.color = ifelse(E(g)$sig, "#C0392B", "grey80"),
     vertex.size = 6, vertex.color = "grey90", vertex.frame.color = "grey50",
     vertex.label.cex = LABEL_CEX * 0.8, vertex.label.color = "black",
     main = sprintf("%s: edge-wise |difference|; red = q < %.2f (nPerm = 5000)",
                    PARTITION, ALPHA_EDGE))
if (SAVE_PLOTS) dev.off()

# ---- 6. PLOT: association agreement between the two groups ----
pe <- edges[edges$present_ige | edges$present_non_ige, ]
if (SAVE_PLOTS) pdf(file.path(FIGD, sprintf("edge_scatter_%s_%s.pdf", DESIGN, PARTITION)),
                    width = 6.5, height = 6)
plot(pe$asso_non_ige, pe$asso_ige, pch = 19, cex = 0.7,
     col = ifelse(!is.na(pe$q_adaptBH) & pe$q_adaptBH < ALPHA_EDGE,
                  "#C0392B", "#00000055"),
     xlab = "association, non-sensitized", ylab = "association, IgE-sensitized",
     main = sprintf("%s: edge associations (red = differential at q<%.2f)",
                    PARTITION, ALPHA_EDGE))
abline(0, 1, lty = 2, col = "grey60")
if (SAVE_PLOTS) dev.off()

# ---- 7. Export the combined supplementary edge table (all 5 partitions) ----
supp <- do.call(rbind, lapply(PARTS_AVAILABLE, function(p)
  read.csv(file.path(NETD, p, sprintf("edges_%s.csv", p)))))
write.csv(supp, file.path(NETD, "network_edges_supplementary.csv"), row.names = FALSE)
message(sprintf("wrote %s/network_edges_supplementary.csv (%d rows, %d partitions)", basename(NETD),
                nrow(supp), length(PARTS_AVAILABLE)))
message("Done - iterate on the SETTINGS block and rerun sections 4-6.")
