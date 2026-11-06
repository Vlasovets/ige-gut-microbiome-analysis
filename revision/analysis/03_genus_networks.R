#!/usr/bin/env Rscript
# Module 03 - genus-level association networks with a real differential test.
#
# Answers R1 Major 2 / R2 Minor 11 & 13. Four things the published analysis did
# that this replaces:
#
# 1. NO STATISTICAL TEST. It used diffnet(diffMethod = "threshold",
#    diffThresh = 0.0, adjust = "none"), which declares every edge whose two
#    association estimates differ AT ALL to be "differential", while the notebook
#    printed "Only edges showing significant differences are displayed". Replaced
#    with diffMethod = "permute" and adaptive-BH adjustment, plus netCompare() for
#    global properties (netCompare was never called at all).
#
# 2. SAMPLE-SIZE CONFOUNDING. All comparisons passed the full 233-control table,
#    so cluster networks at n = 22-24 were compared against n = 233. Edge count is
#    monotone in n at fixed p (measured: 6 edges at n=23 ... 50 at n=221), so the
#    reported "loss of edges in IgE" was confounded with sample size. Here the
#    matched 1:1 design is used and, crucially, matchDesign = c(1,1) makes the
#    permutation test swap labels WITHIN pairs.
#
# 3. AN INERT FILTER, APPLIED PER COMPARISON. filtTax = "highestVar" with
#    highestVar = 50 removed nothing (only 41-43 family taxa existed) and, because
#    it ran inside each netConstruct call, the control network was re-estimated on
#    a different taxa subset each time. Here a FROZEN top-50 genus manifest from
#    00c_define_taxa_sets.py is used, so every network has identical nodes.
#
# 4. NO EXPORTED STATISTICS. The published Rmd contains no write.csv at all, so no
#    edge weights, effect sizes, p-values or FDR ever left the session. This module
#    exports all of them.
#
# Every netConstruct/netAnalyze argument is passed explicitly, including ones left
# at defaults, so the parameter table in the response letter is complete.
#
# Two designs (--design):
#   matched  1:1 matched pairs (pipeline source tables), matchDesign = c(1,1) so
#            permutations swap within pairs. Outputs: results/networks/<part>/.
#   full     all sensitized members of the partition vs all 233 controls, the
#            sample set of the primary DA (02f inputs). Unpaired permutations.
#            Outputs: results/networks_full/<part>/. SPIEC-EASI cannot take
#            covariates, so unlike the DA this network is NOT adjusted for
#            sex/age/BMI. Edge count grows with n, so when n_ige and n_ctrl
#            differ a lot (medium/low_score, clusters) only the per-edge
#            permutation test is interpretable, not edge counts or global props.
#
# Usage:
#   Rscript revision/analysis/03_genus_networks.R --partition any_ige --nperm 1000
#   Rscript revision/analysis/03_genus_networks.R --partition any_ige --nperm 1000 --design full

suppressPackageStartupMessages({
  library(yaml); library(NetCoMi); library(SpiecEasi); library(pulsar)
})

args_raw <- commandArgs(trailingOnly = FALSE)
script <- sub("--file=", "", args_raw[grep("--file=", args_raw)])
REVISION <- normalizePath(file.path(dirname(script), ".."))
.libPaths(c(file.path(REVISION, "env", "rlibs"), .libPaths()))

a <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- which(a == f); if (length(i) && i < length(a)) a[i+1] else d }
PART  <- getarg("--partition", "any_ige")
NPERM <- as.integer(getarg("--nperm", "1000"))
TOPN  <- as.integer(getarg("--top-n", "50"))
CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
DESIGN <- getarg("--design", "matched")
stopifnot(DESIGN %in% c("matched", "full"))

cfg  <- yaml::read_yaml(file.path(REVISION, "config", "config.yaml"))
PIPE <- cfg$paths$pipeline_results
NET  <- cfg$networks
LEVEL <- NET$level; PREV <- NET$prevalence
SEED <- cfg$statistics$seed

TAXA <- file.path(REVISION, "results", "taxa_sets")
OUT  <- file.path(REVISION, "results",
                  if (DESIGN == "full") "networks_full" else "networks", PART)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

message(sprintf("=== Genus networks: %s, %s design (top %d, nperm=%d, cores=%d) ===",
                PART, DESIGN, TOPN, NPERM, CORES))

# ── Frozen taxon set ─────────────────────────────────────────────────────
mf <- file.path(TAXA, sprintf("%s_top%d_cohort.csv", LEVEL, TOPN))
if (!file.exists(mf)) stop("Run 00c_define_taxa_sets.py first: ", mf)
frozen <- read.csv(mf, row.names = 1)
frozen_taxa <- rownames(frozen)
tt <- read.csv(file.path(TAXA, sprintf("tax_table_%s.csv", LEVEL)), row.names = 1)

# ── Data: matched pairs (aligned across blocks) or the full cohort ─────────
src <- if (DESIGN == "full") {
  file.path(REVISION, "results", "full_cohort", "inputs", LEVEL, PREV, PART, "source")
} else file.path(PIPE, LEVEL, PREV, PART, "source")
if (!dir.exists(src)) stop("missing inputs (full design: run 02f first): ", src)
oi <- read.csv(file.path(src, "otu_ige.csv"),     row.names = 1, check.names = FALSE)
oc <- read.csv(file.path(src, "otu_non_ige.csv"), row.names = 1, check.names = FALSE)
mi <- read.csv(file.path(src, "sample_data_ige.csv"),     row.names = 1)
mc <- read.csv(file.path(src, "sample_data_non_ige.csv"), row.names = 1)

if (DESIGN == "matched") {
  # Keep only pairs present in both blocks, then order both by pair_nb. matchDesign
  # requires the i-th sample of group 1 to correspond to the i-th of group 2.
  common_pairs <- sort(intersect(mi$pair_nb, mc$pair_nb))
  mi <- mi[order(match(mi$pair_nb, common_pairs)), , drop = FALSE]
  mc <- mc[order(match(mc$pair_nb, common_pairs)), , drop = FALSE]
  mi <- mi[mi$pair_nb %in% common_pairs, , drop = FALSE]
  mc <- mc[mc$pair_nb %in% common_pairs, , drop = FALSE]
  stopifnot(identical(mi$pair_nb, mc$pair_nb))
} else {
  stopifnot(setequal(rownames(mi), rownames(oi)), setequal(rownames(mc), rownames(oc)),
            !length(intersect(rownames(oi), rownames(oc))))
}

avail <- intersect(frozen_taxa, intersect(colnames(oi), colnames(oc)))
message(sprintf("  frozen taxa: %d | available in this partition: %d",
                length(frozen_taxa), length(avail)))
if (length(avail) < length(frozen_taxa))
  message(sprintf("    missing: %s",
                  paste(setdiff(frozen_taxa, avail), collapse = ", ")))

# netConstruct requires a matrix; rbind of data.frames yields a data.frame and
# fails with "data must inherit from one of the classes: matrix, ...".
otu <- as.matrix(rbind(oi[rownames(mi), avail, drop = FALSE],
                       oc[rownames(mc), avail, drop = FALSE]))
grp <- c(rep("ige", nrow(mi)), rep("non_ige", nrow(mc)))
N_IGE <- nrow(mi); N_CTRL <- nrow(mc)
n_pairs <- if (DESIGN == "matched") length(common_pairs) else NA_integer_
message(sprintf("  n=%d (ige %d, non-ige %d; %s), taxa=%d", nrow(otu), N_IGE, N_CTRL,
                if (DESIGN == "matched") sprintf("%d complete pairs", n_pairs)
                else "unpaired", ncol(otu)))

if (min(N_IGE, N_CTRL) < 15) {
  message("  smaller group < 15: graphical-lasso edge selection is not identifiable here.")
  message("  Writing a stub and exiting; report this partition as exploratory only.")
  write.csv(data.frame(partition = PART, design = DESIGN, n_ige = N_IGE, n_ctrl = N_CTRL,
                       note = "skipped: smaller group < 15, network not identifiable"),
            file.path(OUT, "SKIPPED.csv"), row.names = FALSE)
  quit(status = 0)
}

# ── Network construction: EVERY argument explicit ─────────────────────────
PARAMS <- list(
  measure = "spieceasi",
  measurePar = list(method = "glasso", nlambda = 30, lambda.min.ratio = 0.01,
                    sel.criterion = "stars",
                    # ncores MUST stay 1: pulsar's parallel workers cannot see
                    # SpiecEasi::sparseiCov and fail with "object 'sparseiCov' of
                    # mode 'function' was not found". Serial is fine here anyway -
                    # one network takes ~6 s at p=50.
                    pulsar.params = list(rep.num = 20, thresh = 0.05,
                                         subsample.ratio = 0.8, ncores = 1)),
  # A frozen manifest already fixed the taxa, so no in-call filtering.
  filtTax = "none", filtTaxPar = NULL, filtSamp = "none", filtSampPar = NULL,
  zeroMethod = "none", zeroPar = NULL, normMethod = "none", normPar = NULL,
  # SPIEC-EASI + StARS already selects sparsity; a second sparsification would
  # double-filter. This is why sparsMethod is "none" and not the "t-test" default.
  sparsMethod = "none",
  dissFunc = "signed", weighted = TRUE, jointPrepro = FALSE,
  # matched: permutations swap labels within pairs; full: unpaired permutations
  matchDesign = if (DESIGN == "matched") c(1, 1) else NULL,
  seed = SEED, cores = CORES,
  logFile = NULL,             # <- default is "log.txt"; array tasks would race
  verbose = 0
)
message("\n  Constructing network (SPIEC-EASI / glasso / StARS) ...")
t0 <- Sys.time()
net <- do.call(netConstruct, c(list(data = otu, group = grp), PARAMS))
message(sprintf("    done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
message(sprintf("    matchDesign stored on the object: %s", if (is.null(net$matchDesign)) "none (unpaired)" else
                paste(net$matchDesign, collapse = ",")))

n_edge <- function(m) sum(m[upper.tri(m)] != 0)
message(sprintf("    edges: IgE=%d  non-IgE=%d  (of %d possible)",
                n_edge(net$adjaMat1), n_edge(net$adjaMat2),
                choose(ncol(otu), 2)))

# ── Analysis: explicit arguments, not silent defaults ────────────────────
message("  Analyzing ...")
props <- netAnalyze(net,
  centrLCC = FALSE, weightDeg = FALSE, normDeg = TRUE, normBetw = TRUE,
  normClose = TRUE, normEigen = TRUE,
  clustMethod = "cluster_fast_greedy", hubPar = "eigenvector", hubQuant = 0.95,
  lnormFit = FALSE, graphlet = FALSE, avDissIgnoreInf = FALSE,
  sPathAlgo = "dijkstra", sPathNorm = TRUE, normNatConnect = TRUE,
  weightClustCoef = TRUE, verbose = 0)

# ── Global comparison with a permutation test ────────────────────────────
# storeAssoPerm lets diffnet reuse these permuted fits, halving total compute.
perm_file <- file.path(OUT, "assoPerm")
message(sprintf("  netCompare with %d permutations (this dominates runtime) ...", NPERM))
t0 <- Sys.time()
cmp <- tryCatch(
  netCompare(props, permTest = TRUE, nPerm = NPERM, adjust = "adaptBH",
             testRand = FALSE, gcd = FALSE, cores = CORES, seed = SEED,
             storeAssoPerm = TRUE, fileStoreAssoPerm = perm_file,
             logFile = NULL, verbose = FALSE),
  error = function(e) { message("    netCompare ERROR: ", conditionMessage(e)); NULL })
message(sprintf("    done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ── Per-edge differential test ───────────────────────────────────────────
message(sprintf("  diffnet (permute, %d perms, adaptBH) ...", NPERM))
dn <- tryCatch(
  diffnet(net, diffMethod = "permute", nPerm = NPERM, fisherTrans = TRUE,
          permPvalsMethod = "pseudo", adjust = "adaptBH", alpha = 0.05,
          cores = CORES, seed = SEED, logFile = NULL, verbose = FALSE,
          fileLoadAssoPerm = if (!is.null(cmp)) perm_file else NULL),
  error = function(e) { message("    diffnet ERROR: ", conditionMessage(e)); NULL })

# ── Export edge statistics — none of this existed before ─────────────────
taxa_n <- colnames(net$adjaMat1)
ut <- which(upper.tri(matrix(0, length(taxa_n), length(taxa_n))), arr.ind = TRUE)
edges <- data.frame(
  partition = PART, design = DESIGN, level = LEVEL,
  n_ige = N_IGE, n_ctrl = N_CTRL, n_pairs = n_pairs,
  taxon_i = taxa_n[ut[, 1]], taxon_j = taxa_n[ut[, 2]],
  asso_ige     = net$assoMat1[ut],
  asso_non_ige = net$assoMat2[ut],
  adja_ige     = net$adjaMat1[ut],
  adja_non_ige = net$adjaMat2[ut],
  stringsAsFactors = FALSE)
edges$difference <- edges$asso_ige - edges$asso_non_ige
edges$present_ige     <- edges$adja_ige != 0
edges$present_non_ige <- edges$adja_non_ige != 0
edges$phylum_i <- tt[edges$taxon_i, "Phylum"]
edges$phylum_j <- tt[edges$taxon_j, "Phylum"]

if (!is.null(dn)) {
  key <- paste(edges$taxon_i, edges$taxon_j)
  getmat <- function(m) {
    if (is.null(m)) return(rep(NA_real_, nrow(edges)))
    v <- rep(NA_real_, nrow(edges))
    cx <- intersect(rownames(m), taxa_n)
    mm <- matrix(NA_real_, length(taxa_n), length(taxa_n),
                 dimnames = list(taxa_n, taxa_n))
    mm[cx, cx] <- m[cx, cx]
    mm[ut]
  }
  edges$p_perm <- getmat(dn$pvalsVec_mat)
  if (all(is.na(edges$p_perm)) && !is.null(dn$pAdjustVec))
    edges$p_perm <- getmat(dn$pvalsMat)
  edges$q_adaptBH_all <- getmat(dn$pAdjustMat)

  # ── Restrict the multiplicity to testable edges, and check the floor ─────
  # Two things matter here.
  #
  # (a) Correcting over all choose(p,2) = 1225 node pairs is needlessly harsh:
  #     an edge absent from BOTH networks carries no hypothesis. Restricting to
  #     edges present in at least one group is standard and is pre-specified.
  #
  # (b) A permutation p-value has a FLOOR of 1/(nPerm+1). With m tests, BH can
  #     only ever reach q < alpha if floor * m / 1 < alpha, i.e.
  #         nPerm > m/alpha - 1.
  #     Otherwise the test is structurally incapable of a positive result no
  #     matter how large the true effect. This is easy to miss and would make a
  #     null look like evidence of absence. It is asserted below.
  edges$testable <- edges$present_ige | edges$present_non_ige
  m <- sum(edges$testable)
  # NB: despite the column name, this is PLAIN Benjamini-Hochberg over the
  # testable edges; diffnet's own adaptive-BH values (over all pairs) are kept
  # in q_adaptBH_all. Name kept for compatibility with downstream modules.
  edges$q_adaptBH <- NA_real_
  edges$q_adaptBH[edges$testable] <-
    p.adjust(edges$p_perm[edges$testable], method = "BH")

  floor_p <- 1 / (NPERM + 1)
  min_q_possible <- floor_p * m
  message(sprintf("\n    testable edges (present in >=1 group): %d of %d", m,
                  nrow(edges)))
  message(sprintf("    permutation p floor = 1/(%d+1) = %.5f", NPERM, floor_p))
  message(sprintf("    smallest attainable BH q over %d tests = %.4f", m,
                  min_q_possible))
  if (min_q_possible >= 0.05) {
    message(sprintf(paste0("    *** UNDERPOWERED BY CONSTRUCTION: no edge can reach ",
                           "q<0.05.\n        Need nPerm > %d. Re-run before drawing any ",
                           "conclusion."), ceiling(m / 0.05) - 1))
  } else {
    message("    resolution is adequate: a real edge difference could be detected.")
  }

  n_sig <- sum(edges$q_adaptBH < 0.05, na.rm = TRUE)
  message(sprintf("    differential edges at q<0.05: %d of %d testable", n_sig, m))
  if (n_sig == 0 && min_q_possible < 0.05)
    message("    -> No individual edge survives FDR control, and the test HAD the\n",
            "       resolution to find one. Pre-committed as a result: consistent\n",
            "       with the ~0.5% compositional effect size.")
} else {
  edges$p_perm <- NA_real_; edges$q_adaptBH <- NA_real_
}

write.csv(edges, file.path(OUT, sprintf("edges_%s.csv", PART)), row.names = FALSE)

# Global properties table.
if (!is.null(cmp)) {
  gl <- tryCatch({
    s <- summary(cmp)
    d <- as.data.frame(s$globalProps)
    d$property <- rownames(d); d$partition <- PART; d
  }, error = function(e) NULL)
  if (!is.null(gl)) {
    write.csv(gl, file.path(OUT, sprintf("global_props_%s.csv", PART)),
              row.names = FALSE)
    message("\n  Global network properties (both groups, permutation p):")
    print(head(gl, 12), row.names = FALSE)
  }
}

# Node-level properties for the figure.
nodes <- data.frame(
  partition = PART, design = DESIGN, taxon = taxa_n,
  phylum = tt[taxa_n, "Phylum"],
  degree_ige     = props$centralities$degree1[taxa_n],
  degree_non_ige = props$centralities$degree2[taxa_n],
  clust_ige      = props$clustering$clust1[taxa_n],
  clust_non_ige  = props$clustering$clust2[taxa_n],
  stringsAsFactors = FALSE)
write.csv(nodes, file.path(OUT, sprintf("nodes_%s.csv", PART)), row.names = FALSE)

# Parameter record for the STAR Methods table.
writeLines(c(
  sprintf("partition: %s", PART),
  sprintf("level: %s   prevalence: %s   top_n: %d", LEVEL, PREV, TOPN),
  sprintf("design: %s   n: %d (ige %d, non-ige %d)   n_pairs: %s   taxa: %d",
          DESIGN, nrow(otu), N_IGE, N_CTRL, n_pairs, ncol(otu)),
  sprintf("nPerm: %d   cores: %d   seed: %d", NPERM, CORES, SEED),
  "netConstruct arguments (all explicit):",
  paste0("  ", capture.output(str(PARAMS, max.level = 2))),
  "",
  "SpiecEasi version: ", as.character(packageVersion("SpiecEasi")),
  "NetCoMi version: ", as.character(packageVersion("NetCoMi")),
  "", capture.output(sessionInfo())
), file.path(OUT, sprintf("params_%s.txt", PART)))

message(sprintf("\nWrote edge/node/global tables to %s", OUT))
message("Done.")
