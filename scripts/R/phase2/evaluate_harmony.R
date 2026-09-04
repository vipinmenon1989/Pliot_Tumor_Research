# scripts/R/phase2/evaluate_harmony.R
#
# Phase 2 / Milestone M12 - Pre/post Harmony integration evaluation.
#
# Question: did default Harmony (M11) reduce unwanted sample/dataset-associated
# structure WHILE preserving biologically meaningful structure?
#
# Compares, on identical cells with identical parameters and the same seed:
#     PRE  : reduction 'pca'                dims 1:30   (frozen Phase 1 baseline)
#     POST : reduction 'postint_harmony'    dims 1:30   (M11 output)
#
# Nothing is tuned. Harmony is NOT re-run. No clustering, no FindMarkers, no
# annotation. The M11 and Phase 1 objects are opened read-only and never rewritten.
#
# The Harmony UMAP required for the paired figures is computed in memory and only
# its 2-column embedding is written out; no third multi-gigabyte object is created.

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = +Inf)

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(RANN); library(cluster)
  library(ggplot2); library(patchwork); library(viridis); library(RColorBrewer)
  library(jsonlite); library(digest)
})

source("scripts/R/logging_utils.R")
source("scripts/R/provenance_utils.R")
setup_strict_logging()

STAGE <- "phase2_m12_harmony_evaluation"

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_rds     <- "results/phase2/harmony/phase2_harmony_integrated.rds"
out_dir       <- "results/phase2/harmony/evaluation"
fig_dir       <- "reports/phase2/figures/m12"
group_var     <- "sample_id"
pre_reduction <- "pca"
post_reduction<- "postint_harmony"
pre_umap      <- "umap_preintegration"
n_dims        <- 30L
k_primary     <- 15L
k_robust      <- 50L
random_seed   <- 42L
expected_cells<- 19716L
expected_md5  <- "cf63e84313a91de25fe2e41660f78f06"
expected_sha  <- "6085976f788633b68c07c1ddff49e20e6f6f6727579cb7900854febd2ae6e69d"
validation_mode <- FALSE

i <- 1
while (i <= length(args)) {
  k <- args[i]
  if (k == "--input") { input_rds <- args[i+1]; i <- i+2
  } else if (k == "--out-dir") { out_dir <- args[i+1]; i <- i+2
  } else if (k == "--fig-dir") { fig_dir <- args[i+1]; i <- i+2
  } else if (k == "--group-by") { group_var <- args[i+1]; i <- i+2
  } else if (k == "--dims") { n_dims <- as.integer(args[i+1]); i <- i+2
  } else if (k == "--k-primary") { k_primary <- as.integer(args[i+1]); i <- i+2
  } else if (k == "--k-robust") { k_robust <- as.integer(args[i+1]); i <- i+2
  } else if (k == "--random-seed") { random_seed <- as.integer(args[i+1]); i <- i+2
  } else if (k == "--expected-cells") { expected_cells <- as.integer(args[i+1]); i <- i+2
  } else if (k == "--validation-mode") { validation_mode <- TRUE; i <- i+1
  } else stop(sprintf("Unknown argument: %s", k))
}

t_start <- Sys.time()
warnings_collected <- character(0)
record_warning <- function(m) { warnings_collected <<- c(warnings_collected, m); log_warn(m, stage = STAGE) }
fail <- function(m) { log_error(m, stage = STAGE); stop(m, call. = FALSE) }

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

log_info("=============== M12 PRE/POST HARMONY EVALUATION ===============", stage = STAGE)
log_info(sprintf("Input      : %s", input_rds), stage = STAGE)
log_info(sprintf("PRE  space : reduction '%s' dims 1:%d", pre_reduction, n_dims), stage = STAGE)
log_info(sprintf("POST space : reduction '%s' dims 1:%d", post_reduction, n_dims), stage = STAGE)
log_info(sprintf("Grouping   : %s | k = %d (primary), %d (robustness) | seed = %d",
                 group_var, k_primary, k_robust, random_seed), stage = STAGE)

# ---------------------------------------------------------------------------
# GUARD - input integrity
# ---------------------------------------------------------------------------
if (!file.exists(input_rds)) fail(sprintf("M11 object not found: %s", input_rds))
log_info("Verifying M11 output checksums ...", stage = STAGE)
in_md5 <- digest(input_rds, file = TRUE, algo = "md5")
in_sha <- digest(input_rds, file = TRUE, algo = "sha256")
log_info(sprintf("  md5    = %s", in_md5), stage = STAGE)
log_info(sprintf("  sha256 = %s", in_sha), stage = STAGE)
if (!validation_mode) {
  if (!identical(in_md5, expected_md5)) fail(sprintf("M11 md5 mismatch (expected %s).", expected_md5))
  if (!identical(in_sha, expected_sha)) fail(sprintf("M11 sha256 mismatch (expected %s).", expected_sha))
  log_info("PASS - M11 object matches the checksums recorded in M11_REPORT.md.", stage = STAGE)
}

# ---------------------------------------------------------------------------
# Load (read-only) and validate structure
# ---------------------------------------------------------------------------
log_info("Loading M11 Harmony object (read-only) ...", stage = STAGE)
t0 <- Sys.time(); obj <- readRDS(input_rds)
log_info(sprintf("Loaded in %.1f s: %d features x %d cells, default assay '%s'",
                 as.numeric(difftime(Sys.time(), t0, units="secs")),
                 nrow(obj), ncol(obj), DefaultAssay(obj)), stage = STAGE)
log_system_usage(stage = STAGE)

if (!validation_mode && !identical(as.integer(ncol(obj)), expected_cells))
  fail(sprintf("Cell count %d != expected %d.", ncol(obj), expected_cells))
for (r in c(pre_reduction, post_reduction, pre_umap))
  if (!(r %in% Reductions(obj))) fail(sprintf("Required reduction '%s' absent.", r))
if (!(group_var %in% colnames(obj@meta.data))) fail(sprintf("Grouping variable '%s' absent.", group_var))
if (ncol(Embeddings(obj, post_reduction)) < n_dims)
  fail(sprintf("'%s' has %d dims, need %d.", post_reduction, ncol(Embeddings(obj, post_reduction)), n_dims))

digest_pre_umap_before <- digest(Embeddings(obj, pre_umap), algo = "md5")
digest_pca_before      <- digest(Embeddings(obj, pre_reduction), algo = "md5")
digest_harmony_before  <- digest(Embeddings(obj, post_reduction), algo = "md5")
log_info(sprintf("Phase 1 '%s' digest before: %s", pre_umap, digest_pre_umap_before), stage = STAGE)

# ---------------------------------------------------------------------------
# POST-HARMONY UMAP - matched to the Phase 1 baseline UMAP parameters exactly
#
# Phase 1 built 'umap_preintegration' with RunUMAP(reduction='pca', dims=1:30,
# seed.use=42) and every other parameter at the Seurat default (n.neighbors = 30,
# min.dist = 0.3, metric = 'cosine'). The Harmony UMAP uses the identical call
# with only the input reduction changed. Nothing is optimised for appearance.
# ---------------------------------------------------------------------------
post_umap <- "umap_harmony_m12"
log_info(sprintf("Computing post-Harmony UMAP '%s' with parameters matched to Phase 1 ...", post_umap), stage = STAGE)
set.seed(random_seed)
t0 <- Sys.time()
obj <- RunUMAP(obj, reduction = post_reduction, dims = seq_len(n_dims),
               reduction.name = post_umap, reduction.key = "UMAPHARM_",
               seed.use = random_seed, verbose = FALSE)
umap_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
log_info(sprintf("UMAP computed in %.1f s.", umap_seconds), stage = STAGE)

umap_params <- list(
  pre  = list(reduction = pre_umap, input_reduction = pre_reduction, dims = sprintf("1:%d", n_dims),
              seed = 42, n_neighbors = 30, min_dist = 0.3, metric = "cosine",
              source = "Phase 1 M8 (frozen, reused unmodified)"),
  post = list(reduction = post_umap, input_reduction = post_reduction, dims = sprintf("1:%d", n_dims),
              seed = random_seed, n_neighbors = 30, min_dist = 0.3, metric = "cosine",
              source = "M12, Seurat RunUMAP defaults except the input reduction")
)

if (!identical(digest(Embeddings(obj, pre_umap), algo = "md5"), digest_pre_umap_before))
  fail("Phase 1 UMAP was modified. Aborting.")
if (!identical(digest(Embeddings(obj, pre_reduction), algo = "md5"), digest_pca_before))
  fail("Phase 1 PCA was modified. Aborting.")
if (!identical(digest(Embeddings(obj, post_reduction), algo = "md5"), digest_harmony_before))
  fail("M11 Harmony embedding was modified. Aborting.")
log_info("PASS - Phase 1 PCA/UMAP and the M11 Harmony embedding are unchanged.", stage = STAGE)

# ---------------------------------------------------------------------------
# BIOLOGICAL REFERENCE 1 - canonical broad-lineage program scores
#
# Sanity-check programs only. These are NOT cell-type annotations and are never
# written into any saved object. Gene sets are canonical, textbook-level broad
# lineage markers; they are identical for the PRE and POST comparison, so any
# normalisation quirk affects both sides equally and cannot bias the contrast.
#
# Legacy 'orig.anno' is deliberately NOT used: it derives from the original
# integrated analysis, and M12 forbids validating Harmony with legacy
# integration-derived annotation.
# ---------------------------------------------------------------------------
program_sets <- list(
  Panleukocyte  = c("PTPRC"),
  T_NK          = c("CD3D","CD3E","CD2","IL7R","TRAC","NKG7","GNLY","KLRD1"),
  Myeloid       = c("LYZ","CD68","AIF1","CSF1R","ITGAM","C1QA","C1QB","FCGR3A"),
  B_Plasma      = c("MS4A1","CD79A","CD79B","JCHAIN","MZB1"),
  Mast          = c("TPSAB1","TPSB2","CPA3","MS4A2"),
  Endothelial   = c("PECAM1","VWF","CDH5","CLDN5"),
  Fibroblast    = c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA"),
  Mural         = c("ACTA2","RGS5","PDGFRB","MYH11","NOTCH3"),
  SchwannNC     = c("SOX10","S100B","PLP1","MPZ","NGFR","PMP22","ERBB3"),
  Proliferation = c("MKI67","TOP2A","CCNB1","CDK1","PCNA")
)
avail <- rownames(obj)
gene_report <- do.call(rbind, lapply(names(program_sets), function(p) {
  found <- intersect(program_sets[[p]], avail)
  data.frame(program = p, n_requested = length(program_sets[[p]]), n_found = length(found),
             genes_found = paste(found, collapse = ";"),
             genes_missing = paste(setdiff(program_sets[[p]], avail), collapse = ";"),
             stringsAsFactors = FALSE)
}))
write.table(gene_report, file.path(out_dir, "biological_program_gene_sets.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
log_info("Canonical program gene availability:", stage = STAGE)
for (r in seq_len(nrow(gene_report)))
  log_info(sprintf("  %-14s %d/%d found", gene_report$program[r], gene_report$n_found[r], gene_report$n_requested[r]), stage = STAGE)

usable <- program_sets[gene_report$n_found >= 1]
usable <- lapply(usable, function(g) intersect(g, avail))

# AddModuleScore control-feature count. The package default is ctrl = 100 with
# nbin = 24; that requires >= 100 genes in every expression bin. On the real object
# (29,113 features -> ~1,213 genes/bin) the default is used unchanged. The guard only
# lowers ctrl for small synthetic validation inputs, and the value used is logged.
nbin_used <- 24L
ctrl_used <- min(100L, max(5L, as.integer(floor(length(avail) / nbin_used / 2))))
if (ctrl_used != 100L)
  record_warning(sprintf("AddModuleScore ctrl lowered from the default 100 to %d because the object has only %d features (%d bins). Expected only for synthetic validation input.",
                         ctrl_used, length(avail), nbin_used))
log_info(sprintf("AddModuleScore: nbin = %d, ctrl = %d (package default ctrl = 100)", nbin_used, ctrl_used), stage = STAGE)

set.seed(random_seed)
obj <- AddModuleScore(obj, features = usable, name = "m12prog_", seed = random_seed,
                      assay = DefaultAssay(obj), search = FALSE,
                      nbin = nbin_used, ctrl = ctrl_used)
score_cols <- paste0("m12prog_", seq_along(usable))
score_mat <- as.matrix(obj@meta.data[, score_cols, drop = FALSE])
colnames(score_mat) <- names(usable)

# Dominant program label: z-score each program across cells, take the argmax,
# and leave a cell 'Unassigned' when no program stands out (max z < 1).
z <- scale(score_mat)
best_idx <- max.col(z, ties.method = "first")
best_z <- z[cbind(seq_len(nrow(z)), best_idx)]
program_label <- ifelse(best_z >= 1, colnames(z)[best_idx], "Unassigned")
log_info(sprintf("Program sanity-check labels (NOT annotation): %s",
                 paste(sprintf("%s=%d", names(table(program_label)), as.integer(table(program_label))), collapse = ", ")),
         stage = STAGE)

# ---------------------------------------------------------------------------
# Extract everything needed, then free the Seurat object
# ---------------------------------------------------------------------------
meta <- obj@meta.data
cells <- colnames(obj)
grp <- as.character(meta[[group_var]])
preint_cluster <- if ("preint_recommended_cluster" %in% colnames(meta)) as.character(meta$preint_recommended_cluster) else rep(NA_character_, length(cells))
qc <- meta[, intersect(c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo"), colnames(meta)), drop = FALSE]

emb_pre   <- Embeddings(obj, pre_reduction)[, seq_len(n_dims), drop = FALSE]
emb_post  <- Embeddings(obj, post_reduction)[, seq_len(n_dims), drop = FALSE]
umap_pre_xy  <- Embeddings(obj, pre_umap)
umap_post_xy <- Embeddings(obj, post_umap)

write.table(data.frame(cell = rownames(umap_post_xy), UMAP_1 = umap_post_xy[,1], UMAP_2 = umap_post_xy[,2]),
            file.path(out_dir, "umap_harmony_m12_embedding.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(umap_post_xy, file.path(out_dir, "umap_harmony_m12_embedding.rds"))

rm(obj, z, score_mat); invisible(gc(verbose = FALSE))
log_info("Seurat object released from memory; all metrics computed on extracted matrices.", stage = STAGE)
log_system_usage(stage = STAGE)

n_cells <- length(cells)
grp_levels <- sort(unique(grp))
K <- length(grp_levels)
global_p <- as.numeric(table(factor(grp, levels = grp_levels))) / n_cells
names(global_p) <- grp_levels
global_entropy <- -sum(global_p * log(global_p))
log_info(sprintf("Global %s composition: %s", group_var,
                 paste(sprintf("%s=%.4f", grp_levels, global_p), collapse = ", ")), stage = STAGE)
log_info(sprintf("Global composition entropy = %.4f (log K = %.4f). Neighbourhoods matching global",
                 global_entropy, log(K)), stage = STAGE)
log_info("composition would average this value, NOT log K - imbalance makes log K unreachable.", stage = STAGE)

# ---------------------------------------------------------------------------
# Neighbourhood metrics
# ---------------------------------------------------------------------------
knn_indices <- function(mat, k) {
  nn <- RANN::nn2(data = mat, query = mat, k = k + 1L)$nn.idx
  nn[, -1, drop = FALSE]   # drop self
}

neighbourhood_metrics <- function(nn_idx, labels, levels_vec, global_p) {
  k <- ncol(nn_idx)
  lab_codes <- match(labels, levels_vec)
  nb <- matrix(lab_codes[nn_idx], nrow = nrow(nn_idx), ncol = k)
  counts <- vapply(seq_along(levels_vec), function(j) rowSums(nb == j), numeric(nrow(nb)))
  p <- counts / k
  own <- p[cbind(seq_len(nrow(p)), lab_codes)]
  plogp <- ifelse(p > 0, p * log(p), 0)
  ent <- -rowSums(plogp)
  invsimp <- 1 / rowSums(p^2)
  rich <- rowSums(counts > 0)
  list(same_fraction = own, entropy = ent, inv_simpson = invsimp, richness = rich,
       dominance_ratio = own / global_p[lab_codes], counts = counts)
}

metrics_list <- list()
for (kk in c(k_primary, k_robust)) {
  log_info(sprintf("Computing kNN (k = %d) in the PRE and POST spaces ...", kk), stage = STAGE)
  nn_pre  <- knn_indices(emb_pre,  kk)
  nn_post <- knn_indices(emb_post, kk)
  m_pre  <- neighbourhood_metrics(nn_pre,  grp, grp_levels, global_p)
  m_post <- neighbourhood_metrics(nn_post, grp, grp_levels, global_p)

  # kNN retention: how many of a cell's PRE neighbours survive as POST neighbours
  retention <- vapply(seq_len(n_cells), function(i)
    length(intersect(nn_pre[i, ], nn_post[i, ])) / kk, numeric(1))

  metrics_list[[as.character(kk)]] <- list(pre = m_pre, post = m_post, retention = retention,
                                           nn_pre = if (kk == k_primary) nn_pre else NULL,
                                           nn_post = if (kk == k_primary) nn_post else NULL)
  log_info(sprintf("  k=%d  same-%s fraction: PRE %.4f -> POST %.4f", kk, group_var,
                   mean(m_pre$same_fraction), mean(m_post$same_fraction)), stage = STAGE)
  log_info(sprintf("  k=%d  entropy          : PRE %.4f -> POST %.4f (global-composition reference %.4f)",
                   kk, mean(m_pre$entropy), mean(m_post$entropy), global_entropy), stage = STAGE)
  log_info(sprintf("  k=%d  inverse Simpson  : PRE %.4f -> POST %.4f (max %d)",
                   kk, mean(m_pre$inv_simpson), mean(m_post$inv_simpson), K), stage = STAGE)
  log_info(sprintf("  k=%d  kNN retention    : mean %.4f", kk, mean(retention)), stage = STAGE)
  if (kk != k_primary) rm(nn_pre, nn_post)
  invisible(gc(verbose = FALSE))
}

mp <- metrics_list[[as.character(k_primary)]]

# ---------------------------------------------------------------------------
# WITHIN-SAMPLE local structure preservation
#
# Cross-sample neighbour change is exactly what Harmony is supposed to cause.
# WITHIN-sample neighbour change is potential damage to real biology. This is the
# key over-correction diagnostic and is computed separately per sample.
# ---------------------------------------------------------------------------
log_info("Computing WITHIN-sample local structure preservation ...", stage = STAGE)
within_retention <- rep(NA_real_, n_cells)
within_cluster_coherence_pre <- rep(NA_real_, n_cells)
within_cluster_coherence_post <- rep(NA_real_, n_cells)

for (lv in grp_levels) {
  sel <- which(grp == lv)
  kk <- min(k_primary, length(sel) - 1L)
  nnp  <- RANN::nn2(emb_pre[sel, , drop = FALSE],  emb_pre[sel, , drop = FALSE],  k = kk + 1L)$nn.idx[, -1, drop = FALSE]
  nnq  <- RANN::nn2(emb_post[sel, , drop = FALSE], emb_post[sel, , drop = FALSE], k = kk + 1L)$nn.idx[, -1, drop = FALSE]
  within_retention[sel] <- vapply(seq_along(sel), function(i)
    length(intersect(nnp[i, ], nnq[i, ])) / kk, numeric(1))
  cl <- preint_cluster[sel]
  if (!all(is.na(cl))) {
    within_cluster_coherence_pre[sel]  <- rowMeans(matrix(cl[nnp], nrow = length(sel)) == cl)
    within_cluster_coherence_post[sel] <- rowMeans(matrix(cl[nnq], nrow = length(sel)) == cl)
  }
  log_info(sprintf("  %s (n=%d): within-sample kNN retention %.4f | Phase 1 cluster coherence PRE %.4f -> POST %.4f",
                   lv, length(sel), mean(within_retention[sel]),
                   mean(within_cluster_coherence_pre[sel]), mean(within_cluster_coherence_post[sel])), stage = STAGE)
  rm(nnp, nnq)
}
invisible(gc(verbose = FALSE))

# ---------------------------------------------------------------------------
# Silhouette widths - technical separation and biological-program cohesion
# ---------------------------------------------------------------------------
sil_summary <- list()
compute_sils <- function(mat, space_name) {
  log_info(sprintf("Computing distance matrix for the %s space (n = %d) ...", space_name, nrow(mat)), stage = STAGE)
  d <- dist(mat)
  out <- list()
  for (lab_name in c("technical", "biological")) {
    lab <- if (lab_name == "technical") grp else program_label
    f <- factor(lab)
    if (nlevels(f) < 2) next
    s <- cluster::silhouette(as.integer(f), d)
    per_grp <- tapply(s[, 3], f, mean)
    out[[lab_name]] <- list(overall = mean(s[, 3]), per_group = per_grp, per_cell = s[, 3])
    log_info(sprintf("  %s space, %s label ('%s'): mean silhouette width = %.4f",
                     space_name, lab_name,
                     if (lab_name == "technical") group_var else "canonical program",
                     mean(s[, 3])), stage = STAGE)
  }
  rm(d); invisible(gc(verbose = FALSE))
  out
}
sil_pre  <- compute_sils(emb_pre,  "PRE-Harmony")
sil_post <- compute_sils(emb_post, "POST-Harmony")

# ---------------------------------------------------------------------------
# Sample-restricted population check (over-correction risk)
#
# Populations defined from EXPRESSION (canonical program scores), never from
# legacy annotation. A population concentrated in one sample has no counterpart
# to mix with and is the most fragile under batch correction.
# ---------------------------------------------------------------------------
log_info("Checking sample-restricted populations for over-correction ...", stage = STAGE)
restricted_rows <- list()
for (pg in setdiff(unique(program_label), "Unassigned")) {
  sel <- which(program_label == pg)
  if (length(sel) < 30) next
  tab <- table(factor(grp[sel], levels = grp_levels))
  top_share <- max(tab) / sum(tab)
  restricted_rows[[pg]] <- data.frame(
    program = pg, n_cells = length(sel),
    dominant_sample = grp_levels[which.max(tab)],
    dominant_share = as.numeric(top_share),
    sample_restricted = top_share >= 0.60,
    within_sample_knn_retention = mean(within_retention[sel], na.rm = TRUE),
    global_knn_retention = mean(mp$retention[sel]),
    same_sample_fraction_pre = mean(mp$pre$same_fraction[sel]),
    same_sample_fraction_post = mean(mp$post$same_fraction[sel]),
    program_silhouette_pre = if (!is.null(sil_pre$biological)) as.numeric(sil_pre$biological$per_group[pg]) else NA_real_,
    program_silhouette_post = if (!is.null(sil_post$biological)) as.numeric(sil_post$biological$per_group[pg]) else NA_real_,
    stringsAsFactors = FALSE)
}
restricted_df <- do.call(rbind, restricted_rows)
if (!is.null(restricted_df)) {
  write.table(restricted_df, file.path(out_dir, "sample_restricted_population_check.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  for (r in seq_len(nrow(restricted_df)))
    log_info(sprintf("  %-14s n=%5d  dominant %s (%.1f%%)  within-sample retention %.3f  program silhouette %.3f -> %.3f",
                     restricted_df$program[r], restricted_df$n_cells[r], restricted_df$dominant_sample[r],
                     100*restricted_df$dominant_share[r], restricted_df$within_sample_knn_retention[r],
                     restricted_df$program_silhouette_pre[r], restricted_df$program_silhouette_post[r]), stage = STAGE)
}

# ---------------------------------------------------------------------------
# Machine-readable outputs
# ---------------------------------------------------------------------------
log_info("Writing machine-readable summaries ...", stage = STAGE)

summ <- function(x) c(mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE),
                      min = min(x, na.rm=TRUE), q25 = unname(quantile(x, .25, na.rm=TRUE)),
                      median = median(x, na.rm=TRUE), q75 = unname(quantile(x, .75, na.rm=TRUE)),
                      max = max(x, na.rm=TRUE))

rows <- list()
add_row <- function(metric, k, space, x, desirable, definition) {
  s <- summ(x)
  rows[[length(rows)+1]] <<- data.frame(metric = metric, k = k, space = space,
    mean = s["mean"], sd = s["sd"], min = s["min"], q25 = s["q25"], median = s["median"],
    q75 = s["q75"], max = s["max"], desirable_direction = desirable, definition = definition,
    stringsAsFactors = FALSE)
}
for (kk in c(k_primary, k_robust)) {
  m <- metrics_list[[as.character(kk)]]
  add_row("same_sample_neighbor_fraction", kk, "pre_harmony_pca", m$pre$same_fraction, "lower",
          "Fraction of a cell's k nearest neighbours drawn from its own sample_id.")
  add_row("same_sample_neighbor_fraction", kk, "post_harmony", m$post$same_fraction, "lower",
          "Fraction of a cell's k nearest neighbours drawn from its own sample_id.")
  add_row("neighborhood_entropy", kk, "pre_harmony_pca", m$pre$entropy, "higher (toward global-composition entropy)",
          "Shannon entropy of the sample_id composition of a cell's k nearest neighbours (natural log).")
  add_row("neighborhood_entropy", kk, "post_harmony", m$post$entropy, "higher (toward global-composition entropy)",
          "Shannon entropy of the sample_id composition of a cell's k nearest neighbours (natural log).")
  add_row("neighborhood_inverse_simpson", kk, "pre_harmony_pca", m$pre$inv_simpson, "higher (max = 4)",
          "iLISI-equivalent: 1 / sum(p^2) over sample_id proportions among k exact nearest neighbours.")
  add_row("neighborhood_inverse_simpson", kk, "post_harmony", m$post$inv_simpson, "higher (max = 4)",
          "iLISI-equivalent: 1 / sum(p^2) over sample_id proportions among k exact nearest neighbours.")
  add_row("neighborhood_sample_richness", kk, "pre_harmony_pca", m$pre$richness, "higher (max = 4)",
          "Number of distinct sample_id values represented among a cell's k nearest neighbours.")
  add_row("neighborhood_sample_richness", kk, "post_harmony", m$post$richness, "higher (max = 4)",
          "Number of distinct sample_id values represented among a cell's k nearest neighbours.")
  add_row("sample_dominance_ratio", kk, "pre_harmony_pca", m$pre$dominance_ratio, "toward 1.0",
          "Observed same-sample neighbour fraction divided by that sample's global abundance. 1.0 = neighbourhoods match global composition.")
  add_row("sample_dominance_ratio", kk, "post_harmony", m$post$dominance_ratio, "toward 1.0",
          "Observed same-sample neighbour fraction divided by that sample's global abundance. 1.0 = neighbourhoods match global composition.")
  add_row("knn_retention_global", kk, "pre_to_post", m$retention, "context only",
          "Fraction of a cell's PRE-Harmony k nearest neighbours still among its POST-Harmony k nearest neighbours. A drop is expected and is not by itself damage.")
}
add_row("knn_retention_within_sample", k_primary, "pre_to_post", within_retention, "higher",
        "Same as knn_retention_global but neighbours restricted to the cell's own sample. Measures damage to within-sample biology; cross-sample rearrangement is excluded by construction.")
add_row("phase1_cluster_coherence", k_primary, "pre_harmony_pca", within_cluster_coherence_pre, "preserved",
        "Fraction of a cell's within-sample k nearest neighbours sharing its independent Phase 1 cluster label (preint_recommended_cluster).")
add_row("phase1_cluster_coherence", k_primary, "post_harmony", within_cluster_coherence_post, "preserved",
        "Fraction of a cell's within-sample k nearest neighbours sharing its independent Phase 1 cluster label (preint_recommended_cluster).")
mix_df <- do.call(rbind, rows); rownames(mix_df) <- NULL
write.table(mix_df, file.path(out_dir, "pre_post_mixing_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Per-sample breakdown at k_primary
per_sample <- do.call(rbind, lapply(grp_levels, function(lv) {
  sel <- grp == lv
  data.frame(sample_id = lv, n_cells = sum(sel), global_share = global_p[[lv]],
    same_sample_fraction_pre = mean(mp$pre$same_fraction[sel]),
    same_sample_fraction_post = mean(mp$post$same_fraction[sel]),
    entropy_pre = mean(mp$pre$entropy[sel]), entropy_post = mean(mp$post$entropy[sel]),
    inv_simpson_pre = mean(mp$pre$inv_simpson[sel]), inv_simpson_post = mean(mp$post$inv_simpson[sel]),
    dominance_ratio_pre = mean(mp$pre$dominance_ratio[sel]), dominance_ratio_post = mean(mp$post$dominance_ratio[sel]),
    knn_retention_global = mean(mp$retention[sel]),
    knn_retention_within_sample = mean(within_retention[sel], na.rm = TRUE),
    phase1_cluster_coherence_pre = mean(within_cluster_coherence_pre[sel], na.rm = TRUE),
    phase1_cluster_coherence_post = mean(within_cluster_coherence_post[sel], na.rm = TRUE),
    stringsAsFactors = FALSE)
}))
write.table(per_sample, file.path(out_dir, "neighborhood_mixing_metrics_by_sample.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Silhouette summary
sil_rows <- list()
for (lab_name in c("technical","biological")) {
  if (is.null(sil_pre[[lab_name]]) || is.null(sil_post[[lab_name]])) next
  sil_rows[[length(sil_rows)+1]] <- data.frame(
    label_type = lab_name,
    label = if (lab_name == "technical") group_var else "canonical_program_sanitycheck",
    group = "ALL",
    silhouette_pre = sil_pre[[lab_name]]$overall,
    silhouette_post = sil_post[[lab_name]]$overall,
    change = sil_post[[lab_name]]$overall - sil_pre[[lab_name]]$overall,
    desirable_direction = if (lab_name == "technical") "decrease toward 0" else "preserved (large drop = lost biology)",
    stringsAsFactors = FALSE)
  gs <- union(names(sil_pre[[lab_name]]$per_group), names(sil_post[[lab_name]]$per_group))
  for (g in gs) sil_rows[[length(sil_rows)+1]] <- data.frame(
    label_type = lab_name,
    label = if (lab_name == "technical") group_var else "canonical_program_sanitycheck",
    group = g,
    silhouette_pre = as.numeric(sil_pre[[lab_name]]$per_group[g]),
    silhouette_post = as.numeric(sil_post[[lab_name]]$per_group[g]),
    change = as.numeric(sil_post[[lab_name]]$per_group[g]) - as.numeric(sil_pre[[lab_name]]$per_group[g]),
    desirable_direction = if (lab_name == "technical") "decrease toward 0" else "preserved",
    stringsAsFactors = FALSE)
}
sil_df <- do.call(rbind, sil_rows)
write.table(sil_df, file.path(out_dir, "technical_silhouette_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Biological preservation summary
bio_rows <- data.frame(
  measure = c("within_sample_knn_retention_mean",
              "phase1_cluster_coherence_pre_mean",
              "phase1_cluster_coherence_post_mean",
              "phase1_cluster_coherence_relative_change",
              "program_silhouette_pre", "program_silhouette_post", "program_silhouette_change",
              "global_knn_retention_mean"),
  value = c(mean(within_retention, na.rm = TRUE),
            mean(within_cluster_coherence_pre, na.rm = TRUE),
            mean(within_cluster_coherence_post, na.rm = TRUE),
            (mean(within_cluster_coherence_post, na.rm=TRUE) - mean(within_cluster_coherence_pre, na.rm=TRUE)) /
              mean(within_cluster_coherence_pre, na.rm=TRUE),
            if (!is.null(sil_pre$biological)) sil_pre$biological$overall else NA_real_,
            if (!is.null(sil_post$biological)) sil_post$biological$overall else NA_real_,
            if (!is.null(sil_pre$biological)) sil_post$biological$overall - sil_pre$biological$overall else NA_real_,
            mean(mp$retention)),
  interpretation = c(
    "Fraction of within-sample nearest neighbours preserved by Harmony. High = within-sample biology intact.",
    "Agreement between within-sample neighbourhoods and independent Phase 1 clusters, PRE-Harmony.",
    "Same, POST-Harmony.",
    "Relative change in Phase 1 cluster coherence. Near zero = independent pre-integration structure retained.",
    "Cohesion of canonical broad-lineage programs, PRE-Harmony.",
    "Cohesion of canonical broad-lineage programs, POST-Harmony.",
    "Change in program cohesion. A large negative value would indicate lineage structure was collapsed.",
    "Fraction of global nearest neighbours preserved. Expected to fall; falling is not by itself damage."),
  stringsAsFactors = FALSE)
write.table(bio_rows, file.path(out_dir, "biological_preservation_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Program x sample composition
prog_comp <- as.data.frame(table(program = program_label, sample_id = grp), stringsAsFactors = FALSE)
write.table(prog_comp, file.path(out_dir, "program_by_sample_composition.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Evaluation parameters
eval_params <- data.frame(
  parameter = c("input_object","input_md5","input_sha256","cells_evaluated","subsampling",
                "pre_space","post_space","dims","k_primary","k_robust","random_seed",
                "pre_umap_reduction","post_umap_reduction","umap_params_matched",
                "knn_method","silhouette_method","distance_metric",
                "biological_reference_1","biological_reference_2","legacy_orig_anno_used",
                "module_score_nbin","module_score_ctrl",
                "lisi_package","kBET_package","harmony_retuned"),
  value = c(input_rds, in_md5, in_sha, n_cells, "NONE - all cells evaluated",
            sprintf("%s dims 1:%d", pre_reduction, n_dims), sprintf("%s dims 1:%d", post_reduction, n_dims),
            n_dims, k_primary, k_robust, random_seed,
            pre_umap, post_umap, "yes - identical dims, seed, n.neighbors, min.dist, metric",
            "RANN::nn2 exact k-nearest neighbours", "cluster::silhouette on full pairwise distances",
            "euclidean", "independent Phase 1 per-sample clusters (preint_recommended_cluster)",
            "canonical broad-lineage program scores (Seurat AddModuleScore)",
            "NO - excluded as legacy integration-derived annotation",
            nbin_used, ctrl_used,
            "not installed - inverse Simpson computed natively on exact kNN",
            "not installed - replaced by dominance-ratio vs global composition",
            "NO - M11 default result evaluated as-is"),
  stringsAsFactors = FALSE)
write.table(eval_params, file.path(out_dir, "harmony_evaluation_parameters.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
log_info("Generating matched pre/post figures ...", stage = STAGE)

pal_sample <- setNames(RColorBrewer::brewer.pal(max(3, K), "Set1")[seq_len(K)], grp_levels)
base_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.key.size = unit(0.4, "cm"))

df_pre  <- data.frame(x = umap_pre_xy[,1],  y = umap_pre_xy[,2],  sample_id = grp,
                      program = program_label, preint = preint_cluster, stringsAsFactors = FALSE)
df_post <- data.frame(x = umap_post_xy[,1], y = umap_post_xy[,2], sample_id = grp,
                      program = program_label, preint = preint_cluster, stringsAsFactors = FALSE)
for (q in colnames(qc)) { df_pre[[q]] <- qc[[q]]; df_post[[q]] <- qc[[q]] }

# deterministic plotting order so neither panel is favoured by overplotting
set.seed(random_seed); ord <- sample.int(n_cells)
df_pre <- df_pre[ord, ]; df_post <- df_post[ord, ]

figure_rows <- list()
save_fig <- function(plot, name, w, h, method, params) {
  pdf_p <- file.path(fig_dir, paste0(name, ".pdf"))
  png_p <- file.path(fig_dir, paste0(name, ".png"))
  ggsave(pdf_p, plot, width = w, height = h, units = "in", device = "pdf")
  ggsave(png_p, plot, width = w, height = h, units = "in", dpi = 200, device = "png")
  for (p in c(pdf_p, png_p))
    figure_rows[[length(figure_rows)+1]] <<- data.frame(
      figure_path = p, phase = "phase2", milestone = "M12", dataset = "combined",
      processing_stage = "pre_post_harmony_evaluation", analysis_method = method,
      parameters = params, input_object_checksum = in_md5,
      generating_script = "scripts/R/phase2/evaluate_harmony.R",
      stringsAsFactors = FALSE)
  log_info(sprintf("  wrote %s(.pdf/.png)", name), stage = STAGE)
}

umap_panel <- function(df, title, colour_by, discrete = TRUE, pal = NULL, legend = TRUE) {
  p <- ggplot(df, aes(x = x, y = y, colour = .data[[colour_by]])) +
    geom_point(size = 0.18, alpha = 0.75, stroke = 0) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2", colour = colour_by) + base_theme
  if (discrete) {
    if (!is.null(pal)) p <- p + scale_colour_manual(values = pal)
    p <- p + guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
  } else p <- p + scale_colour_viridis_c(option = "viridis")
  if (!legend) p <- p + theme(legend.position = "none")
  p
}

# 01 - sample_id (also the dataset, patient and Harmony grouping variable)
p01 <- umap_panel(df_pre,  "PRE-Harmony  (pca 1:30 -> umap_preintegration)", "sample_id", TRUE, pal_sample) |
       umap_panel(df_post, "POST-Harmony (postint_harmony 1:30 -> umap_harmony_m12)", "sample_id", TRUE, pal_sample)
p01 <- p01 + plot_layout(guides = "collect") +
  plot_annotation(title = "Matched pre/post Harmony UMAP coloured by sample_id",
                  subtitle = "sample_id is simultaneously the dataset, the patient and the Harmony grouping variable (all four are the same 4-level field)")
save_fig(p01, "01_pre_post_umap_by_sample_id", 13, 5.5, "UMAP_scatter_paired",
         sprintf("dims=1:%d;seed=%d;n.neighbors=30;min.dist=0.3;metric=cosine;point=0.18", n_dims, random_seed))

# 02 - faceted composition/density
df_facet <- rbind(cbind(df_pre, space = "PRE-Harmony"), cbind(df_post, space = "POST-Harmony"))
df_facet$space <- factor(df_facet$space, levels = c("PRE-Harmony","POST-Harmony"))
p02 <- ggplot(df_facet, aes(x = x, y = y)) +
  geom_point(data = transform(df_facet, sample_id = NULL), colour = "grey88", size = 0.12, stroke = 0) +
  geom_point(aes(colour = sample_id), size = 0.16, alpha = 0.8, stroke = 0) +
  scale_colour_manual(values = pal_sample) +
  facet_grid(space ~ sample_id, scales = "free") +
  labs(title = "Per-sample occupancy of the shared embedding, pre vs post Harmony",
       subtitle = "Grey = all cells; colour = the faceted sample. Shows whether any sample dominates a region.",
       x = "UMAP 1", y = "UMAP 2") +
  base_theme + theme(legend.position = "none")
save_fig(p02, "02_pre_post_umap_faceted_by_sample_id", 14, 7, "UMAP_facet_density",
         sprintf("dims=1:%d;seed=%d;facet=sample_id x space", n_dims, random_seed))

# 03 - QC covariates
qc_plots <- list()
for (q in colnames(qc)) {
  qc_plots[[paste0(q, "_pre")]]  <- umap_panel(df_pre,  paste0("PRE  - ", q),  q, FALSE)
  qc_plots[[paste0(q, "_post")]] <- umap_panel(df_post, paste0("POST - ", q), q, FALSE)
}
p03 <- wrap_plots(qc_plots, ncol = 2) +
  plot_annotation(title = "Matched pre/post Harmony UMAP coloured by technical QC covariates")
save_fig(p03, "03_pre_post_umap_by_qc_covariates", 11, 4.2 * length(colnames(qc)), "UMAP_scatter_continuous",
         sprintf("dims=1:%d;seed=%d;metrics=%s", n_dims, random_seed, paste(colnames(qc), collapse=",")))

# 04 - independent Phase 1 clusters (provenance/biological reference 1)
p04 <- umap_panel(df_pre, "PRE-Harmony", "preint", TRUE, NULL, legend = FALSE) |
       umap_panel(df_post, "POST-Harmony", "preint", TRUE, NULL, legend = FALSE)
p04 <- p04 + plot_annotation(
  title = "Independent Phase 1 per-sample clusters (preint_recommended_cluster) in both spaces",
  subtitle = "54 sample-specific clusters derived independently before integration; legend omitted for readability. Used as a non-legacy biological reference.")
save_fig(p04, "04_pre_post_umap_by_phase1_independent_clusters", 13, 5.5, "UMAP_scatter_paired",
         sprintf("dims=1:%d;seed=%d;label=preint_recommended_cluster;legend=suppressed", n_dims, random_seed))

# 05-08 - metric distributions
metric_density <- function(pre_v, post_v, title, xlab, vline = NULL) {
  d <- rbind(data.frame(value = pre_v,  space = "PRE-Harmony"),
             data.frame(value = post_v, space = "POST-Harmony"))
  d$space <- factor(d$space, levels = c("PRE-Harmony","POST-Harmony"))
  p <- ggplot(d, aes(x = value, fill = space, colour = space)) +
    geom_density(alpha = 0.35, adjust = 1.2) +
    scale_fill_manual(values = c("PRE-Harmony" = "#3B6EA5", "POST-Harmony" = "#C1452B")) +
    scale_colour_manual(values = c("PRE-Harmony" = "#3B6EA5", "POST-Harmony" = "#C1452B")) +
    labs(title = title, x = xlab, y = "density") + base_theme
  if (!is.null(vline)) p <- p + geom_vline(xintercept = vline, linetype = "dashed", colour = "grey30")
  p
}
p05 <- metric_density(mp$pre$entropy, mp$post$entropy,
        sprintf("Neighbourhood sample_id entropy (k = %d)", k_primary),
        "Shannon entropy (nats)", vline = global_entropy) +
  labs(caption = sprintf("Dashed line = global-composition entropy %.4f, the value expected if neighbourhoods matched overall sample abundances. log(K) = %.4f is unreachable given 3.3x imbalance.",
                         global_entropy, log(K)))
save_fig(p05, "05_pre_post_neighborhood_entropy", 8, 5, "density_metric",
         sprintf("k=%d;metric=shannon_entropy;reference=global_composition_entropy", k_primary))

p06 <- metric_density(mp$pre$same_fraction, mp$post$same_fraction,
        sprintf("Same-sample nearest-neighbour fraction (k = %d)", k_primary),
        "fraction of k neighbours from the same sample_id")
save_fig(p06, "06_pre_post_same_sample_neighbor_fraction", 8, 5, "density_metric",
         sprintf("k=%d;metric=same_sample_fraction", k_primary))

p07 <- metric_density(mp$pre$inv_simpson, mp$post$inv_simpson,
        sprintf("Neighbourhood inverse Simpson index / iLISI-equivalent (k = %d)", k_primary),
        "inverse Simpson (1 = single sample, 4 = all four equally)")
save_fig(p07, "07_pre_post_inverse_simpson", 8, 5, "density_metric",
         sprintf("k=%d;metric=inverse_simpson;max=%d", k_primary, K))

sil_plot_df <- sil_df[sil_df$group != "ALL", ]
sil_long <- rbind(
  data.frame(label_type = sil_plot_df$label_type, group = sil_plot_df$group,
             space = "PRE-Harmony", value = sil_plot_df$silhouette_pre),
  data.frame(label_type = sil_plot_df$label_type, group = sil_plot_df$group,
             space = "POST-Harmony", value = sil_plot_df$silhouette_post))
sil_long$space <- factor(sil_long$space, levels = c("PRE-Harmony","POST-Harmony"))
sil_long$label_type <- factor(sil_long$label_type, levels = c("technical","biological"),
  labels = c("TECHNICAL: sample_id (want ~0)", "BIOLOGICAL: canonical program (want preserved)"))
p08 <- ggplot(sil_long, aes(x = reorder(group, value), y = value, fill = space)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_hline(yintercept = 0, colour = "grey30") +
  coord_flip() + facet_wrap(~ label_type, scales = "free_y") +
  scale_fill_manual(values = c("PRE-Harmony" = "#3B6EA5", "POST-Harmony" = "#C1452B")) +
  labs(title = "Silhouette width by group, pre vs post Harmony",
       subtitle = "Technical separation should fall toward 0; biological program cohesion should be retained",
       x = NULL, y = "mean silhouette width") + base_theme
save_fig(p08, "08_pre_post_technical_silhouette", 11, 6, "silhouette_barplot",
         sprintf("dims=1:%d;metric=euclidean;labels=sample_id,canonical_program", n_dims))

# 09 - kNN retention, global vs within-sample
ret_df <- rbind(data.frame(value = mp$retention, type = "Global kNN retention"),
                data.frame(value = within_retention, type = "Within-sample kNN retention"))
p09 <- ggplot(ret_df, aes(x = value, fill = type, colour = type)) +
  geom_density(alpha = 0.35, adjust = 1.2) +
  scale_fill_manual(values = c("Global kNN retention" = "#7A7A7A", "Within-sample kNN retention" = "#2E8B57")) +
  scale_colour_manual(values = c("Global kNN retention" = "#7A7A7A", "Within-sample kNN retention" = "#2E8B57")) +
  labs(title = sprintf("Nearest-neighbour retention from the pre- to the post-Harmony space (k = %d)", k_primary),
       subtitle = "Global retention is expected to fall (that is the correction). A fall in WITHIN-sample retention would indicate damage to real biology.",
       x = "fraction of pre-Harmony neighbours retained", y = "density") + base_theme
save_fig(p09, "09_knn_retention_global_vs_within_sample", 9, 5, "density_metric",
         sprintf("k=%d;metric=knn_retention", k_primary))

# 10 - biological preservation summary
bio_plot <- data.frame(
  measure = factor(c("Phase 1 cluster coherence","Phase 1 cluster coherence",
                     "Canonical program silhouette","Canonical program silhouette"),
                   levels = c("Phase 1 cluster coherence","Canonical program silhouette")),
  space = factor(c("PRE-Harmony","POST-Harmony","PRE-Harmony","POST-Harmony"),
                 levels = c("PRE-Harmony","POST-Harmony")),
  value = c(mean(within_cluster_coherence_pre, na.rm=TRUE), mean(within_cluster_coherence_post, na.rm=TRUE),
            if (!is.null(sil_pre$biological)) sil_pre$biological$overall else NA_real_,
            if (!is.null(sil_post$biological)) sil_post$biological$overall else NA_real_))
p10 <- ggplot(bio_plot, aes(x = space, y = value, fill = space)) +
  geom_col(width = 0.6) + geom_text(aes(label = sprintf("%.4f", value)), vjust = -0.4, size = 3.4) +
  facet_wrap(~ measure, scales = "free_y") +
  scale_fill_manual(values = c("PRE-Harmony" = "#3B6EA5", "POST-Harmony" = "#C1452B")) +
  labs(title = "Biological preservation summary",
       subtitle = "Left: agreement of within-sample neighbourhoods with independent Phase 1 clusters. Right: cohesion of canonical broad-lineage programs.",
       x = NULL, y = NULL) + base_theme + theme(legend.position = "none")
save_fig(p10, "10_biological_preservation_summary", 9, 5, "barplot_summary",
         sprintf("k=%d;references=preint_recommended_cluster,canonical_programs", k_primary))

# 11 - program sanity-check labels
p11 <- umap_panel(df_pre, "PRE-Harmony", "program", TRUE) |
       umap_panel(df_post, "POST-Harmony", "program", TRUE)
p11 <- p11 + plot_layout(guides = "collect") + plot_annotation(
  title = "Canonical broad-lineage program sanity-check labels in both spaces",
  subtitle = "Dominant canonical program per cell (z >= 1). A preservation sanity check only - NOT a cell-type annotation, which is M15.")
save_fig(p11, "11_pre_post_umap_by_canonical_program", 13, 5.5, "UMAP_scatter_paired",
         sprintf("dims=1:%d;seed=%d;labels=argmax_z_program(threshold=1)", n_dims, random_seed))

# 12 - local dominance vs global abundance
dom_df <- rbind(data.frame(value = mp$pre$dominance_ratio,  space = "PRE-Harmony",  sample_id = grp),
                data.frame(value = mp$post$dominance_ratio, space = "POST-Harmony", sample_id = grp))
dom_df$space <- factor(dom_df$space, levels = c("PRE-Harmony","POST-Harmony"))
p12 <- ggplot(dom_df, aes(x = sample_id, y = value, fill = space)) +
  geom_boxplot(outlier.size = 0.15, outlier.alpha = 0.2, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(values = c("PRE-Harmony" = "#3B6EA5", "POST-Harmony" = "#C1452B")) +
  labs(title = "Local sample dominance relative to global abundance",
       subtitle = "Observed same-sample neighbour fraction / that sample's global share. Dashed line at 1.0 = neighbourhoods match global composition; >1 = locally over-represented.",
       x = NULL, y = "dominance ratio") + base_theme
save_fig(p12, "12_local_sample_dominance_vs_global", 9, 5.5, "boxplot_metric",
         sprintf("k=%d;metric=dominance_ratio;reference=global_composition", k_primary))

fig_df <- do.call(rbind, figure_rows)
write.table(fig_df, file.path(out_dir, "figure_index_m12.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------
log_info("Recording provenance ...", stage = STAGE)
prov_rec <- record_provenance(
  step_name = STAGE,
  inputs = list(m11_harmony_rds = input_rds),
  outputs = list(
    pre_post_mixing_summary = file.path(out_dir, "pre_post_mixing_summary.tsv"),
    neighborhood_mixing_metrics_by_sample = file.path(out_dir, "neighborhood_mixing_metrics_by_sample.tsv"),
    technical_silhouette_summary = file.path(out_dir, "technical_silhouette_summary.tsv"),
    biological_preservation_summary = file.path(out_dir, "biological_preservation_summary.tsv"),
    sample_restricted_population_check = file.path(out_dir, "sample_restricted_population_check.tsv"),
    biological_program_gene_sets = file.path(out_dir, "biological_program_gene_sets.tsv"),
    program_by_sample_composition = file.path(out_dir, "program_by_sample_composition.tsv"),
    harmony_evaluation_parameters = file.path(out_dir, "harmony_evaluation_parameters.tsv"),
    umap_harmony_m12_embedding = file.path(out_dir, "umap_harmony_m12_embedding.rds"),
    figure_index_m12 = file.path(out_dir, "figure_index_m12.tsv")),
  parameters = list(
    milestone = "M12", cells_evaluated = n_cells, subsampling = "none",
    pre_space = sprintf("%s 1:%d", pre_reduction, n_dims),
    post_space = sprintf("%s 1:%d", post_reduction, n_dims),
    k_primary = k_primary, k_robust = k_robust, random_seed = random_seed,
    umap_parameters = umap_params, harmony_retuned = FALSE,
    input_md5 = in_md5, input_sha256 = in_sha,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", "NOT_IN_SLURM")),
  dataset = "combined")
save_provenance_json(prov_rec, file.path(out_dir, "prov_m12_evaluation.json"))

# Headline JSON
headline <- list(
  milestone = "M12", phase = "phase2",
  cells_evaluated = n_cells, subsampling = "none - all cells",
  k_primary = k_primary, k_robust = k_robust,
  global_composition = as.list(global_p), global_composition_entropy = global_entropy, log_K = log(K),
  technical_mixing = list(
    same_sample_fraction_pre = mean(mp$pre$same_fraction),
    same_sample_fraction_post = mean(mp$post$same_fraction),
    entropy_pre = mean(mp$pre$entropy), entropy_post = mean(mp$post$entropy),
    inv_simpson_pre = mean(mp$pre$inv_simpson), inv_simpson_post = mean(mp$post$inv_simpson),
    dominance_ratio_pre = mean(mp$pre$dominance_ratio), dominance_ratio_post = mean(mp$post$dominance_ratio),
    silhouette_sample_pre = if (!is.null(sil_pre$technical)) sil_pre$technical$overall else NA_real_,
    silhouette_sample_post = if (!is.null(sil_post$technical)) sil_post$technical$overall else NA_real_),
  biological_preservation = list(
    within_sample_knn_retention = mean(within_retention, na.rm = TRUE),
    global_knn_retention = mean(mp$retention),
    phase1_cluster_coherence_pre = mean(within_cluster_coherence_pre, na.rm = TRUE),
    phase1_cluster_coherence_post = mean(within_cluster_coherence_post, na.rm = TRUE),
    program_silhouette_pre = if (!is.null(sil_pre$biological)) sil_pre$biological$overall else NA_real_,
    program_silhouette_post = if (!is.null(sil_post$biological)) sil_post$biological$overall else NA_real_),
  phase1_reference_values = list(same_dataset_fraction = 0.9613, entropy = 0.0966,
                                 note = "Phase 1 M8 computed these on the same pca 1:30 space with k=15."),
  warnings = warnings_collected,
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "NOT_IN_SLURM"),
  git_commit = tryCatch(trimws(system("git rev-parse HEAD", intern = TRUE)), error = function(e) NA_character_),
  elapsed_seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
  session_info = capture.output(sessionInfo()))
write_json(headline, file.path(out_dir, "m12_headline_metrics.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)

log_info("---------------------------------------------------------------", stage = STAGE)
log_info(sprintf("M12 COMPLETE in %.1f s (%.2f min). Warnings: %d",
                 as.numeric(difftime(Sys.time(), t_start, units="secs")),
                 as.numeric(difftime(Sys.time(), t_start, units="mins")), length(warnings_collected)), stage = STAGE)
log_info("Harmony was NOT re-run or re-tuned. No clustering, markers or annotation were produced.", stage = STAGE)
log_system_usage(stage = STAGE)
