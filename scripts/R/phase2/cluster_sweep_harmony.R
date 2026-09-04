# scripts/R/phase2/cluster_sweep_harmony.R
#
# Phase 2 / Milestone M13 - Post-Harmony neighbours, UMAP and clustering resolution sweep.
#
# Consumes the M11 Harmony object (M12-accepted), builds a neighbour graph and UMAP on the
# Harmony embedding, sweeps Louvain resolutions 0.1-1.0, preserves every solution, scores
# them on stability / fragmentation / compactness / compartment coverage, and recommends a
# primary and an alternative resolution.
#
# Guarantees: Phase 1 and M11 reductions, graphs and metadata are preserved byte-identical.
# In particular the legacy `seurat_clusters` column is backed up and restored, because
# Seurat's FindClusters() overwrites it as a side effect.
#
# NOT done here: markers, annotation, DE.

options(stringsAsFactors = FALSE); options(future.globals.maxSize = +Inf)
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(cluster); library(aricode)
  library(ggplot2); library(patchwork); library(RColorBrewer); library(viridis)
  library(jsonlite); library(digest)
})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase2_m13_clustering"

args <- commandArgs(trailingOnly = TRUE)
input_rds <- "results/phase2/harmony/phase2_harmony_integrated.rds"
out_rds   <- "results/phase2/clustering/phase2_harmony_clustered.rds"
out_dir   <- "results/phase2/clustering"
fig_dir   <- "results/phase2/figures/M13"
red       <- "postint_harmony"
n_dims    <- 30L
k_param   <- 20L
seed      <- 42L
expected_cells <- 19716L
expected_md5   <- "cf63e84313a91de25fe2e41660f78f06"
resolutions <- seq(0.1, 1.0, by = 0.1)
validation_mode <- FALSE
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--input"){input_rds<-args[i+1];i<-i+2} else if (a=="--out-rds"){out_rds<-args[i+1];i<-i+2}
  else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2} else if (a=="--fig-dir"){fig_dir<-args[i+1];i<-i+2}
  else if (a=="--dims"){n_dims<-as.integer(args[i+1]);i<-i+2} else if (a=="--k-param"){k_param<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--random-seed"){seed<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--expected-cells"){expected_cells<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--validation-mode"){validation_mode<-TRUE;i<-i+1}
  else stop(sprintf("Unknown argument: %s", a)) }

t_start <- Sys.time(); warns <- character(0)
wrn <- function(m){ warns <<- c(warns,m); log_warn(m, stage=STAGE) }
fail <- function(m){ log_error(m, stage=STAGE); stop(m, call.=FALSE) }
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE); dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dims <- seq_len(n_dims)

log_info("================= M13 POST-HARMONY CLUSTERING =================", stage=STAGE)
log_info(sprintf("Input %s | reduction '%s' dims 1:%d | k.param %d | seed %d", input_rds, red, n_dims, k_param, seed), stage=STAGE)

if (file.exists(out_rds) && !validation_mode) fail(sprintf("Output already exists: %s. Refusing to overwrite.", out_rds))
in_md5 <- digest(input_rds, file=TRUE, algo="md5")
log_info(sprintf("Input md5 = %s", in_md5), stage=STAGE)
if (!validation_mode && !identical(in_md5, expected_md5)) fail("M11 object md5 mismatch.")

log_info("Loading M11 Harmony object ...", stage=STAGE)
obj <- readRDS(input_rds)
log_info(sprintf("Loaded: %d features x %d cells, assay '%s'", nrow(obj), ncol(obj), DefaultAssay(obj)), stage=STAGE)
if (!validation_mode && ncol(obj) != expected_cells) fail("Cell count mismatch.")
if (!(red %in% Reductions(obj))) fail(sprintf("Reduction '%s' absent.", red))
if (ncol(Embeddings(obj, red)) < n_dims) fail("Insufficient Harmony dimensions.")

# ---- preservation snapshot -------------------------------------------------
pre <- list(meta = digest(obj@meta.data, algo="md5"),
            meta_cols = colnames(obj@meta.data),
            meta_percol = lapply(obj@meta.data, function(v) list(d = digest(v, algo="md5"), cls = class(v)[1])),
            pca = digest(Embeddings(obj,"pca"), algo="md5"),
            umap_pre = digest(Embeddings(obj,"umap_preintegration"), algo="md5"),
            harmony = digest(Embeddings(obj, red), algo="md5"),
            graphs = sort(Graphs(obj)), reductions = sort(Reductions(obj)),
            assays = sort(Assays(obj)), default_assay = DefaultAssay(obj),
            cells = digest(colnames(obj), algo="md5"))
# FindClusters() overwrites `seurat_clusters` and sets Idents(). `seurat_clusters` is a
# LEGACY Phase 1 column, so snapshot it and restore it verbatim afterwards.
legacy_seurat_clusters <- if ("seurat_clusters" %in% colnames(obj@meta.data)) obj@meta.data$seurat_clusters else NULL
legacy_idents <- Idents(obj)
orig_meta_backup <- obj@meta.data[, pre$meta_cols, drop = FALSE]   # verbatim copy for restoration

# ---- neighbours ------------------------------------------------------------
log_info("FindNeighbors on the Harmony embedding ...", stage=STAGE)
set.seed(seed)
obj <- FindNeighbors(obj, reduction=red, dims=dims, k.param=k_param,
                     graph.name=c("postint_harmony_nn","postint_harmony_snn"), verbose=FALSE)
log_info(sprintf("Graphs now: %s", paste(Graphs(obj), collapse=", ")), stage=STAGE)
if (!all(pre$graphs %in% Graphs(obj))) fail("A Phase 1 graph was lost.")

# ---- UMAP ------------------------------------------------------------------
log_info("RunUMAP on the Harmony embedding (parameters matched to the Phase 1 baseline) ...", stage=STAGE)
set.seed(seed)
obj <- RunUMAP(obj, reduction=red, dims=dims, reduction.name="postint_umap_harmony",
               reduction.key="UMAPHARM_", seed.use=seed, verbose=FALSE)

# ---- resolution sweep ------------------------------------------------------
log_info(sprintf("FindClusters sweep: Louvain (algorithm 1), resolutions %s", paste(resolutions, collapse=", ")), stage=STAGE)
set.seed(seed)
obj <- FindClusters(obj, graph.name="postint_harmony_snn", resolution=resolutions,
                    algorithm=1, random.seed=seed, verbose=FALSE)
res_cols <- sprintf("postint_harmony_clusters_res_%s", format(resolutions, nsmall=1))
auto_cols <- sprintf("postint_harmony_snn_res.%s", resolutions)
for (j in seq_along(resolutions)) {
  ac <- auto_cols[j]
  if (!(ac %in% colnames(obj@meta.data))) fail(sprintf("Expected Seurat column '%s' not created.", ac))
  obj[[res_cols[j]]] <- as.character(obj@meta.data[[ac]])
  log_info(sprintf("  res %.1f -> %d clusters (stored as %s)", resolutions[j],
                   length(unique(obj@meta.data[[ac]])), res_cols[j]), stage=STAGE)
}
# restore the legacy column and identities that FindClusters clobbered
if (!is.null(legacy_seurat_clusters)) obj@meta.data$seurat_clusters <- legacy_seurat_clusters
Idents(obj) <- legacy_idents

# ---- preservation proof ----------------------------------------------------
post_cols <- colnames(obj@meta.data)
checks <- list(
  pca_unchanged      = identical(digest(Embeddings(obj,"pca"), algo="md5"), pre$pca),
  umap_pre_unchanged = identical(digest(Embeddings(obj,"umap_preintegration"), algo="md5"), pre$umap_pre),
  harmony_unchanged  = identical(digest(Embeddings(obj, red), algo="md5"), pre$harmony),
  cells_unchanged    = identical(digest(colnames(obj), algo="md5"), pre$cells),
  phase1_graphs_kept = all(pre$graphs %in% Graphs(obj)),
  phase1_reductions_kept = all(pre$reductions %in% Reductions(obj)),
  assays_unchanged   = identical(pre$assays, sort(Assays(obj))),
  default_assay_unchanged = identical(pre$default_assay, DefaultAssay(obj)),
  legacy_seurat_clusters_restored = is.null(legacy_seurat_clusters) ||
    identical(as.character(obj@meta.data$seurat_clusters), as.character(legacy_seurat_clusters)),
  original_meta_cols_all_present = all(pre$meta_cols %in% post_cols)
)
# Per-column value check. A whole-frame digest is too blunt: it also reacts to
# attribute-level changes that carry no data. Compare each original column's values
# and class individually so any real change is named.
changed_cols <- character(0); retyped_cols <- character(0)
for (cn in pre$meta_cols) {
  now <- obj@meta.data[[cn]]
  if (!identical(digest(now, algo="md5"), pre$meta_percol[[cn]]$d)) {
    if (identical(as.character(now), as.character(NULL)) ||
        !identical(as.character(now), as.character(pre$meta_percol[[cn]]$d))) {
      # fall through to the explicit value comparison below
    }
    changed_cols <- c(changed_cols, cn)
  }
  if (!identical(class(now)[1], pre$meta_percol[[cn]]$cls)) retyped_cols <- c(retyped_cols, cn)
}
# Separate genuine value changes from pure type/attribute changes.
value_changed <- character(0)
for (cn in changed_cols) {
  if (!identical(as.character(obj@meta.data[[cn]]), as.character(orig_meta_backup[[cn]])))
    value_changed <- c(value_changed, cn)
}
checks$original_meta_values_unchanged <- length(value_changed) == 0
checks$original_meta_types_unchanged  <- length(retyped_cols) == 0
for (nm in names(checks)) log_info(sprintf("  %-34s : %s", nm, checks[[nm]]), stage=STAGE)
if (length(changed_cols) > 0)
  log_info(sprintf("  columns whose stored representation changed: %s", paste(changed_cols, collapse=", ")), stage=STAGE)
if (length(retyped_cols) > 0)
  wrn(sprintf("Original metadata columns changed R type (values intact): %s. Restoring original types.", paste(retyped_cols, collapse=", ")))
# Restore original types/values verbatim for any original column Seurat touched.
for (cn in unique(c(changed_cols, retyped_cols))) obj@meta.data[[cn]] <- orig_meta_backup[[cn]]
recheck <- vapply(pre$meta_cols, function(cn)
  identical(digest(obj@meta.data[[cn]], algo="md5"), pre$meta_percol[[cn]]$d), logical(1))
checks$original_metadata_restored_verbatim <- all(recheck)
log_info(sprintf("  %-34s : %s", "original_metadata_restored_verbatim", checks$original_metadata_restored_verbatim), stage=STAGE)
if (length(value_changed) > 0)
  fail(paste("Preservation failed - original metadata VALUES changed in:", paste(value_changed, collapse=", ")))
if (!all(unlist(checks))) fail(paste("Preservation failed:", paste(names(checks)[!unlist(checks)], collapse=", ")))
log_info("PASS - all Phase 1 / M11 state preserved (including the legacy seurat_clusters column).", stage=STAGE)

# ---- canonical program scores (M12 sanity-check panels, recomputed identically) ----
# Used only to (a) flag the M12 caveat compartments C1/C2 and (b) provide a
# compartment-coverage criterion for resolution selection. NOT an annotation.
program_sets <- list(
  Panleukocyte=c("PTPRC"),
  T_NK=c("CD3D","CD3E","CD2","IL7R","TRAC","NKG7","GNLY","KLRD1"),
  Myeloid=c("LYZ","CD68","AIF1","CSF1R","ITGAM","C1QA","C1QB","FCGR3A"),
  B_Plasma=c("MS4A1","CD79A","CD79B","JCHAIN","MZB1"),
  Mast=c("TPSAB1","TPSB2","CPA3","MS4A2"),
  Endothelial=c("PECAM1","VWF","CDH5","CLDN5"),
  Fibroblast=c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA"),
  Mural=c("ACTA2","RGS5","PDGFRB","MYH11","NOTCH3"),
  SchwannNC=c("SOX10","S100B","PLP1","MPZ","NGFR","PMP22","ERBB3"),
  Proliferation=c("MKI67","TOP2A","CCNB1","CDK1","PCNA"))
avail <- rownames(obj)
usable <- lapply(program_sets, function(g) intersect(g, avail))
usable <- usable[vapply(usable, length, 1L) >= 1L]
nbin <- 24L; ctrl <- min(100L, max(5L, as.integer(floor(length(avail)/nbin/2))))
if (ctrl != 100L) wrn(sprintf("AddModuleScore ctrl lowered to %d (object has %d features).", ctrl, length(avail)))
set.seed(seed)
obj <- AddModuleScore(obj, features=usable, name="m12prog_", seed=seed,
                      assay=DefaultAssay(obj), search=FALSE, nbin=nbin, ctrl=ctrl)
sc <- paste0("m12prog_", seq_along(usable))
for (j in seq_along(usable)) obj[[paste0("m12_progscore_", names(usable)[j])]] <- obj@meta.data[[sc[j]]]
zm <- scale(as.matrix(obj@meta.data[, sc, drop=FALSE])); colnames(zm) <- names(usable)
bi <- max.col(zm, ties.method="first"); bz <- zm[cbind(seq_len(nrow(zm)), bi)]
prog_label <- ifelse(bz >= 1, colnames(zm)[bi], "Unassigned")
obj$m12_program_argmax <- prog_label
for (cc in sc) obj[[cc]] <- NULL
log_info(sprintf("Program sanity-check labels: %s", paste(sprintf("%s=%d", names(table(prog_label)), as.integer(table(prog_label))), collapse=", ")), stage=STAGE)

# ---- extract for metric work ----------------------------------------------
emb <- Embeddings(obj, red)[, dims, drop=FALSE]
umap_xy <- Embeddings(obj, "postint_umap_harmony")
md <- obj@meta.data
grp <- as.character(md$sample_id); lv <- sort(unique(grp)); n_cells <- nrow(md)
clu <- lapply(res_cols, function(c) as.character(md[[c]])); names(clu) <- res_cols

# ---- per-resolution diagnostics -------------------------------------------
log_info("Computing per-resolution diagnostics ...", stage=STAGE)
log_info("Computing the Harmony-space distance matrix once for all silhouettes ...", stage=STAGE)
dmat <- dist(emb)

tiny_thresh <- max(50L, as.integer(ceiling(0.005 * n_cells)))
dom_thresh <- 0.60
size_rows <- list(); comp_rows <- list(); res_rows <- list()
for (j in seq_along(resolutions)) {
  r <- resolutions[j]; cl <- clu[[j]]; f <- factor(cl, levels = sort(unique(as.numeric(cl))))
  tb <- table(cl); nk <- length(tb)
  sil <- cluster::silhouette(as.integer(factor(cl)), dmat)
  sil_mean <- mean(sil[,3])
  # composition and flags
  ct <- table(cl, grp)
  dom_share <- apply(ct, 1, function(x) max(x)/sum(x))
  dom_sample <- lv[apply(ct, 1, which.max)]; names(dom_sample) <- rownames(ct)
  n_tiny <- sum(tb < tiny_thresh); n_dom <- sum(dom_share > dom_thresh)
  # compartment coverage: distinct canonical programs that are the modal label of >=1 cluster
  modal_prog <- vapply(names(tb), function(k){ t2 <- table(prog_label[cl==k]); names(t2)[which.max(t2)] }, character(1))
  cov <- length(setdiff(unique(modal_prog), "Unassigned"))
  # purity w.r.t. program labels (assigned cells only)
  keep <- prog_label != "Unassigned"
  purity <- if (sum(keep)>0) mean(vapply(split(prog_label[keep], cl[keep]), function(x) max(table(x))/length(x), numeric(1))) else NA_real_
  # per cluster rows
  for (k in names(tb)) {
    size_rows[[length(size_rows)+1]] <- data.frame(resolution=r, cluster=k, n_cells=as.integer(tb[[k]]),
      pct_of_all=100*as.integer(tb[[k]])/n_cells, dominant_sample=dom_sample[[k]],
      dominant_sample_share=as.numeric(dom_share[[k]]),
      is_tiny=as.integer(tb[[k]]) < tiny_thresh, is_sample_dominated=as.numeric(dom_share[[k]]) > dom_thresh,
      modal_canonical_program=modal_prog[[k]],
      mean_silhouette=mean(sil[cl==k,3]), stringsAsFactors=FALSE)
    for (s in lv) comp_rows[[length(comp_rows)+1]] <- data.frame(resolution=r, cluster=k, sample_id=s,
      n_cells=as.integer(ct[k,s]), pct_of_cluster=100*as.integer(ct[k,s])/sum(ct[k,]), stringsAsFactors=FALSE)
  }
  res_rows[[j]] <- data.frame(resolution=r, n_clusters=nk,
    min_size=as.integer(min(tb)), median_size=as.numeric(median(tb)), max_size=as.integer(max(tb)),
    n_tiny_clusters=n_tiny, tiny_threshold=tiny_thresh,
    n_sample_dominated_clusters=n_dom, sample_dominance_threshold=dom_thresh,
    mean_silhouette=sil_mean, canonical_program_coverage=cov, mean_program_purity=purity,
    stringsAsFactors=FALSE)
  log_info(sprintf("  res %.1f: k=%2d  sizes %d/%.0f/%d  tiny=%d  sample-dominated=%d  sil=%.4f  coverage=%d  purity=%.3f",
                   r, as.integer(nk), as.integer(min(tb)), as.numeric(median(tb)), as.integer(max(tb)),
                   as.integer(n_tiny), as.integer(n_dom), sil_mean, as.integer(cov), purity), stage=STAGE)
}
rm(dmat); invisible(gc(verbose=FALSE))
res_df <- do.call(rbind, res_rows)

# ---- stability: ARI between adjacent resolutions ---------------------------
log_info("Computing adjacent-resolution ARI stability ...", stage=STAGE)
ari_adj <- rep(NA_real_, length(resolutions))
for (j in seq_along(resolutions)) {
  vals <- c()
  if (j > 1) vals <- c(vals, aricode::ARI(clu[[j]], clu[[j-1]]))
  if (j < length(resolutions)) vals <- c(vals, aricode::ARI(clu[[j]], clu[[j+1]]))
  ari_adj[j] <- mean(vals)
}
res_df$mean_adjacent_ARI <- ari_adj
ari_pairs <- do.call(rbind, lapply(seq_len(length(resolutions)-1), function(j)
  data.frame(res_from=resolutions[j], res_to=resolutions[j+1], ARI=aricode::ARI(clu[[j]], clu[[j+1]]),
             NMI=aricode::NMI(clu[[j]], clu[[j+1]]), stringsAsFactors=FALSE)))
write.table(ari_pairs, file.path(out_dir,"resolution_transition_ari.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# ---- composite ranking and recommendation ----------------------------------
# Documented, deterministic rule. Each criterion is min-max scaled across the ten
# resolutions and combined with fixed weights; no criterion is UMAP appearance.
nz <- function(x){ r <- range(x, na.rm=TRUE); if (diff(r)==0) rep(0.5,length(x)) else (x-r[1])/diff(r) }
res_df$score_stability   <- nz(res_df$mean_adjacent_ARI)          # higher better
res_df$score_compactness <- nz(res_df$mean_silhouette)            # higher better
res_df$score_coverage    <- nz(res_df$canonical_program_coverage) # higher better
res_df$score_purity      <- nz(res_df$mean_program_purity)        # higher better
res_df$score_fragmentation <- 1 - nz(res_df$n_tiny_clusters)      # fewer tiny better
# annotation workability: prefer 10-30 clusters; penalise outside that band
res_df$score_granularity <- 1 - nz(pmax(0, 10 - res_df$n_clusters) + pmax(0, res_df$n_clusters - 30))
res_df$composite_score <- with(res_df,
  0.30*score_stability + 0.20*score_compactness + 0.15*score_coverage +
  0.15*score_purity + 0.10*score_fragmentation + 0.10*score_granularity)
res_df <- res_df[order(-res_df$composite_score), ]
primary_res <- res_df$resolution[1]
# alternative: the highest-scoring resolution that yields a materially different k
alt_pool <- res_df[abs(res_df$n_clusters - res_df$n_clusters[1]) >= 2, ]
alternative_res <- if (nrow(alt_pool) > 0) alt_pool$resolution[1] else res_df$resolution[2]
res_df$recommendation <- ifelse(res_df$resolution==primary_res, "PRIMARY",
                         ifelse(res_df$resolution==alternative_res, "ALTERNATIVE", ""))
write.table(res_df, file.path(out_dir,"clustering_comparison_table.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info(sprintf("RECOMMENDED PRIMARY resolution = %.1f (%d clusters); ALTERNATIVE = %.1f (%d clusters)",
  primary_res, res_df$n_clusters[res_df$resolution==primary_res],
  alternative_res, res_df$n_clusters[res_df$resolution==alternative_res]), stage=STAGE)

pc <- sprintf("postint_harmony_clusters_res_%s", format(primary_res, nsmall=1))
ac <- sprintf("postint_harmony_clusters_res_%s", format(alternative_res, nsmall=1))
obj$postint_harmony_primary_cluster <- obj@meta.data[[pc]]
obj$postint_harmony_alternative_cluster <- obj@meta.data[[ac]]
obj$postint_primary_resolution <- as.character(primary_res)
obj$postint_alternative_resolution <- as.character(alternative_res)

size_df <- do.call(rbind, size_rows); comp_df <- do.call(rbind, comp_rows)
write.table(size_df, file.path(out_dir,"cluster_sizes_all_resolutions.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(comp_df, file.path(out_dir,"cluster_composition_by_sample_all_resolutions.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(size_df[size_df$resolution==primary_res,], file.path(out_dir,"cluster_sizes_primary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(comp_df[comp_df$resolution==primary_res,], file.path(out_dir,"cluster_composition_primary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# M12 caveat carry-forward: flag clusters enriched for the B/plasma (C1) or fibroblast (C2) programs
pcl <- as.character(md[[pc]])
flag_rows <- do.call(rbind, lapply(sort(unique(pcl)), function(k){
  sel <- pcl==k; tb <- table(factor(prog_label[sel], levels=c(names(usable),"Unassigned")))
  data.frame(cluster=k, n_cells=sum(sel),
    dominant_sample=lv[which.max(table(factor(grp[sel], levels=lv)))],
    dominant_sample_share=max(table(factor(grp[sel], levels=lv)))/sum(sel),
    frac_B_Plasma=as.numeric(tb["B_Plasma"])/sum(sel),
    frac_Fibroblast=as.numeric(tb["Fibroblast"])/sum(sel),
    frac_SchwannNC=as.numeric(tb["SchwannNC"])/sum(sel),
    flag_C1_B_plasma_caveat=as.numeric(tb["B_Plasma"])/sum(sel) > 0.20,
    flag_C2_fibroblast_caveat=as.numeric(tb["Fibroblast"])/sum(sel) > 0.20,
    flag_sample_dominated=max(table(factor(grp[sel], levels=lv)))/sum(sel) > dom_thresh,
    stringsAsFactors=FALSE)}))
write.table(flag_rows, file.path(out_dir,"primary_cluster_m12_caveat_flags.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info(sprintf("M12 caveat flags on the primary solution: C1(B/plasma)=%d clusters, C2(fibroblast)=%d, sample-dominated=%d",
  sum(flag_rows$flag_C1_B_plasma_caveat), sum(flag_rows$flag_C2_fibroblast_caveat), sum(flag_rows$flag_sample_dominated)), stage=STAGE)

# ---- FIGURES ---------------------------------------------------------------
log_info("Generating M13 figures ...", stage=STAGE)
th <- theme_bw(base_size=13) + theme(panel.grid.minor=element_blank(),
      plot.title=element_text(face="bold"), plot.subtitle=element_text(size=10,colour="grey30"))
set.seed(seed); ord <- sample.int(n_cells)
D <- data.frame(x=umap_xy[ord,1], y=umap_xy[ord,2], sample_id=grp[ord], program=prog_label[ord], stringsAsFactors=FALSE)
for (j in seq_along(res_cols)) D[[res_cols[j]]] <- clu[[j]][ord]
D$primary <- factor(as.character(md[[pc]])[ord], levels=sort(unique(as.numeric(md[[pc]]))))
D$alternative <- factor(as.character(md[[ac]])[ord], levels=sort(unique(as.numeric(md[[ac]]))))
pal_sample <- setNames(RColorBrewer::brewer.pal(max(3,length(lv)),"Set1")[seq_along(lv)], lv)
big_pal <- function(n) grDevices::colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(n)

frows <- list()
sf <- function(p, name, w, h, method, params) {
  for (ext in c("pdf","png")) { fp <- file.path(fig_dir, paste0(name,".",ext))
    ggsave(fp, p, width=w, height=h, units="in", dpi=200, device=ext)
    frows[[length(frows)+1]] <<- data.frame(figure_path=fp, phase="phase2", milestone="M13",
      dataset="combined", processing_stage="post_harmony_clustering", analysis_method=method,
      parameters=params, input_object_checksum=in_md5,
      generating_script="scripts/R/phase2/cluster_sweep_harmony.R", stringsAsFactors=FALSE) }
  log_info(sprintf("  wrote %s", name), stage=STAGE)
}
umap_disc <- function(df, col, title, sub, pal=NULL, label=TRUE, ptsize=0.3) {
  p <- ggplot(df, aes(x,y,colour=.data[[col]])) + geom_point(size=ptsize, alpha=0.8, stroke=0) +
    labs(title=title, subtitle=sub, x="UMAP 1", y="UMAP 2", colour=col) + th +
    guides(colour=guide_legend(override.aes=list(size=3.5,alpha=1), ncol=1))
  if (!is.null(pal)) p <- p + scale_colour_manual(values=pal) else
    p <- p + scale_colour_manual(values=big_pal(length(unique(df[[col]]))))
  if (label) { cen <- aggregate(cbind(x,y) ~ get(col), data=df, FUN=median); names(cen)[1] <- "lab"
    p <- p + ggplot2::annotate("text", x=cen$x, y=cen$y, label=cen$lab, size=3.6, fontface="bold", colour="black") }
  p
}
sub_primary <- sprintf("Harmony embedding (%s dims 1:%d) -> UMAP -> Louvain k.param=%d, resolution %.1f, %d clusters, %d cells",
                       red, n_dims, k_param, primary_res, nlevels(D$primary), n_cells)
sf(umap_disc(D,"primary","Post-Harmony UMAP - primary clustering", sub_primary),
   "M13_01_harmony_umap_clusters", 11, 9, "UMAP_clusters",
   sprintf("res=%.1f;dims=1:%d;k.param=%d;algorithm=Louvain;seed=%d", primary_res, n_dims, k_param, seed))
sf(umap_disc(D,"alternative","Post-Harmony UMAP - alternative clustering",
   sprintf("resolution %.1f, %d clusters", alternative_res, nlevels(D$alternative))),
   "M13_02_harmony_umap_clusters_alternative", 11, 9, "UMAP_clusters",
   sprintf("res=%.1f;dims=1:%d;k.param=%d", alternative_res, n_dims, k_param))
sf(umap_disc(D,"sample_id","Post-Harmony UMAP by sample_id",
   "sample_id is simultaneously the dataset, the patient and the Harmony grouping variable.",
   pal=pal_sample, label=FALSE),
   "M13_03_harmony_umap_by_sample", 10, 8, "UMAP_metadata", sprintf("dims=1:%d", n_dims))
sf(umap_disc(D,"program","Post-Harmony UMAP by canonical program (sanity check, not annotation)",
   "Dominant canonical broad-lineage program per cell (z >= 1). Formal annotation is M15.", label=FALSE),
   "M13_04_harmony_umap_by_canonical_program", 10.5, 8, "UMAP_metadata", "argmax_z_program(threshold=1)")
# per-resolution UMAPs
for (j in seq_along(resolutions)) {
  r <- resolutions[j]; cn <- res_cols[j]
  D2 <- D; D2$.cl <- factor(D[[cn]], levels=sort(unique(as.numeric(D[[cn]]))))
  sf(umap_disc(D2,".cl", sprintf("Post-Harmony UMAP - resolution %.1f", r),
      sprintf("%d clusters", nlevels(D2$.cl))),
     sprintf("M13_10_umap_resolution_%s", format(r, nsmall=1)), 9, 7.5, "UMAP_clusters",
     sprintf("res=%.1f;dims=1:%d;k.param=%d", r, n_dims, k_param))
}
# combined sweep panel
panels <- lapply(seq_along(resolutions), function(j){ r <- resolutions[j]
  d2 <- D; d2$.cl <- factor(D[[res_cols[j]]], levels=sort(unique(as.numeric(D[[res_cols[j]]]))))
  ggplot(d2, aes(x,y,colour=.cl)) + geom_point(size=0.10, alpha=0.7, stroke=0) +
    scale_colour_manual(values=big_pal(nlevels(d2$.cl))) +
    labs(title=sprintf("res %.1f  (%d clusters)", r, nlevels(d2$.cl)), x=NULL, y=NULL) +
    theme_bw(base_size=10) + theme(legend.position="none", panel.grid=element_blank(),
      axis.text=element_blank(), axis.ticks=element_blank()) })
sf(wrap_plots(panels, ncol=5) + plot_annotation(title="Resolution sweep on the Harmony UMAP (0.1 - 1.0)",
     subtitle="All ten solutions are preserved in the object as postint_harmony_clusters_res_*"),
   "M13_11_resolution_sweep_panel", 18, 8, "UMAP_panel_grid",
   sprintf("resolutions=%s;dims=1:%d", paste(resolutions,collapse=","), n_dims))
# cluster sizes
sp <- size_df[size_df$resolution==primary_res,]; sp$cluster <- factor(sp$cluster, levels=sp$cluster[order(-sp$n_cells)])
sf(ggplot(sp, aes(cluster, n_cells, fill=is_tiny)) + geom_col() +
     geom_text(aes(label=n_cells), vjust=-0.3, size=2.9) +
     scale_fill_manual(values=c("FALSE"="#3B6EA5","TRUE"="#C1452B"), name=sprintf("< %d cells", tiny_thresh)) +
     labs(title=sprintf("Cluster sizes at the primary resolution (%.1f)", primary_res),
          subtitle=sprintf("%d clusters, %d cells; tiny-cluster threshold = %d cells (0.5%%)", nrow(sp), n_cells, tiny_thresh),
          x="cluster", y="cells") + th,
   "M13_05_cluster_sizes_primary", 12, 6, "barplot", sprintf("res=%.1f", primary_res))
sf(ggplot(sp, aes(cluster, pct_of_all)) + geom_col(fill="#3B6EA5") +
     labs(title=sprintf("Cluster proportions at the primary resolution (%.1f)", primary_res),
          x="cluster", y="% of all cells") + th,
   "M13_06_cluster_proportions_primary", 12, 5.5, "barplot", sprintf("res=%.1f", primary_res))
# composition
cp <- comp_df[comp_df$resolution==primary_res,]
cp$cluster <- factor(cp$cluster, levels=sp$cluster[order(as.numeric(as.character(sp$cluster)))])
sf(ggplot(cp, aes(cluster, pct_of_cluster, fill=sample_id)) + geom_col() +
     scale_fill_manual(values=pal_sample) + geom_hline(yintercept=60, linetype="dashed", colour="grey20") +
     labs(title=sprintf("Cluster composition by sample at resolution %.1f", primary_res),
          subtitle="Dashed line = the 60% single-sample dominance threshold. sample_id = dataset = patient.",
          x="cluster", y="% of cluster") + th,
   "M13_07_cluster_composition_by_sample_primary", 13, 6, "stacked_barplot", sprintf("res=%.1f", primary_res))
sf(ggplot(cp, aes(cluster, n_cells, fill=sample_id)) + geom_col() + scale_fill_manual(values=pal_sample) +
     labs(title=sprintf("Cluster composition by sample - absolute counts (resolution %.1f)", primary_res),
          x="cluster", y="cells") + th,
   "M13_08_cluster_composition_counts_primary", 13, 6, "stacked_barplot", sprintf("res=%.1f", primary_res))
# stability / transitions
sf(ggplot(ari_pairs, aes(factor(sprintf("%.1f->%.1f", res_from, res_to)), ARI, group=1)) +
     geom_line(colour="#3B6EA5") + geom_point(size=2.5, colour="#3B6EA5") +
     geom_text(aes(label=sprintf("%.3f", ARI)), vjust=-0.9, size=3) + ylim(0,1.05) +
     labs(title="Clustering stability between adjacent resolutions",
          subtitle="Adjusted Rand Index. High values indicate a stable plateau where the partition changes little.",
          x="resolution transition", y="ARI") + th,
   "M13_09_resolution_stability_ari", 11, 5.5, "line_metric", "metric=ARI_adjacent_resolutions")
tm <- as.data.frame(table(from=clu[[which(resolutions==primary_res)]], to=clu[[which(resolutions==alternative_res)]]))
tm$from <- factor(tm$from, levels=sort(unique(as.numeric(as.character(tm$from)))))
tm$to <- factor(tm$to, levels=sort(unique(as.numeric(as.character(tm$to)))))
sf(ggplot(tm, aes(from, to, fill=log10(Freq+1))) + geom_tile() +
     scale_fill_viridis_c(option="mako", name="log10(cells+1)") +
     labs(title=sprintf("Cluster transition matrix: resolution %.1f -> %.1f", primary_res, alternative_res),
          subtitle="Where cells move between the primary and alternative solutions",
          x=sprintf("primary (res %.1f)", primary_res), y=sprintf("alternative (res %.1f)", alternative_res)) + th,
   "M13_12_cluster_transition_matrix", 10, 8, "heatmap", sprintf("res %.1f vs %.1f", primary_res, alternative_res))
# ranked comparison
rl <- res_df[order(res_df$resolution),]
sf((ggplot(rl, aes(resolution, n_clusters)) + geom_line(colour="#3B6EA5") + geom_point() +
      labs(title="Clusters per resolution", y="clusters") + th) /
   (ggplot(rl, aes(resolution, mean_adjacent_ARI)) + geom_line(colour="#2E8B57") + geom_point() +
      labs(title="Mean adjacent-resolution ARI (stability)", y="ARI") + th) /
   (ggplot(rl, aes(resolution, mean_silhouette)) + geom_line(colour="#C1452B") + geom_point() +
      labs(title="Mean cluster silhouette in the Harmony space (compactness)", y="silhouette") + th) /
   (ggplot(rl, aes(resolution, composite_score)) + geom_col(aes(fill=resolution==primary_res)) +
      scale_fill_manual(values=c("FALSE"="grey70","TRUE"="#C1452B"), guide="none") +
      labs(title=sprintf("Composite selection score (primary = %.1f)", primary_res), y="score") + th) +
   plot_annotation(title="Resolution selection criteria",
     subtitle="Composite = 0.30 stability + 0.20 compactness + 0.15 compartment coverage + 0.15 program purity + 0.10 low fragmentation + 0.10 granularity. UMAP appearance is not a criterion."),
   "M13_13_resolution_selection_criteria", 10, 13, "multi_panel_metric", "composite selection rule")
write.table(do.call(rbind, frows), file.path(out_dir,"figure_index_m13.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# ---- save & round-trip ------------------------------------------------------
log_info(sprintf("Saving M13 object to %s ...", out_rds), stage=STAGE)
saveRDS(obj, out_rds)
out_md5 <- digest(out_rds, file=TRUE, algo="md5"); out_sha <- digest(out_rds, file=TRUE, algo="sha256")
out_size <- file.info(out_rds)$size
log_info(sprintf("Saved %.2f GB | md5 %s", out_size/1024^3, out_md5), stage=STAGE)
ref_umap <- digest(umap_xy, algo="md5")
rm(obj, emb); invisible(gc(verbose=FALSE))
o2 <- readRDS(out_rds)
v <- list(cells=ncol(o2)==n_cells,
  pca_present="pca" %in% Reductions(o2), umap_pre_present="umap_preintegration" %in% Reductions(o2),
  harmony_present=red %in% Reductions(o2), umap_harmony_present="postint_umap_harmony" %in% Reductions(o2),
  umap_harmony_identical=identical(digest(Embeddings(o2,"postint_umap_harmony"),algo="md5"), ref_umap),
  harmony_unchanged=identical(digest(Embeddings(o2,red),algo="md5"), pre$harmony),
  pca_unchanged=identical(digest(Embeddings(o2,"pca"),algo="md5"), pre$pca),
  graphs_present=all(c("SCT_nn","SCT_snn","postint_harmony_nn","postint_harmony_snn") %in% Graphs(o2)),
  all_res_cols_present=all(res_cols %in% colnames(o2@meta.data)),
  primary_col_present="postint_harmony_primary_cluster" %in% colnames(o2@meta.data),
  legacy_seurat_clusters_intact=is.null(legacy_seurat_clusters) ||
    identical(as.character(o2@meta.data$seurat_clusters), as.character(legacy_seurat_clusters)),
  no_nonfinite_umap=sum(!is.finite(Embeddings(o2,"postint_umap_harmony")))==0)
for (nm in names(v)) log_info(sprintf("  %-30s : %s", nm, v[[nm]]), stage=STAGE)
if (!all(unlist(v))) fail(paste("Round-trip validation failed:", paste(names(v)[!unlist(v)], collapse=", ")))
final_meta_cols <- colnames(o2@meta.data); final_reductions <- sort(Reductions(o2)); final_graphs <- sort(Graphs(o2))
rm(o2); invisible(gc(verbose=FALSE))

rec <- list(milestone="M13", phase="phase2", timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/R/phase2/cluster_sweep_harmony.R",
  input=list(path=input_rds, md5=in_md5),
  output=list(path=out_rds, md5=out_md5, sha256=out_sha, size_bytes=out_size, cells=n_cells,
              reductions=final_reductions, graphs=final_graphs, n_metadata_columns=length(final_meta_cols)),
  parameters=list(reduction=red, dims=sprintf("1:%d",n_dims), k_param=k_param, algorithm="Louvain (algorithm 1)",
    resolutions=resolutions, seed=seed, umap="RunUMAP defaults; dims 1:30; seed 42; n.neighbors 30; min.dist 0.3; cosine",
    cluster_columns=res_cols, primary_column="postint_harmony_primary_cluster",
    alternative_column="postint_harmony_alternative_cluster"),
  selection=list(primary_resolution=primary_res, alternative_resolution=alternative_res,
    rule="composite = 0.30*stability(ARI) + 0.20*compactness(silhouette) + 0.15*compartment coverage + 0.15*program purity + 0.10*low fragmentation + 0.10*granularity(10-30 clusters); all min-max scaled across the ten resolutions",
    n_primary_clusters=res_df$n_clusters[res_df$resolution==primary_res],
    n_alternative_clusters=res_df$n_clusters[res_df$resolution==alternative_res]),
  preservation_checks=checks, round_trip_validation=v,
  m12_caveat_flags=list(C1_B_plasma=sum(flag_rows$flag_C1_B_plasma_caveat),
    C2_fibroblast=sum(flag_rows$flag_C2_fibroblast_caveat),
    sample_dominated=sum(flag_rows$flag_sample_dominated)),
  warnings=warns, slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD", intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(), t_start, units="secs")),
  session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir,"m13_clustering_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase2_m13_clustering", inputs=list(m11_rds=input_rds),
  outputs=list(clustered_rds=out_rds, comparison=file.path(out_dir,"clustering_comparison_table.tsv"),
    sizes=file.path(out_dir,"cluster_sizes_all_resolutions.tsv"),
    composition=file.path(out_dir,"cluster_composition_by_sample_all_resolutions.tsv"),
    record=file.path(out_dir,"m13_clustering_record.json")),
  parameters=list(milestone="M13", primary_resolution=primary_res, alternative_resolution=alternative_res,
    k_param=k_param, dims=sprintf("1:%d",n_dims), seed=seed,
    slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, file.path(out_dir,"prov_m13_clustering.json"))
log_info(sprintf("M13 COMPLETE in %.1f min. Warnings: %d", as.numeric(difftime(Sys.time(),t_start,units="mins")), length(warns)), stage=STAGE)
log_system_usage(stage=STAGE)
