# scripts/R/phase2/discover_markers_harmony.R
#
# Phase 2 / Milestone M14 - cluster marker discovery on the M13 primary (and alternative)
# post-Harmony clustering.
#
# SCIENTIFIC SCOPE: these are exploratory CLUSTER-CHARACTERISATION markers used to support
# annotation in M15. They are NOT condition-level differential expression. Cells are not
# independent biological replicates, and no between-sample / between-patient /
# between-condition claim may be made from these tables.
#
# Harmony coordinates are NEVER used for testing. Testing uses the SCT expression assay.
# Seurat v5 layer handling is inspected explicitly, not assumed: the SCT assay carries
# FOUR SCTransform models (M10 Finding 1), so PrepSCTFindMarkers() is mandatory.

options(stringsAsFactors = FALSE); options(future.globals.maxSize = +Inf)
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(ggplot2); library(patchwork); library(viridis); library(RColorBrewer); library(dplyr)
  library(jsonlite); library(digest)
})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase2_m14_markers"

args <- commandArgs(trailingOnly = TRUE)
input_rds <- "results/phase2/clustering/phase2_harmony_clustered.rds"
out_dir   <- "results/phase2/markers"
fig_dir   <- "results/phase2/figures/M14"
assay_use <- "SCT"; layer_use <- "data"
test_use  <- "wilcox"; min_pct <- 0.25; logfc_thr <- 0.25; only_pos <- TRUE
seed <- 42L; expected_cells <- 19716L; validation_mode <- FALSE
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--input"){input_rds<-args[i+1];i<-i+2} else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2}
  else if (a=="--fig-dir"){fig_dir<-args[i+1];i<-i+2} else if (a=="--assay"){assay_use<-args[i+1];i<-i+2}
  else if (a=="--min-pct"){min_pct<-as.numeric(args[i+1]);i<-i+2}
  else if (a=="--logfc"){logfc_thr<-as.numeric(args[i+1]);i<-i+2}
  else if (a=="--random-seed"){seed<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--expected-cells"){expected_cells<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--validation-mode"){validation_mode<-TRUE;i<-i+1}
  else stop(sprintf("Unknown argument: %s", a)) }

t_start <- Sys.time(); warns <- character(0)
wrn <- function(m){ warns <<- c(warns,m); log_warn(m, stage=STAGE) }
fail <- function(m){ log_error(m, stage=STAGE); stop(m, call.=FALSE) }
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE); dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
set.seed(seed)

log_info("=================== M14 MARKER DISCOVERY ===================", stage=STAGE)
in_md5 <- digest(input_rds, file=TRUE, algo="md5")
log_info(sprintf("Input %s (md5 %s)", input_rds, in_md5), stage=STAGE)
obj <- readRDS(input_rds)
log_info(sprintf("Loaded: %d features x %d cells", nrow(obj), ncol(obj)), stage=STAGE)
if (!validation_mode && ncol(obj) != expected_cells) fail("Cell count mismatch.")

# ---------------------------------------------------------------------------
# Explicit Seurat v5 assay / layer inspection - do not guess
# ---------------------------------------------------------------------------
log_info("Inspecting assay and layer structure before any testing ...", stage=STAGE)
assay_rows <- list()
for (a in Assays(obj)) {
  ao <- obj[[a]]
  lys <- tryCatch(SeuratObject::Layers(ao), error=function(e) NA_character_)
  nm <- if (inherits(ao,"SCTAssay")) length(levels(ao)) else NA_integer_
  assay_rows[[a]] <- data.frame(assay=a, class=class(ao)[1], is_default=identical(a,DefaultAssay(obj)),
    n_features=nrow(ao), n_cells=ncol(ao), layers=paste(lys, collapse=";"),
    n_sct_models=nm, stringsAsFactors=FALSE)
  log_info(sprintf("  %s (%s%s): layers = %s%s", a, class(ao)[1],
    if (identical(a,DefaultAssay(obj))) ", DEFAULT" else "", paste(lys, collapse=";"),
    if (!is.na(nm)) sprintf(" | SCT models = %d", nm) else ""), stage=STAGE)
}
write.table(do.call(rbind, assay_rows), file.path(out_dir,"assay_layer_inspection.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

DefaultAssay(obj) <- assay_use
n_models <- if (inherits(obj[[assay_use]],"SCTAssay")) length(levels(obj[[assay_use]])) else NA_integer_
prep_needed <- isTRUE(!is.na(n_models) && n_models > 1)
if (prep_needed) {
  log_info(sprintf("SCT assay carries %d models -> PrepSCTFindMarkers() is REQUIRED and will be run.", n_models), stage=STAGE)
  t0 <- Sys.time(); set.seed(seed)
  obj <- PrepSCTFindMarkers(obj, assay=assay_use, verbose=FALSE)
  log_info(sprintf("PrepSCTFindMarkers completed in %.1f s (counts recorrected to a common sequencing depth).",
                   as.numeric(difftime(Sys.time(),t0,units="secs"))), stage=STAGE)
} else log_info("Single SCT model - PrepSCTFindMarkers not required.", stage=STAGE)

primary_col <- "postint_harmony_primary_cluster"
alt_col     <- "postint_harmony_alternative_cluster"
if (!(primary_col %in% colnames(obj@meta.data))) fail("Primary cluster column missing - run M13 first.")
primary_res <- unique(obj$postint_primary_resolution)[1]
alt_res     <- unique(obj$postint_alternative_resolution)[1]
ord_lv <- function(x) sort(unique(as.numeric(as.character(x))))
obj$.m14_primary <- factor(as.character(obj@meta.data[[primary_col]]), levels=ord_lv(obj@meta.data[[primary_col]]))
obj$.m14_alt     <- factor(as.character(obj@meta.data[[alt_col]]),     levels=ord_lv(obj@meta.data[[alt_col]]))
log_info(sprintf("Primary resolution %s -> %d clusters; alternative %s -> %d clusters",
                 primary_res, nlevels(obj$.m14_primary), alt_res, nlevels(obj$.m14_alt)), stage=STAGE)

marker_params <- list(assay=assay_use, layer=layer_use,
  normalisation="SCTransform residuals; PrepSCTFindMarkers applied (4 models)",
  test.use=test_use, min.pct=min_pct, logfc.threshold=logfc_thr, only.pos=only_pos,
  recorrect_umi=FALSE, random.seed=seed,
  presto_accelerated=requireNamespace("presto", quietly=TRUE),
  scope="cluster characterisation for annotation - NOT condition-level differential expression")
write.table(data.frame(parameter=names(marker_params),
  value=vapply(marker_params, function(x) paste(as.character(x), collapse=";"), character(1))),
  file.path(out_dir,"marker_parameters.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

run_markers <- function(ident_col, tag) {
  log_info(sprintf("FindAllMarkers on %s (%s) ...", ident_col, tag), stage=STAGE)
  Idents(obj) <- obj@meta.data[[ident_col]]
  t0 <- Sys.time(); set.seed(seed)
  m <- FindAllMarkers(obj, assay=assay_use, slot=layer_use, test.use=test_use,
                      min.pct=min_pct, logfc.threshold=logfc_thr, only.pos=only_pos,
                      recorrect_umi=FALSE, random.seed=seed, verbose=FALSE)
  log_info(sprintf("  %d marker rows across %d clusters in %.1f s", nrow(m),
                   length(unique(m$cluster)), as.numeric(difftime(Sys.time(),t0,units="secs"))), stage=STAGE)
  m
}
mk_primary <- run_markers(".m14_primary", sprintf("primary, res %s", primary_res))
mk_alt     <- run_markers(".m14_alt",     sprintf("alternative, res %s", alt_res))

tidy_markers <- function(m, res_label) {
  m$cluster <- as.character(m$cluster)
  m <- m[, c("cluster","gene","avg_log2FC","pct.1","pct.2","p_val","p_val_adj")]
  m$pct_diff <- m$pct.1 - m$pct.2
  m$resolution <- res_label
  m[order(m$cluster, -m$avg_log2FC), ]
}
mk_primary <- tidy_markers(mk_primary, primary_res)
mk_alt     <- tidy_markers(mk_alt, alt_res)
write.table(mk_primary, file.path(out_dir,"all_cluster_markers.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(mk_alt, file.path(out_dir,"all_cluster_markers_alternative_resolution.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

sig <- mk_primary[mk_primary$p_val_adj < 0.05 & mk_primary$avg_log2FC >= logfc_thr, ]
write.table(sig, file.path(out_dir,"filtered_cluster_markers.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
top_n <- function(df, n) do.call(rbind, lapply(split(df, df$cluster), function(d) head(d[order(-d$avg_log2FC),], n)))
for (n in c(10,20,50)) {
  tn <- top_n(sig, n)
  write.table(tn, file.path(out_dir, sprintf("top%d_markers_per_cluster.tsv", n)), sep="\t", row.names=FALSE, quote=FALSE)
}
top10 <- top_n(sig, 10)
write.table(top10, file.path(out_dir,"top_markers_per_cluster.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

grp <- as.character(obj$sample_id); lv <- sort(unique(grp))
pcl <- as.character(obj@meta.data[[primary_col]])
prog <- if ("m12_program_argmax" %in% colnames(obj@meta.data)) as.character(obj$m12_program_argmax) else rep(NA_character_, ncol(obj))
summ <- do.call(rbind, lapply(levels(obj$.m14_primary), function(k) {
  sel <- pcl == k; s <- sig[sig$cluster == k, ]
  ct <- table(factor(grp[sel], levels=lv))
  pt <- table(prog[sel]); mp <- if (length(pt)) names(pt)[which.max(pt)] else NA_character_
  data.frame(cluster=k, n_cells=sum(sel), pct_of_all=100*sum(sel)/length(sel),
    n_significant_markers=nrow(s),
    top10_markers=paste(head(s$gene[order(-s$avg_log2FC)],10), collapse=","),
    max_avg_log2FC=if (nrow(s)) max(s$avg_log2FC) else NA_real_,
    dominant_sample=lv[which.max(ct)], dominant_sample_pct=100*max(ct)/sum(ct),
    modal_canonical_program=mp,
    stringsAsFactors=FALSE) }))
write.table(summ, file.path(out_dir,"marker_summary_by_cluster.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info("Per-cluster marker summary:", stage=STAGE)
for (r in seq_len(nrow(summ)))
  log_info(sprintf("  C%-3s n=%5d  markers=%4d  %-14s %.0f%% %s | %s", summ$cluster[r], summ$n_cells[r],
    summ$n_significant_markers[r], summ$dominant_sample[r], summ$dominant_sample_pct[r],
    summ$modal_canonical_program[r], substr(summ$top10_markers[r],1,70)), stage=STAGE)

# ---------------------------------------------------------------------------
# FIGURES
# ---------------------------------------------------------------------------
log_info("Generating M14 figures ...", stage=STAGE)
frows <- list()
sf <- function(p, name, w, h, method, params) {
  for (ext in c("pdf","png")) { fp <- file.path(fig_dir, paste0(name,".",ext))
    ggsave(fp, p, width=w, height=h, units="in", dpi=200, device=ext, limitsize=FALSE)
    frows[[length(frows)+1]] <<- data.frame(figure_path=fp, phase="phase2", milestone="M14",
      dataset="combined", processing_stage="cluster_marker_discovery", analysis_method=method,
      parameters=params, input_object_checksum=in_md5,
      generating_script="scripts/R/phase2/discover_markers_harmony.R", stringsAsFactors=FALSE) }
  log_info(sprintf("  wrote %s", name), stage=STAGE)
}
th <- theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"))
mp <- sprintf("assay=%s;layer=%s;test=%s;min.pct=%.2f;logfc=%.2f;only.pos=%s;res=%s",
              assay_use, layer_use, test_use, min_pct, logfc_thr, only_pos, primary_res)

# cluster UMAP (labelled)
uxy <- Embeddings(obj, "postint_umap_harmony")
set.seed(seed); ordv <- sample.int(ncol(obj))
D <- data.frame(x=uxy[ordv,1], y=uxy[ordv,2], cl=factor(pcl[ordv], levels=levels(obj$.m14_primary)))
bigpal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(nlevels(D$cl))
cen <- aggregate(cbind(x,y) ~ cl, data=D, FUN=median)
sf(ggplot(D, aes(x,y,colour=cl)) + geom_point(size=0.3, alpha=0.8, stroke=0) +
     scale_colour_manual(values=bigpal) +
     ggplot2::annotate("text", x=cen$x, y=cen$y, label=cen$cl, size=4, fontface="bold") +
     guides(colour=guide_legend(override.aes=list(size=3.5,alpha=1), ncol=1)) +
     labs(title=sprintf("Post-Harmony primary clusters (resolution %s)", primary_res),
          subtitle=sprintf("%d clusters, %d cells - the identities used for marker discovery",
                           nlevels(D$cl), ncol(obj)), x="UMAP 1", y="UMAP 2", colour="cluster") + th,
   "M14_01_primary_cluster_umap", 11, 9, "UMAP_clusters", mp)

# heatmap of top 5 markers per cluster (scaled SCT data on a balanced cell subsample)
top5 <- top_n(sig, 5); hm_genes <- unique(top5$gene)
set.seed(seed)
cells_hm <- unlist(lapply(split(colnames(obj), pcl), function(cc) if (length(cc) > 120) sample(cc,120) else cc))
obj_hm <- subset(obj, cells = cells_hm)
Idents(obj_hm) <- factor(as.character(obj_hm@meta.data[[primary_col]]), levels=levels(obj$.m14_primary))
obj_hm <- ScaleData(obj_hm, features=hm_genes, assay=assay_use, verbose=FALSE)
sf(DoHeatmap(obj_hm, features=hm_genes, assay=assay_use, size=3.2, group.colors=bigpal) +
     scale_fill_viridis_c(option="magma", na.value="white") +
     theme(axis.text.y=element_text(size=5.2)) +
     labs(title=sprintf("Top 5 markers per cluster (resolution %s)", primary_res),
          subtitle="Scaled SCT expression; up to 120 cells sampled per cluster for legibility"),
   "M14_02_marker_heatmap_top5", 16, 20, "heatmap_markers", paste0(mp,";top_n=5;cells_per_cluster<=120"))
rm(obj_hm); invisible(gc(verbose=FALSE))

# dot plot of top 3 markers per cluster
top3 <- top_n(sig, 3); dp_genes <- unique(top3$gene)
Idents(obj) <- obj$.m14_primary
sf(DotPlot(obj, features=dp_genes, assay=assay_use, cluster.idents=FALSE) +
     scale_colour_viridis_c(option="viridis") +
     theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6.2)) +
     labs(title=sprintf("Top 3 markers per cluster (resolution %s)", primary_res),
          subtitle="Dot size = fraction of cells expressing; colour = scaled mean expression",
          x=NULL, y="cluster"),
   "M14_03_marker_dotplot_top3", max(16, 0.16*length(dp_genes)), 9, "dotplot_markers", paste0(mp,";top_n=3"))

# feature plots: the single strongest marker of each cluster
best <- do.call(rbind, lapply(split(sig, sig$cluster), function(d) head(d[order(-d$avg_log2FC),],1)))
fp_genes <- unique(best$gene)
sf(FeaturePlot(obj, features=fp_genes, reduction="postint_umap_harmony", ncol=4,
               order=TRUE, pt.size=0.15, combine=TRUE) &
     scale_colour_viridis_c(option="viridis") & theme_bw(base_size=10) &
     theme(panel.grid=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()),
   "M14_04_featureplots_top_marker_per_cluster",
   14, 3.3*ceiling(length(fp_genes)/4), "featureplot", paste0(mp,";one gene per cluster"))

# violins for the same gene set
sf(VlnPlot(obj, features=head(fp_genes,12), group.by=primary_col, assay=assay_use,
           pt.size=0, ncol=3, combine=TRUE) &
     theme_bw(base_size=9) & theme(legend.position="none", axis.text.x=element_text(size=6)),
   "M14_05_violin_top_markers", 15, 12, "violin_markers", paste0(mp,";first 12 cluster-top genes"))

# markers per cluster barplot
sf(ggplot(summ, aes(factor(cluster, levels=levels(obj$.m14_primary)), n_significant_markers)) +
     geom_col(fill="#3B6EA5") + geom_text(aes(label=n_significant_markers), vjust=-0.3, size=2.9) +
     labs(title="Significant markers per cluster",
          subtitle=sprintf("Wilcoxon, adj. p < 0.05, log2FC >= %.2f, only.pos = TRUE. Cluster characterisation only - NOT condition-level DE.", logfc_thr),
          x="cluster", y="significant markers") + th,
   "M14_06_markers_per_cluster", 12, 5.5, "barplot", mp)
write.table(do.call(rbind, frows), file.path(out_dir,"figure_index_m14.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

rec <- list(milestone="M14", phase="phase2", timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/R/phase2/discover_markers_harmony.R",
  input=list(path=input_rds, md5=in_md5, cells=ncol(obj)),
  assay_layer_inspection=lapply(assay_rows, as.list),
  prep_sct_find_markers=list(required=prep_needed, n_sct_models=n_models, applied=prep_needed),
  marker_parameters=marker_params,
  clustering=list(primary_resolution=primary_res, n_primary_clusters=nlevels(obj$.m14_primary),
                  alternative_resolution=alt_res, n_alternative_clusters=nlevels(obj$.m14_alt)),
  results=list(n_marker_rows_primary=nrow(mk_primary), n_significant_primary=nrow(sig),
               n_marker_rows_alternative=nrow(mk_alt),
               clusters_with_zero_markers=sum(summ$n_significant_markers == 0)),
  statistical_scope="Exploratory cluster-characterisation markers for annotation. NOT condition-level differential expression. Cells are not independent biological replicates.",
  warnings=warns, slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD", intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(), t_start, units="secs")),
  session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir,"m14_marker_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase2_m14_markers", inputs=list(m13_rds=input_rds),
  outputs=list(all_markers=file.path(out_dir,"all_cluster_markers.tsv"),
    filtered=file.path(out_dir,"filtered_cluster_markers.tsv"),
    top=file.path(out_dir,"top_markers_per_cluster.tsv"),
    summary=file.path(out_dir,"marker_summary_by_cluster.tsv"),
    record=file.path(out_dir,"m14_marker_record.json")),
  parameters=c(marker_params, list(milestone="M14", primary_resolution=primary_res,
    slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"))), dataset="combined")
save_provenance_json(prov, file.path(out_dir,"prov_m14_markers.json"))
log_info(sprintf("M14 COMPLETE in %.1f min. Warnings: %d",
                 as.numeric(difftime(Sys.time(),t_start,units="mins")), length(warns)), stage=STAGE)
log_system_usage(stage=STAGE)
