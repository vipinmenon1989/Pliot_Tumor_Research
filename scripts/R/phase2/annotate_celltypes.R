# scripts/R/phase2/annotate_celltypes.R
#
# Phase 2 / Milestone M15 - literature-grounded hierarchical cell-type annotation.
#
# Applies the researcher-auditable annotation map in config/phase2/annotation_map_M15.tsv
# (built from the M14 cluster markers plus the primary MPNST / immunology literature cited
# in that file) to the M13 primary clustering, writes ANNOTATION_EVIDENCE.tsv, and produces
# the annotation figure suite.
#
# The map is data, not code: every label, its marker evidence, its conflicting evidence and
# its citation live in one reviewable TSV. No label is inferred inside this script.
#
# NO CNV inference. Malignant calls are conservative and explicitly labelled as such.

options(stringsAsFactors=FALSE); options(future.globals.maxSize=+Inf)
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix)
  library(ggplot2); library(patchwork); library(viridis); library(RColorBrewer); library(dplyr)
  library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase2_m15_annotation"

args <- commandArgs(trailingOnly=TRUE)
input_rds <- "results/phase2/clustering/phase2_harmony_clustered.rds"
map_tsv   <- "config/phase2/annotation_map_M15.tsv"
out_rds   <- "results/phase2/annotation/phase2_harmony_annotated.rds"
out_dir   <- "results/phase2/annotation"
fig_dir   <- "results/phase2/figures/M15"
seed <- 42L; expected_cells <- 19716L; validation_mode <- FALSE
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--input"){input_rds<-args[i+1];i<-i+2} else if (a=="--map"){map_tsv<-args[i+1];i<-i+2}
  else if (a=="--out-rds"){out_rds<-args[i+1];i<-i+2} else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2}
  else if (a=="--fig-dir"){fig_dir<-args[i+1];i<-i+2}
  else if (a=="--random-seed"){seed<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--expected-cells"){expected_cells<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--validation-mode"){validation_mode<-TRUE;i<-i+1} else stop(sprintf("Unknown argument: %s", a)) }

t_start <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
fail <- function(m){log_error(m,stage=STAGE); stop(m,call.=FALSE)}
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE); dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
set.seed(seed)

log_info("================ M15 CELL-TYPE ANNOTATION ================", stage=STAGE)
in_md5 <- digest(input_rds, file=TRUE, algo="md5")
amap <- read.delim(map_tsv, sep="\t", check.names=FALSE, colClasses="character")
log_info(sprintf("Annotation map: %s (%d clusters, md5 %s)", map_tsv, nrow(amap), digest(map_tsv, file=TRUE, algo="md5")), stage=STAGE)
obj <- readRDS(input_rds)
log_info(sprintf("Loaded %d features x %d cells (md5 %s)", nrow(obj), ncol(obj), in_md5), stage=STAGE)
if (!validation_mode && ncol(obj) != expected_cells) fail("Cell count mismatch.")

pcl <- as.character(obj$postint_harmony_primary_cluster)
clusters <- sort(unique(as.numeric(pcl)))
missing <- setdiff(as.character(clusters), amap$cluster)
if (length(missing)) fail(sprintf("Annotation map lacks clusters: %s", paste(missing,collapse=",")))
extra <- setdiff(amap$cluster, as.character(clusters))
if (length(extra)) wrn(sprintf("Annotation map has entries for absent clusters: %s", paste(extra,collapse=",")))
idx <- match(pcl, amap$cluster)

obj$postint_celltype_level1 <- amap$level1[idx]
obj$postint_celltype_level2 <- amap$level2[idx]
obj$postint_celltype_level3 <- amap$level3[idx]
obj$postint_annotation_confidence <- amap$confidence[idx]
obj$postint_annotation_initial <- amap$level2[idx]      # M15 version; M16 may refine level2
obj$postint_annotation_source_cluster <- pcl
log_info("Annotation columns written: postint_celltype_level1/2/3, postint_annotation_confidence, postint_annotation_initial, postint_annotation_source_cluster", stage=STAGE)
for (lv in c(1,2,3)) {
  tb <- sort(table(obj@meta.data[[sprintf("postint_celltype_level%d", lv)]]), decreasing=TRUE)
  log_info(sprintf("Level %d: %s", lv, paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse=" | ")), stage=STAGE)
}
tbc <- table(obj$postint_annotation_confidence)
log_info(sprintf("Confidence: %s", paste(sprintf("%s=%d", names(tbc), as.integer(tbc)), collapse=" | ")), stage=STAGE)

# ---- ANNOTATION_EVIDENCE.tsv ----------------------------------------------
grp <- as.character(obj$sample_id); lv_s <- sort(unique(grp))
mk <- read.delim("results/phase2/markers/filtered_cluster_markers.tsv", sep="\t")
ev <- do.call(rbind, lapply(amap$cluster, function(k) {
  sel <- pcl == k; s <- mk[mk$cluster == as.integer(k), ]
  s <- s[order(-s$avg_log2FC), ]
  ct <- table(factor(grp[sel], levels=lv_s))
  r <- amap[amap$cluster == k, ]
  data.frame(cluster=k, n_cells=sum(sel), pct_of_all=round(100*sum(sel)/length(sel),3),
    proposed_level1=r$level1, proposed_level2=r$level2, proposed_level3=r$level3,
    positive_markers=r$positive_markers, negative_markers=r$negative_markers,
    supporting_pathways=r$supporting_pathways, conflicting_evidence=r$conflicting_evidence,
    confidence=r$confidence, source=r$source, doi_or_url=r$doi_or_url, notes=r$notes,
    data_derived_top15_markers=paste(head(s$gene,15), collapse=";"),
    n_significant_markers=nrow(s),
    dominant_sample=lv_s[which.max(ct)], dominant_sample_pct=round(100*max(ct)/sum(ct),2),
    n_samples_contributing=sum(ct > 0),
    sample_counts=paste(sprintf("%s=%d", lv_s, as.integer(ct)), collapse=";"),
    stringsAsFactors=FALSE) }))
write.table(ev, file.path(out_dir,"ANNOTATION_EVIDENCE.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info(sprintf("Wrote %s (%d clusters)", file.path(out_dir,"ANNOTATION_EVIDENCE.tsv"), nrow(ev)), stage=STAGE)

cmap <- data.frame(cluster=ev$cluster, n_cells=ev$n_cells, level1=ev$proposed_level1,
  level2=ev$proposed_level2, level3=ev$proposed_level3, confidence=ev$confidence,
  dominant_sample=ev$dominant_sample, dominant_sample_pct=ev$dominant_sample_pct, stringsAsFactors=FALSE)
write.table(cmap, file.path(out_dir,"cluster_to_annotation_map.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# ---- FIGURES ---------------------------------------------------------------
log_info("Generating M15 figures ...", stage=STAGE)
frows <- list()
sf <- function(p,name,w,h,method,params){ for (ext in c("pdf","png")) {
    fp <- file.path(fig_dir, paste0(name,".",ext))
    ggsave(fp,p,width=w,height=h,units="in",dpi=200,device=ext,limitsize=FALSE)
    frows[[length(frows)+1]] <<- data.frame(figure_path=fp, phase="phase2", milestone="M15",
      dataset="combined", processing_stage="celltype_annotation", analysis_method=method,
      parameters=params, input_object_checksum=in_md5,
      generating_script="scripts/R/phase2/annotate_celltypes.R", stringsAsFactors=FALSE) }
  log_info(sprintf("  wrote %s", name), stage=STAGE) }
th <- theme_bw(base_size=13)+theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"),
     plot.subtitle=element_text(size=10,colour="grey30"))
uxy <- Embeddings(obj,"postint_umap_harmony"); set.seed(seed); ordv <- sample.int(ncol(obj))
D <- data.frame(x=uxy[ordv,1], y=uxy[ordv,2],
  L1=obj$postint_celltype_level1[ordv], L2=obj$postint_celltype_level2[ordv],
  L3=obj$postint_celltype_level3[ordv], conf=obj$postint_annotation_confidence[ordv],
  cl=pcl[ordv], sample_id=grp[ordv], stringsAsFactors=FALSE)
D$conf <- factor(D$conf, levels=c("High","Moderate","Low","Uncertain"))

pal_L1 <- c("Malignant / tumour"="#B2182B","T/NK"="#2166AC","Myeloid"="#F4A582","B/Plasma"="#8073AC",
            "Endothelial"="#1B7837","Fibroblast/Stromal"="#D9A441","Other"="#8C8C8C","Uncertain"="#BEBEBE")
pal_conf <- c(High="#1A9850", Moderate="#FEE08B", Low="#F46D43", Uncertain="#9E9E9E")
lab_umap <- function(df,col,title,sub,pal=NULL,labels=TRUE,ncol_leg=1){
  p <- ggplot(df, aes(x,y,colour=.data[[col]]))+geom_point(size=0.32,alpha=0.82,stroke=0)+
    labs(title=title,subtitle=sub,x="UMAP 1",y="UMAP 2",colour=NULL)+th+
    guides(colour=guide_legend(override.aes=list(size=4,alpha=1), ncol=ncol_leg))
  p <- if (!is.null(pal)) p+scale_colour_manual(values=pal, na.value="grey80") else
       p+scale_colour_manual(values=grDevices::colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(length(unique(df[[col]]))))
  if (labels){ cen <- aggregate(cbind(x,y)~get(col), data=df, FUN=median); names(cen)[1]<-"lab"
    p <- p+ggplot2::annotate("text",x=cen$x,y=cen$y,label=cen$lab,size=3.2,fontface="bold") }
  p }

sf(lab_umap(D,"L1","Broad cell-type compartments (Level 1)",
    sprintf("%d cells, %d clusters at Harmony resolution 1.0. Malignant calls are conservative: no CNV evidence was used.", ncol(obj), length(clusters)),
    pal=pal_L1, labels=FALSE),
   "M15_01_harmony_umap_broad_celltypes", 11, 8.5, "UMAP_annotation", "level1;res=1.0")
sf(lab_umap(D,"L2","Detailed cell types (Level 2)",
    "Assigned from M14 cluster markers plus the primary literature cited in ANNOTATION_EVIDENCE.tsv",
    labels=FALSE, ncol_leg=1),
   "M15_02_harmony_umap_detailed_celltypes", 12.5, 8.5, "UMAP_annotation", "level2;res=1.0")
sf(lab_umap(D,"L3","Cell states (Level 3)","States are assigned only where marker evidence supports them", labels=FALSE),
   "M15_03_harmony_umap_celltype_states", 13.5, 8.5, "UMAP_annotation", "level3;res=1.0")
sf(lab_umap(D,"conf","Annotation confidence",
    "High / Moderate / Low / Uncertain. Low and Uncertain regions must not be interpreted as established cell types.",
    pal=pal_conf, labels=FALSE),
   "M15_04_harmony_umap_annotation_confidence", 10.5, 8.5, "UMAP_annotation", "confidence;res=1.0")
sf(lab_umap(D,"cl","M13 primary clusters with annotation overlay",
    "Cluster identifiers used throughout M14/M15; see ANNOTATION_EVIDENCE.tsv for the per-cluster evidence", labels=TRUE)+
     theme(legend.position="none"),
   "M15_05_harmony_umap_clusters_labelled", 10.5, 8.5, "UMAP_clusters", "primary clusters;res=1.0")

# canonical marker panels
panels <- list(
  `MPNST / Schwann-lineage`=c("SOX10","S100B","PLP1","MPZ","NGFR","PMP22","ERBB3","L1CAM","APOD","GFRA3","SHH","CDKN2A"),
  `T / NK`=c("CD3D","CD3E","CD2","IL7R","TRAC","GZMA","GZMK","NKG7","GNLY","KLRD1"),
  `Myeloid`=c("LYZ","CD68","AIF1","CSF1R","C1QA","CD163","FOLR2","IL1B","S100A8","FCN1","CD1C","FCER1A","LILRA4","CLEC4C"),
  `B / Plasma`=c("MS4A1","CD79A","POU2AF1","MZB1","JCHAIN","XBP1","FCRL5"),
  `Fibroblast / stromal`=c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA","PI16","COMP","COL10A1"),
  `Mural`=c("ACTA2","RGS5","PDGFRB","MYH11","NOTCH3"),
  `Endothelial`=c("PECAM1","VWF","CDH5","CLDN5","ACKR1","ESM1","DLL4"),
  `Proliferation`=c("MKI67","TOP2A","UBE2C","CDK1"))
avail <- rownames(obj)
panels <- lapply(panels, function(g) intersect(g, avail))
dot_genes <- unique(unlist(panels))
Idents(obj) <- factor(pcl, levels=as.character(clusters))
lab_for_cluster <- setNames(sprintf("C%s: %s", amap$cluster, amap$level2), amap$cluster)
sf(DotPlot(obj, features=dot_genes, assay="SCT", cluster.idents=FALSE)+
     scale_colour_viridis_c(option="viridis")+
     scale_y_discrete(labels=function(x) lab_for_cluster[x])+
     theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7.5), axis.text.y=element_text(size=8))+
     labs(title="Canonical lineage markers across the M13 primary clusters",
          subtitle="Marker panels for Schwann-lineage/MPNST, T/NK, myeloid, B/plasma, fibroblast, mural, endothelial and proliferation. Sources in ANNOTATION_EVIDENCE.tsv.",
          x=NULL, y=NULL),
   "M15_06_canonical_marker_dotplot", max(15, 0.22*length(dot_genes)), 9.5, "dotplot_canonical", "canonical marker panels")

avg <- AverageExpression(obj, features=dot_genes, assays="SCT", layer="data", group.by="ident")$SCT
zz <- t(scale(t(as.matrix(avg))))
hm <- reshape2::melt(zz, varnames=c("gene","cluster"), value.name="z")
hm$gene <- factor(hm$gene, levels=rev(dot_genes))
hm$cluster <- factor(as.character(hm$cluster), levels=as.character(clusters))
sf(ggplot(hm, aes(cluster, gene, fill=z))+geom_tile()+
     scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0, name="z-score")+
     scale_x_discrete(labels=function(x) lab_for_cluster[x])+
     theme_minimal(base_size=10)+
     theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=8), axis.text.y=element_text(size=7),
           panel.grid=element_blank())+
     labs(title="Canonical marker heatmap (z-scored mean SCT expression per cluster)",
          subtitle="Rows grouped by lineage panel; see ANNOTATION_EVIDENCE.tsv for the assignment evidence", x=NULL, y=NULL),
   "M15_07_canonical_marker_heatmap", 13, 14, "heatmap_canonical", "z-scored AverageExpression")

for (nm in names(panels)) {
  g <- panels[[nm]]; if (!length(g)) next
  safe <- gsub("[^A-Za-z0-9]+","_", nm)
  sf(FeaturePlot(obj, features=g, reduction="postint_umap_harmony", ncol=4, order=TRUE, pt.size=0.14, combine=TRUE) &
       scale_colour_viridis_c(option="viridis") & theme_bw(base_size=9) &
       theme(panel.grid=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()),
     sprintf("M15_08_featureplots_%s", safe), 14, 3.2*ceiling(length(g)/4), "featureplot_panel",
     sprintf("panel=%s;genes=%s", nm, paste(g, collapse=",")))
}

# cluster -> annotation relationship
rel <- as.data.frame(table(cluster=pcl, level2=obj$postint_celltype_level2))
rel <- rel[rel$Freq > 0, ]; rel$cluster <- factor(rel$cluster, levels=as.character(clusters))
sf(ggplot(rel, aes(cluster, level2, size=Freq, colour=Freq))+geom_point()+
     scale_colour_viridis_c(option="rocket", direction=-1, name="cells")+scale_size_area(max_size=9, name="cells")+
     labs(title="How the M13 primary clusters map onto the M15 annotations",
          subtitle="One-to-one by construction: annotation is assigned per cluster, so each cluster contributes to exactly one label",
          x="M13 primary cluster", y="Level 2 cell type")+th+
     theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)),
   "M15_09_cluster_to_annotation_map", 12, 8, "dotplot_mapping", "cluster x level2")

conf_df <- as.data.frame(table(level1=obj$postint_celltype_level1, confidence=obj$postint_annotation_confidence))
conf_df <- conf_df[conf_df$Freq>0,]; conf_df$confidence <- factor(conf_df$confidence, levels=c("High","Moderate","Low","Uncertain"))
sf(ggplot(conf_df, aes(reorder(level1,Freq), Freq, fill=confidence))+geom_col()+coord_flip()+
     scale_fill_manual(values=pal_conf)+
     labs(title="Cells per broad compartment, split by annotation confidence",
          subtitle="Low and Uncertain cells are retained, never forced into a compartment", x=NULL, y="cells")+th,
   "M15_10_annotation_confidence_by_compartment", 10, 6, "barplot", "level1 x confidence")
write.table(do.call(rbind, frows), file.path(out_dir,"figure_index_m15.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# ---- save ------------------------------------------------------------------
log_info(sprintf("Saving annotated object to %s ...", out_rds), stage=STAGE)
Idents(obj) <- factor(obj$postint_celltype_level2)
saveRDS(obj, out_rds)
out_md5 <- digest(out_rds,file=TRUE,algo="md5"); out_sha <- digest(out_rds,file=TRUE,algo="sha256")
out_size <- file.info(out_rds)$size
log_info(sprintf("Saved %.2f GB | md5 %s", out_size/1024^3, out_md5), stage=STAGE)
n_cells <- ncol(obj); rm(obj); invisible(gc(verbose=FALSE))
o2 <- readRDS(out_rds)
v <- list(cells=ncol(o2)==n_cells,
  level1=all(!is.na(o2$postint_celltype_level1)), level2=all(!is.na(o2$postint_celltype_level2)),
  level3=all(!is.na(o2$postint_celltype_level3)), confidence=all(!is.na(o2$postint_annotation_confidence)),
  clusters_preserved=all(sprintf("postint_harmony_clusters_res_%s", format(seq(0.1,1.0,by=0.1),nsmall=1)) %in% colnames(o2@meta.data)),
  reductions=all(c("pca","umap_preintegration","postint_harmony","postint_umap_harmony") %in% Reductions(o2)))
for (nm in names(v)) log_info(sprintf("  %-20s : %s", nm, v[[nm]]), stage=STAGE)
if (!all(unlist(v))) fail(paste("Round-trip validation failed:", paste(names(v)[!unlist(v)], collapse=", ")))
rm(o2); invisible(gc(verbose=FALSE))

rec <- list(milestone="M15", phase="phase2", timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/R/phase2/annotate_celltypes.R", annotation_map=map_tsv,
  annotation_map_md5=digest(map_tsv,file=TRUE,algo="md5"),
  input=list(path=input_rds, md5=in_md5), output=list(path=out_rds, md5=out_md5, sha256=out_sha, size_bytes=out_size),
  annotation_columns=c("postint_celltype_level1","postint_celltype_level2","postint_celltype_level3",
    "postint_annotation_confidence","postint_annotation_initial","postint_annotation_source_cluster"),
  cnv_inference_performed=FALSE,
  malignancy_policy="Malignant calls integrate cluster markers, Schwann/neural-crest lineage biology, MPNST literature, absence of convincing immune/stromal identity and sample provenance. No CNV evidence was used (out of Phase 2 scope). Conservative labels (MPNST-like malignant, candidate malignant, Schwann-lineage tumour-like) are used throughout.",
  round_trip_validation=v, warnings=warns,
  slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(),t_start,units="secs")), session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir,"m15_annotation_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase2_m15_annotation", inputs=list(m13_rds=input_rds, annotation_map=map_tsv),
  outputs=list(annotated_rds=out_rds, evidence=file.path(out_dir,"ANNOTATION_EVIDENCE.tsv"),
    cluster_map=file.path(out_dir,"cluster_to_annotation_map.tsv"), record=file.path(out_dir,"m15_annotation_record.json")),
  parameters=list(milestone="M15", seed=seed, cnv_inference=FALSE,
    slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, file.path(out_dir,"prov_m15_annotation.json"))
log_info(sprintf("M15 COMPLETE in %.1f min. Warnings: %d", as.numeric(difftime(Sys.time(),t_start,units="mins")), length(warns)), stage=STAGE)
log_system_usage(stage=STAGE)
