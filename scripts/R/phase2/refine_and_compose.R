# scripts/R/phase2/refine_and_compose.R
#
# Phase 2 / Milestone M16 - annotation refinement and descriptive composition.
#
# (1) REFINEMENT. Reviews the M15 annotation computationally against canonical compartment
#     marker scores and applies only documented, conservative rules. The M15 labels are
#     PRESERVED in *_initial columns; refinements are written to *_refined columns.
# (2) COMPOSITION. Counts and proportions of each cell type by sample / patient / dataset
#     (one 4-level field) and by condition where one exists (none does here).
#
# DESCRIPTIVE ONLY. No inferential condition-level test is performed. Cells are not
# independent biological replicates.

options(stringsAsFactors=FALSE); options(future.globals.maxSize=+Inf)
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix)
  library(ggplot2); library(patchwork); library(viridis); library(RColorBrewer)
  library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase2_m16_refinement_composition"

args <- commandArgs(trailingOnly=TRUE)
input_rds <- "results/phase2/annotation/phase2_harmony_annotated.rds"
out_rds   <- "results/phase2/composition/phase2_harmony_refined.rds"
out_dir   <- "results/phase2/composition"
fig_dir   <- "results/phase2/figures/M16"
seed <- 42L; expected_cells <- 19716L; validation_mode <- FALSE
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--input"){input_rds<-args[i+1];i<-i+2} else if (a=="--out-rds"){out_rds<-args[i+1];i<-i+2}
  else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2} else if (a=="--fig-dir"){fig_dir<-args[i+1];i<-i+2}
  else if (a=="--random-seed"){seed<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--expected-cells"){expected_cells<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--validation-mode"){validation_mode<-TRUE;i<-i+1} else stop(sprintf("Unknown argument: %s",a)) }

t_start <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
fail <- function(m){log_error(m,stage=STAGE); stop(m,call.=FALSE)}
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE); dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
set.seed(seed)

log_info("========== M16 ANNOTATION REFINEMENT AND COMPOSITION ==========", stage=STAGE)
in_md5 <- digest(input_rds,file=TRUE,algo="md5")
obj <- readRDS(input_rds)
log_info(sprintf("Loaded %d features x %d cells (md5 %s)", nrow(obj), ncol(obj), in_md5), stage=STAGE)
if (!validation_mode && ncol(obj) != expected_cells) fail("Cell count mismatch.")
for (cc in c("postint_celltype_level1","postint_celltype_level2","postint_celltype_level3",
             "postint_annotation_confidence","postint_harmony_primary_cluster"))
  if (!(cc %in% colnames(obj@meta.data))) fail(sprintf("Required M15 column missing: %s", cc))

pcl <- as.character(obj$postint_harmony_primary_cluster)
clusters <- as.character(sort(unique(as.numeric(pcl))))

# ---------------------------------------------------------------------------
# (1) COMPUTATIONAL REVIEW - canonical compartment scores per cluster
# ---------------------------------------------------------------------------
log_info("Scoring canonical compartment panels for the annotation review ...", stage=STAGE)
compartment_panels <- list(
  `Malignant / tumour` = c("SOX10","S100B","PLP1","MPZ","NGFR","PMP22","ERBB3","L1CAM","GFRA3","CRYAB"),
  `T/NK`               = c("CD3D","CD3E","CD3G","CD2","IL7R","TRAC","LCK","NKG7","GNLY","KLRD1","GZMA","GZMK"),
  `Myeloid`            = c("LYZ","CD68","AIF1","CSF1R","ITGAM","C1QA","C1QB","CD163","FCN1","S100A9","CD1C","FCER1A"),
  `B/Plasma`           = c("MS4A1","CD79A","CD79B","POU2AF1","MZB1","JCHAIN","XBP1","FCRL5"),
  `Endothelial`        = c("PECAM1","VWF","CDH5","CLDN5","ACKR1","FLT1","ESM1"),
  `Fibroblast/Stromal` = c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA","ACTA2","RGS5","NOTCH3","PDGFRB"))
avail <- rownames(obj)
panels <- lapply(compartment_panels, function(g) intersect(g, avail))
panels <- panels[vapply(panels, length, 1L) >= 3L]
nbin <- 24L; ctrl <- min(100L, max(5L, as.integer(floor(length(avail)/nbin/2))))
set.seed(seed)
obj <- AddModuleScore(obj, features=panels, name="m16comp_", seed=seed, assay=DefaultAssay(obj),
                      search=FALSE, nbin=nbin, ctrl=ctrl)
sc <- paste0("m16comp_", seq_along(panels))
score_mat <- as.matrix(obj@meta.data[, sc, drop=FALSE]); colnames(score_mat) <- names(panels)
for (cc in sc) obj[[cc]] <- NULL
zc <- scale(score_mat)

cl_scores <- t(vapply(clusters, function(k) colMeans(zc[pcl==k, , drop=FALSE]), numeric(ncol(zc))))
colnames(cl_scores) <- colnames(zc)
assigned_l1 <- vapply(clusters, function(k) unique(obj$postint_celltype_level1[pcl==k])[1], character(1))
review <- do.call(rbind, lapply(seq_along(clusters), function(j) {
  k <- clusters[j]; s <- cl_scores[j, ]; top <- names(s)[which.max(s)]
  a <- assigned_l1[j]
  rank_assigned <- if (a %in% names(s)) rank(-s)[[a]] else NA_real_
  data.frame(cluster=k, assigned_level1=a, top_scoring_compartment=top,
    top_score=as.numeric(max(s)), assigned_score=if (a %in% names(s)) as.numeric(s[[a]]) else NA_real_,
    rank_of_assigned=rank_assigned,
    scores=paste(sprintf("%s=%.2f", names(s), s), collapse=";"),
    stringsAsFactors=FALSE) }))

# Documented, conservative refinement rules. Applied in order; each is logged.
#  R1 Assigned compartment is not scored (Uncertain / Other) -> no change; these are
#     deliberate non-assignments from M15.
#  R2 Assigned compartment is the top-scoring panel                 -> CONFIRMED.
#  R3 Assigned compartment ranks 2nd and is within 0.25 z of the top -> CONFIRMED (weak
#     separation between adjacent compartments is expected).
#  R4 Assigned compartment ranks 3rd or worse, OR is more than 0.5 z below the top
#     -> DOWNGRADE CONFIDENCE one step and record the conflict. The LABEL IS NOT CHANGED:
#     the M15 label rests on specific marker genes and cited literature, which is stronger
#     evidence than a coarse panel mean.
step_down <- c(High="Moderate", Moderate="Low", Low="Low", Uncertain="Uncertain")
review$rule <- NA_character_; review$outcome <- NA_character_
conf_initial <- vapply(clusters, function(k) unique(obj$postint_annotation_confidence[pcl==k])[1], character(1))
review$confidence_initial <- conf_initial; review$confidence_refined <- conf_initial
for (j in seq_len(nrow(review))) {
  a <- review$assigned_level1[j]
  if (!(a %in% colnames(cl_scores))) { review$rule[j] <- "R1"; review$outcome[j] <- "not_scored_deliberate_non_assignment"; next }
  if (review$rank_of_assigned[j] == 1) { review$rule[j] <- "R2"; review$outcome[j] <- "confirmed_top_scoring"; next }
  gap <- review$top_score[j] - review$assigned_score[j]
  if (review$rank_of_assigned[j] == 2 && gap <= 0.25) { review$rule[j] <- "R3"; review$outcome[j] <- "confirmed_close_second"; next }
  review$rule[j] <- "R4"; review$outcome[j] <- "confidence_downgraded_label_retained"
  review$confidence_refined[j] <- step_down[[review$confidence_initial[j]]]
}
write.table(review, file.path(out_dir,"annotation_refinement_review.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info("Annotation review outcomes:", stage=STAGE)
for (j in seq_len(nrow(review)))
  log_info(sprintf("  C%-3s %-22s top=%-20s rank=%s  %s -> %s  [%s]", review$cluster[j],
    substr(review$assigned_level1[j],1,22), substr(review$top_scoring_compartment[j],1,20),
    ifelse(is.na(review$rank_of_assigned[j]),"NA",format(review$rank_of_assigned[j])),
    review$confidence_initial[j], review$confidence_refined[j], review$outcome[j]), stage=STAGE)
n_down <- sum(review$outcome=="confidence_downgraded_label_retained")
log_info(sprintf("Refinement: %d clusters confirmed, %d confidence-downgraded, %d deliberately unscored. NO LABEL WAS CHANGED.",
  sum(review$outcome %in% c("confirmed_top_scoring","confirmed_close_second")), n_down,
  sum(review$outcome=="not_scored_deliberate_non_assignment")), stage=STAGE)

# Preserve M15, write refined
obj$postint_celltype_level1_initial <- obj$postint_celltype_level1
obj$postint_celltype_level2_initial <- obj$postint_celltype_level2
obj$postint_celltype_level3_initial <- obj$postint_celltype_level3
obj$postint_annotation_confidence_initial <- obj$postint_annotation_confidence
idx <- match(pcl, review$cluster)
obj$postint_celltype_level1_refined <- obj$postint_celltype_level1      # labels unchanged by design
obj$postint_celltype_level2_refined <- obj$postint_celltype_level2
obj$postint_celltype_level3_refined <- obj$postint_celltype_level3
obj$postint_annotation_confidence_refined <- review$confidence_refined[idx]
obj$postint_annotation_refined <- obj$postint_celltype_level2
obj$postint_annotation_review_rule <- review$rule[idx]
obj$postint_annotation_review_outcome <- review$outcome[idx]
log_info(sprintf("Confidence after refinement: %s",
  paste(sprintf("%s=%d", names(table(obj$postint_annotation_confidence_refined)),
                as.integer(table(obj$postint_annotation_confidence_refined))), collapse=" | ")), stage=STAGE)

# ---------------------------------------------------------------------------
# (2) COMPOSITION - descriptive only
# ---------------------------------------------------------------------------
log_info("Building composition tables (descriptive; no inferential testing) ...", stage=STAGE)
md <- obj@meta.data
grp <- as.character(md$sample_id); lv <- sort(unique(grp))
L1 <- md$postint_celltype_level1_refined; L2 <- md$postint_celltype_level2_refined
n_cells <- nrow(md)
NOTE <- "DESCRIPTIVE ONLY. Cells are not independent biological replicates; no inferential condition-level or differential-abundance test was performed. sample_id is simultaneously the dataset and the patient. No biological-condition field exists in this dataset."

counts_l1 <- as.data.frame(table(annotation_level1=L1))
counts_l1$pct_of_all <- round(100*counts_l1$Freq/n_cells,3); names(counts_l1)[2] <- "n_cells"
counts_l2 <- as.data.frame(table(annotation_level2=L2))
counts_l2$pct_of_all <- round(100*counts_l2$Freq/n_cells,3); names(counts_l2)[2] <- "n_cells"
counts_l1$note <- NOTE; counts_l2$note <- NOTE
write.table(counts_l1, file.path(out_dir,"cell_counts_by_annotation_level1.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(counts_l2, file.path(out_dir,"cell_counts_by_annotation.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

comp_tab <- function(lab, by, byname) {
  t2 <- table(lab, by); df <- as.data.frame(t2)
  names(df) <- c("annotation", byname, "n_cells")
  tot_by <- colSums(t2); tot_lab <- rowSums(t2)
  df$pct_of_group <- round(100*df$n_cells/tot_by[as.character(df[[byname]])],3)
  df$pct_of_annotation <- round(100*df$n_cells/tot_lab[as.character(df$annotation)],3)
  df$note <- NOTE; df }
for (nm in c("sample","patient","dataset")) {
  write.table(comp_tab(L2, grp, nm), file.path(out_dir, sprintf("cell_proportions_by_%s.tsv", nm)),
              sep="\t", row.names=FALSE, quote=FALSE)
  write.table(comp_tab(L1, grp, nm), file.path(out_dir, sprintf("cell_proportions_level1_by_%s.tsv", nm)),
              sep="\t", row.names=FALSE, quote=FALSE)
}
writeLines(c("cell_proportions_by_condition.tsv was NOT generated.","",
 "No biological-condition field exists anywhere in this dataset. The M10 audit of all 93",
 "metadata columns found no treatment arm, disease stage, anatomical site, primary/metastatic",
 "status, NF1 status, tumour grade, age or sex. Fabricating a condition variable to satisfy a",
 "filename would be scientifically invalid.","",
 "See reports/phase2/HARMONY_VARIABLE_DECISION.md section 2."),
 file.path(out_dir,"cell_proportions_by_condition_NOT_APPLICABLE.txt"))

# is any cell type driven by one sample?
drv <- do.call(rbind, lapply(sort(unique(L2)), function(a) {
  sel <- L2==a; ct <- table(factor(grp[sel], levels=lv))
  data.frame(annotation_level2=a, n_cells=sum(sel),
    dominant_sample=lv[which.max(ct)], dominant_sample_pct=round(100*max(ct)/sum(ct),2),
    n_samples_contributing=sum(ct>0),
    single_sample_driven=(max(ct)/sum(ct)) > 0.80,
    sample_counts=paste(sprintf("%s=%d", lv, as.integer(ct)), collapse=";"),
    note=NOTE, stringsAsFactors=FALSE) }))
write.table(drv, file.path(out_dir,"celltype_sample_dependence.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info(sprintf("Cell types driven >80%% by a single sample/patient: %d of %d",
                 sum(drv$single_sample_driven), nrow(drv)), stage=STAGE)

# ---------------------------------------------------------------------------
# FIGURES
# ---------------------------------------------------------------------------
log_info("Generating M16 figures ...", stage=STAGE)
frows <- list()
sf <- function(p,name,w,h,method,params){ for (ext in c("pdf","png")) {
    fp <- file.path(fig_dir, paste0(name,".",ext))
    ggsave(fp,p,width=w,height=h,units="in",dpi=200,device=ext,limitsize=FALSE)
    frows[[length(frows)+1]] <<- data.frame(figure_path=fp, phase="phase2", milestone="M16",
      dataset="combined", processing_stage="annotation_refinement_composition",
      analysis_method=method, parameters=params, input_object_checksum=in_md5,
      generating_script="scripts/R/phase2/refine_and_compose.R", stringsAsFactors=FALSE) }
  log_info(sprintf("  wrote %s", name), stage=STAGE) }
th <- theme_bw(base_size=13)+theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"),
      plot.subtitle=element_text(size=10,colour="grey30"))
pal_L1 <- c("Malignant / tumour"="#B2182B","T/NK"="#2166AC","Myeloid"="#F4A582","B/Plasma"="#8073AC",
            "Endothelial"="#1B7837","Fibroblast/Stromal"="#D9A441","Other"="#8C8C8C","Uncertain"="#BEBEBE")
pal_conf <- c(High="#1A9850", Moderate="#FEE08B", Low="#F46D43", Uncertain="#9E9E9E")
pal_sample <- setNames(RColorBrewer::brewer.pal(max(3,length(lv)),"Set1")[seq_along(lv)], lv)
uxy <- Embeddings(obj,"postint_umap_harmony"); set.seed(seed); ov <- sample.int(n_cells)
D <- data.frame(x=uxy[ov,1], y=uxy[ov,2], L1=L1[ov], L2=L2[ov],
  conf=factor(md$postint_annotation_confidence_refined[ov], levels=c("High","Moderate","Low","Uncertain")),
  sample_id=grp[ov], stringsAsFactors=FALSE)
umapf <- function(col,title,sub,pal=NULL,ncl=1){
  p <- ggplot(D,aes(x,y,colour=.data[[col]]))+geom_point(size=0.32,alpha=0.82,stroke=0)+
    labs(title=title,subtitle=sub,x="UMAP 1",y="UMAP 2",colour=NULL)+th+
    guides(colour=guide_legend(override.aes=list(size=4,alpha=1),ncol=ncl))
  if (!is.null(pal)) p+scale_colour_manual(values=pal,na.value="grey80") else
    p+scale_colour_manual(values=grDevices::colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(length(unique(D[[col]])))) }
sf(umapf("L1","Final annotated UMAP - broad compartments (Level 1, refined)",
   sprintf("%d cells. Malignant calls are conservative and CNV-free; Uncertain cells are retained, not forced.", n_cells), pal_L1),
   "M16_01_final_umap_broad_celltypes", 11, 8.5, "UMAP_annotation", "level1_refined")
sf(umapf("L2","Final annotated UMAP - detailed cell types (Level 2, refined)",
   "Labels unchanged from M15; the M16 review adjusted confidence only", ncl=1),
   "M16_02_final_umap_detailed_celltypes", 12.5, 8.5, "UMAP_annotation", "level2_refined")
sf(umapf("conf","Annotation confidence after the M16 review","High / Moderate / Low / Uncertain", pal_conf),
   "M16_03_final_umap_annotation_confidence", 10.5, 8.5, "UMAP_annotation", "confidence_refined")

cl1 <- comp_tab(L1, grp, "sample"); cl2 <- comp_tab(L2, grp, "sample")
sf(ggplot(cl1, aes(sample, pct_of_group, fill=annotation))+geom_col()+scale_fill_manual(values=pal_L1)+
   labs(title="Broad compartment composition by sample (= patient = dataset)",
        subtitle=paste("DESCRIPTIVE ONLY - no differential-abundance test was performed.", "Four samples, one per patient; no condition variable exists."),
        x=NULL,y="% of sample")+th,
   "M16_04_composition_broad_by_sample", 9, 6.5, "stacked_barplot", "level1 x sample")
sf(ggplot(cl2, aes(sample, pct_of_group, fill=annotation))+geom_col()+
   scale_fill_manual(values=grDevices::colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(length(unique(cl2$annotation))))+
   labs(title="Detailed cell-type composition by sample (= patient = dataset)",
        subtitle="DESCRIPTIVE ONLY - no differential-abundance test was performed.", x=NULL,y="% of sample")+th,
   "M16_05_composition_detailed_by_sample", 11, 7.5, "stacked_barplot", "level2 x sample")
sf(ggplot(cl1, aes(sample, n_cells, fill=annotation))+geom_col()+scale_fill_manual(values=pal_L1)+
   labs(title="Broad compartment counts by sample", subtitle="Absolute cell numbers", x=NULL,y="cells")+th,
   "M16_06_composition_counts_by_sample", 9, 6.5, "stacked_barplot", "level1 counts x sample")
co <- counts_l2[order(-counts_l2$n_cells),]; co$annotation_level2 <- factor(co$annotation_level2, levels=co$annotation_level2)
sf(ggplot(co, aes(annotation_level2, n_cells))+geom_col(fill="#3B6EA5")+
   geom_text(aes(label=sprintf("%d (%.1f%%)", n_cells, pct_of_all)), hjust=-0.05, size=3.1)+
   coord_flip(clip="off")+expand_limits(y=max(co$n_cells)*1.25)+
   labs(title="Cell-type abundance overview (Level 2)", subtitle=sprintf("%d cells total", n_cells), x=NULL,y="cells")+th,
   "M16_07_celltype_abundance_overview", 10, 8, "barplot", "level2 counts")
sf(ggplot(cl2, aes(annotation, pct_of_annotation, fill=sample))+geom_col()+coord_flip()+
   scale_fill_manual(values=pal_sample)+geom_hline(yintercept=80, linetype="dashed", colour="grey20")+
   labs(title="Sample / patient contribution to each cell type",
        subtitle="Dashed line at 80% marks cell types driven almost entirely by one patient - these must not be read as general biology",
        x=NULL, y="% of that cell type")+th,
   "M16_08_sample_contribution_per_celltype", 11, 8, "stacked_barplot", "level2 x sample share")
sf(ggplot(cl1, aes(annotation, pct_of_annotation, fill=sample))+geom_col()+coord_flip()+
   scale_fill_manual(values=pal_sample)+geom_hline(yintercept=80, linetype="dashed", colour="grey20")+
   labs(title="Patient contribution to each broad compartment", subtitle="Same view at Level 1", x=NULL, y="% of that compartment")+th,
   "M16_09_patient_contribution_per_compartment", 10, 6, "stacked_barplot", "level1 x sample share")
cdf <- as.data.frame(table(level1=L1, confidence=md$postint_annotation_confidence_refined))
cdf <- cdf[cdf$Freq>0,]; cdf$confidence <- factor(cdf$confidence, levels=c("High","Moderate","Low","Uncertain"))
sf(ggplot(cdf, aes(reorder(level1,Freq), Freq, fill=confidence))+geom_col()+coord_flip()+
   scale_fill_manual(values=pal_conf)+
   labs(title="Annotation confidence by compartment after the M16 review", x=NULL, y="cells")+th,
   "M16_10_annotation_confidence_after_review", 10, 6, "barplot", "level1 x confidence_refined")
write.table(do.call(rbind, frows), file.path(out_dir,"figure_index_m16.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

log_info(sprintf("Saving refined object to %s ...", out_rds), stage=STAGE)
saveRDS(obj, out_rds)
out_md5 <- digest(out_rds,file=TRUE,algo="md5"); out_sha <- digest(out_rds,file=TRUE,algo="sha256")
log_info(sprintf("Saved %.2f GB | md5 %s", file.info(out_rds)$size/1024^3, out_md5), stage=STAGE)
rm(obj); invisible(gc(verbose=FALSE))
o2 <- readRDS(out_rds)
v <- list(cells=ncol(o2)==n_cells,
  initial_preserved=all(c("postint_celltype_level1_initial","postint_celltype_level2_initial",
                          "postint_annotation_confidence_initial") %in% colnames(o2@meta.data)),
  refined_present=all(c("postint_celltype_level1_refined","postint_celltype_level2_refined",
                        "postint_annotation_confidence_refined") %in% colnames(o2@meta.data)),
  no_missing_labels=all(!is.na(o2$postint_celltype_level2_refined)),
  clusters_preserved=all(sprintf("postint_harmony_clusters_res_%s", format(seq(0.1,1.0,by=0.1),nsmall=1)) %in% colnames(o2@meta.data)),
  reductions=all(c("pca","umap_preintegration","postint_harmony","postint_umap_harmony") %in% Reductions(o2)))
for (nm in names(v)) log_info(sprintf("  %-20s : %s", nm, v[[nm]]), stage=STAGE)
if (!all(unlist(v))) fail(paste("Round-trip validation failed:", paste(names(v)[!unlist(v)],collapse=", ")))
rm(o2); invisible(gc(verbose=FALSE))

rec <- list(milestone="M16", phase="phase2", timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/R/phase2/refine_and_compose.R", input=list(path=input_rds, md5=in_md5),
  output=list(path=out_rds, md5=out_md5, sha256=out_sha, size_bytes=file.info(out_rds)$size),
  refinement=list(rules="R1 unscored/deliberate; R2 assigned = top-scoring panel -> confirmed; R3 assigned ranks 2nd within 0.25 z -> confirmed; R4 otherwise -> confidence downgraded one step, LABEL RETAINED",
    labels_changed=0, confidence_downgraded=n_down,
    rationale="M15 labels rest on specific marker genes and cited literature, which is stronger evidence than a coarse compartment panel mean; the review therefore adjusts confidence rather than overriding labels."),
  composition=list(statistical_scope="DESCRIPTIVE ONLY - no inferential condition-level or differential-abundance test. Cells are not independent biological replicates.",
    condition_field_exists=FALSE,
    celltypes_single_sample_driven=sum(drv$single_sample_driven)),
  preserved_columns=c("postint_celltype_level1_initial","postint_celltype_level2_initial",
    "postint_celltype_level3_initial","postint_annotation_confidence_initial"),
  round_trip_validation=v, warnings=warns,
  slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(),t_start,units="secs")), session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir,"m16_composition_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase2_m16_refinement_composition", inputs=list(m15_rds=input_rds),
  outputs=list(refined_rds=out_rds, review=file.path(out_dir,"annotation_refinement_review.tsv"),
    counts=file.path(out_dir,"cell_counts_by_annotation.tsv"),
    record=file.path(out_dir,"m16_composition_record.json")),
  parameters=list(milestone="M16", seed=seed, inferential_testing=FALSE,
    slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, file.path(out_dir,"prov_m16_composition.json"))
log_info(sprintf("M16 COMPLETE in %.1f min. Warnings: %d", as.numeric(difftime(Sys.time(),t_start,units="mins")), length(warns)), stage=STAGE)
log_system_usage(stage=STAGE)
