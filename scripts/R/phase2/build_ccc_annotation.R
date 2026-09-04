# scripts/R/phase2/build_ccc_annotation.R
#
# Phase 2 / Milestone M15A - CCC-oriented annotation layer (approved amendment, 2026-09-03).
#
# Builds `annotation_ccc` from the DETAILED M15/M16 annotation, collapsing only the
# evidence-supported malignant MPNST populations into a single `MPNST-Tumor` identity while
# immune, stromal and endothelial populations retain biologically meaningful identities.
#
# HARD RULES ENFORCED HERE
#   * The detailed annotation is the SOURCE OF TRUTH and is never modified or deleted.
#   * "non-immune => tumour" is NOT used. The collapse is driven entirely by the auditable
#     map in config/phase2/ccc_annotation_map.tsv, keyed on the detailed Level 2 label.
#   * Ambiguous, uncertain and low-quality populations are never absorbed into MPNST-Tumor.
#   * Every cluster gets an entry in CCC_ANNOTATION_MAPPING.tsv, including every exclusion.
#   * Marker evidence comes from the SCT expression assay, never from Harmony coordinates.
#
# Harmony, M12 metrics, M13 clustering and M14 marker discovery are NOT rerun.

options(stringsAsFactors=FALSE); options(future.globals.maxSize=+Inf)
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix)
  library(ggplot2); library(patchwork); library(viridis); library(RColorBrewer); library(reshape2)
  library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase2_m15a_ccc_annotation"

args <- commandArgs(trailingOnly=TRUE)
input_rds <- "results/phase2/composition/phase2_harmony_refined.rds"
map_tsv   <- "config/phase2/ccc_annotation_map.tsv"
out_rds   <- "results/phase2/annotation/phase2_harmony_ccc.rds"
out_dir   <- "results/phase2/annotation"
comp_dir  <- "results/phase2/composition"
fig_dir   <- "results/phase2/figures/CCC_annotation"
seed <- 42L; expected_cells <- 19716L; validation_mode <- FALSE
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--input"){input_rds<-args[i+1];i<-i+2} else if (a=="--map"){map_tsv<-args[i+1];i<-i+2}
  else if (a=="--out-rds"){out_rds<-args[i+1];i<-i+2} else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2}
  else if (a=="--comp-dir"){comp_dir<-args[i+1];i<-i+2} else if (a=="--fig-dir"){fig_dir<-args[i+1];i<-i+2}
  else if (a=="--random-seed"){seed<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--expected-cells"){expected_cells<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--validation-mode"){validation_mode<-TRUE;i<-i+1} else stop(sprintf("Unknown argument: %s",a)) }

t_start <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
fail <- function(m){log_error(m,stage=STAGE); stop(m,call.=FALSE)}
for (d in c(out_dir, comp_dir, fig_dir)) dir.create(d, recursive=TRUE, showWarnings=FALSE)
set.seed(seed)

log_info("=========== M15A CCC-ORIENTED ANNOTATION LAYER ===========", stage=STAGE)
in_md5 <- digest(input_rds, file=TRUE, algo="md5")
map_md5 <- digest(map_tsv, file=TRUE, algo="md5")
cmap <- read.delim(map_tsv, sep="\t", colClasses="character", check.names=FALSE)
log_info(sprintf("Input %s (md5 %s)", input_rds, in_md5), stage=STAGE)
log_info(sprintf("CCC map %s (md5 %s, %d rules)", map_tsv, map_md5, nrow(cmap)), stage=STAGE)
obj <- readRDS(input_rds)
log_info(sprintf("Loaded %d features x %d cells", nrow(obj), ncol(obj)), stage=STAGE)
if (!validation_mode && ncol(obj) != expected_cells) fail("Cell count mismatch.")

DET_L1 <- "postint_celltype_level1_refined"; DET_L2 <- "postint_celltype_level2_refined"
DET_L3 <- "postint_celltype_level3_refined"; DET_CONF <- "postint_annotation_confidence_refined"
for (cc in c(DET_L1,DET_L2,DET_L3,DET_CONF,"postint_harmony_primary_cluster","sample_id"))
  if (!(cc %in% colnames(obj@meta.data))) fail(sprintf("Required column missing: %s", cc))

md <- obj@meta.data
pcl <- as.character(md$postint_harmony_primary_cluster)
clusters <- as.character(sort(unique(as.numeric(pcl))))
det_l2 <- as.character(md[[DET_L2]])
grp <- as.character(md$sample_id); lv <- sort(unique(grp)); n_cells <- nrow(md)

# preserve the detailed annotation exactly as it was (digest proof after)
det_digests <- sapply(c(DET_L1,DET_L2,DET_L3,DET_CONF,
  "postint_celltype_level1","postint_celltype_level2","postint_celltype_level3",
  "postint_annotation_confidence","postint_celltype_level1_initial",
  "postint_celltype_level2_initial","postint_annotation_confidence_initial"),
  function(c) if (c %in% colnames(md)) digest(md[[c]], algo="md5") else NA_character_)

# ---------------------------------------------------------------------------
# 1. Apply the map (keyed on the DETAILED Level 2 label - never on "not immune")
# ---------------------------------------------------------------------------
missing_rules <- setdiff(unique(det_l2), cmap$level2_detailed)
if (length(missing_rules)) fail(sprintf("CCC map has no rule for detailed labels: %s", paste(missing_rules, collapse=" | ")))
idx <- match(det_l2, cmap$level2_detailed)
ccc <- cmap$annotation_ccc[idx]
ccc_comp <- cmap$ccc_compartment[idx]
ccc_tumor <- cmap$collapsed_to_mpnst_tumor[idx] == "TRUE"
ccc_ready <- cmap$ccc_ready[idx] == "TRUE"

# ---------------------------------------------------------------------------
# 2. Sub-resolve the single T/NK parent cluster by canonical marker gating.
#    Restricted to that cluster only. Not applied anywhere else, and no
#    unsupported subtype is invented: cells that do not resolve stay T-cell-other.
# ---------------------------------------------------------------------------
tnk_cells <- which(ccc == "T-NK-parent")
tnk_counts <- c(`CD8-T`=0L, `CD4-T`=0L, NK=0L, `T-cell-other`=0L)
if (length(tnk_cells)) {
  log_info(sprintf("Sub-resolving the T/NK parent cluster (%d cells) by canonical marker gating ...", length(tnk_cells)), stage=STAGE)
  gg <- function(g) { g <- intersect(g, rownames(obj))
    if (!length(g)) return(rep(0, length(tnk_cells)))
    m <- GetAssayData(obj, assay="SCT", layer="data")[g, colnames(obj)[tnk_cells], drop=FALSE]
    Matrix::colSums(m) }
  cd3 <- gg(c("CD3D","CD3E","CD3G")); cd8 <- gg(c("CD8A","CD8B"))
  cd4 <- gg(c("CD4","IL7R","CCR7")); nkg <- gg(c("NKG7","GNLY","KLRD1","NCAM1","KLRF1"))
  sub <- rep("T-cell-other", length(tnk_cells))
  sub[cd3 == 0 & nkg > 0]            <- "NK"
  sub[cd3 > 0 & cd8 > 0]             <- "CD8-T"
  sub[cd3 > 0 & cd8 == 0 & cd4 > 0]  <- "CD4-T"
  ccc[tnk_cells] <- sub
  ccc_comp[tnk_cells] <- "Immune"; ccc_tumor[tnk_cells] <- FALSE; ccc_ready[tnk_cells] <- TRUE
  tnk_counts <- table(factor(sub, levels=names(tnk_counts)))
  log_info(sprintf("  gating: %s", paste(sprintf("%s=%d", names(tnk_counts), as.integer(tnk_counts)), collapse=", ")), stage=STAGE)
  log_info("  gate logic: NK = CD3(D/E/G)==0 & (NKG7|GNLY|KLRD1|NCAM1|KLRF1)>0; CD8-T = CD3>0 & (CD8A|CD8B)>0; CD4-T = CD3>0 & CD8==0 & (CD4|IL7R|CCR7)>0; else T-cell-other.", stage=STAGE)
  wrn("CD4 transcript detection is sparse in droplet scRNA-seq, so CD4-T relies partly on IL7R/CCR7 and its boundary with T-cell-other is soft. Reported but flagged.")
}
if (any(ccc == "T-NK-parent")) fail("T/NK parent label survived sub-resolution.")

obj$annotation_ccc <- ccc
obj$annotation_ccc_compartment <- ccc_comp
obj$annotation_ccc_is_tumor <- ccc_tumor
obj$annotation_ccc_ccc_ready <- ccc_ready

tb <- sort(table(ccc), decreasing=TRUE)
log_info(sprintf("annotation_ccc (%d identities): %s", length(tb),
  paste(sprintf("%s=%d(%.1f%%)", names(tb), as.integer(tb), 100*as.integer(tb)/n_cells), collapse=" | ")), stage=STAGE)
log_info(sprintf("MPNST-Tumor: %d cells (%.2f%%) | compartments: %s", sum(ccc=="MPNST-Tumor"),
  100*sum(ccc=="MPNST-Tumor")/n_cells,
  paste(sprintf("%s=%d", names(table(ccc_comp)), as.integer(table(ccc_comp))), collapse=", ")), stage=STAGE)

# detailed annotation must be untouched
for (cc in names(det_digests)) if (!is.na(det_digests[[cc]]) &&
  !identical(digest(obj@meta.data[[cc]], algo="md5"), det_digests[[cc]]))
  fail(sprintf("Detailed annotation column '%s' was modified. Aborting.", cc))
log_info("PASS - every detailed annotation column is byte-identical.", stage=STAGE)

# ---------------------------------------------------------------------------
# 3. CCC_ANNOTATION_MAPPING.tsv - one row per cluster, every decision explained
# ---------------------------------------------------------------------------
ev_path <- file.path(out_dir,"ANNOTATION_EVIDENCE.tsv")
ev <- read.delim(ev_path, sep="\t", colClasses="character", check.names=FALSE)
mapping <- do.call(rbind, lapply(clusters, function(k) {
  sel <- pcl == k; e <- ev[ev$cluster == k, ]; r <- cmap[cmap$level2_detailed == unique(det_l2[sel])[1], ]
  cc_here <- sort(table(ccc[sel]), decreasing=TRUE)
  ct <- table(factor(grp[sel], levels=lv))
  data.frame(original_cluster=k, n_cells=sum(sel),
    annotation_level1=unique(md[[DET_L1]][sel])[1], annotation_level2=unique(det_l2[sel])[1],
    annotation_level3=unique(md[[DET_L3]][sel])[1], annotation_confidence=unique(md[[DET_CONF]][sel])[1],
    annotation_ccc=paste(sprintf("%s(%d)", names(cc_here), as.integer(cc_here)), collapse=";"),
    ccc_compartment=unique(ccc_comp[sel])[1],
    collapsed_to_mpnst_tumor=toupper(as.character(all(ccc_tumor[sel]))),
    ccc_ready=toupper(as.character(all(ccc_ready[sel]))),
    positive_evidence=if (nrow(e)) e$positive_markers[1] else NA_character_,
    conflicting_evidence=if (nrow(e)) e$conflicting_evidence[1] else NA_character_,
    mapping_rationale=r$mapping_rationale,
    source=if (nrow(e)) e$source[1] else NA_character_,
    doi_or_url=if (nrow(e)) e$doi_or_url[1] else NA_character_,
    data_derived_top15_markers=if (nrow(e)) e$data_derived_top15_markers[1] else NA_character_,
    dominant_sample=lv[which.max(ct)], dominant_sample_pct=round(100*max(ct)/sum(ct),2),
    n_samples_contributing=sum(ct>0),
    notes=if (nrow(e)) e$notes[1] else NA_character_, stringsAsFactors=FALSE) }))
write.table(mapping, file.path(out_dir,"CCC_ANNOTATION_MAPPING.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info(sprintf("Clusters collapsed into MPNST-Tumor: %s", paste(mapping$original_cluster[mapping$collapsed_to_mpnst_tumor=="TRUE"], collapse=", ")), stage=STAGE)
log_info(sprintf("Clusters explicitly NOT collapsed: %s", paste(mapping$original_cluster[mapping$collapsed_to_mpnst_tumor=="FALSE"], collapse=", ")), stage=STAGE)

# augment ANNOTATION_EVIDENCE.tsv (append columns; nothing deleted)
ev$annotation_ccc <- vapply(ev$cluster, function(k){ s <- pcl==k
  n <- sort(table(ccc[s]), decreasing=TRUE); paste(names(n), collapse=";") }, character(1))
ev$collapsed_to_mpnst_tumor <- vapply(ev$cluster, function(k) toupper(as.character(all(ccc_tumor[pcl==k]))), character(1))
ev$ccc_compartment <- vapply(ev$cluster, function(k) unique(ccc_comp[pcl==k])[1], character(1))
ev$ccc_mapping_rationale <- mapping$mapping_rationale[match(ev$cluster, mapping$original_cluster)]
write.table(ev, ev_path, sep="\t", row.names=FALSE, quote=FALSE)
log_info("ANNOTATION_EVIDENCE.tsv augmented with the CCC traceability columns (nothing deleted).", stage=STAGE)

# ---------------------------------------------------------------------------
# 4. Summary / size audit / tumour composition / representation
# ---------------------------------------------------------------------------
NOTE <- "DESCRIPTIVE ONLY. Cells are not independent biological replicates; no inferential differential-abundance or condition-level test was performed. sample_id is simultaneously the dataset and the patient. No biological-condition field exists in this dataset."
ids <- names(sort(table(ccc), decreasing=TRUE))
summ <- do.call(rbind, lapply(ids, function(a){ s <- ccc==a; ct <- table(factor(grp[s], levels=lv))
  data.frame(annotation_ccc=a, n_cells=sum(s), percent_cells=round(100*sum(s)/n_cells,3),
    n_samples=sum(ct>0), n_patients=sum(ct>0),
    min_cells_per_sample=min(ct), max_cells_per_sample=max(ct),
    note=NOTE, stringsAsFactors=FALSE) }))
write.table(summ, file.path(out_dir,"CCC_ANNOTATION_SUMMARY.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

MIN_TOTAL <- 100L; MIN_PER_SAMPLE <- 10L; MIN_SAMPLES <- 2L
audit <- do.call(rbind, lapply(ids, function(a){ s <- ccc==a; ct <- table(factor(grp[s], levels=lv))
  ready_flag <- all(ccc_ready[s])
  ok_total <- sum(s) >= MIN_TOTAL; ok_samples <- sum(ct>0) >= MIN_SAMPLES
  ok_per <- sum(ct >= MIN_PER_SAMPLE) >= MIN_SAMPLES
  status <- if (!ready_flag) "EXCLUDE - flagged not CCC-ready by the annotation map" else
            if (ok_total && ok_samples && ok_per) "READY" else
            if (ok_total && ok_samples) "READY_WITH_CAUTION - thin in some samples" else
            if (ok_total) "CAUTION - present in too few samples for sample-aware inference" else
            "TOO_RARE - below 100 cells total"
  data.frame(annotation_ccc=a, total_cells=sum(s), samples_present=sum(ct>0), patients_present=sum(ct>0),
    median_cells_per_sample=as.numeric(median(ct[ct>0])),
    minimum_cells_per_sample=as.integer(min(ct[ct>0])), maximum_cells_per_sample=as.integer(max(ct)),
    samples_with_ge10_cells=sum(ct>=MIN_PER_SAMPLE),
    ccc_readiness=status,
    notes=sprintf("thresholds: >=%d cells total, present in >=%d samples with >=%d cells each; per-sample counts %s",
                  MIN_TOTAL, MIN_SAMPLES, MIN_PER_SAMPLE, paste(sprintf("%s=%d", lv, as.integer(ct)), collapse=";")),
    stringsAsFactors=FALSE) }))
write.table(audit, file.path(out_dir,"CCC_POPULATION_SIZE_AUDIT.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info("CCC population size audit:", stage=STAGE)
for (r in seq_len(nrow(audit)))
  log_info(sprintf("  %-30s n=%5d  samples=%d  min/sample=%4d  %s", audit$annotation_ccc[r],
    audit$total_cells[r], audit$samples_present[r], audit$minimum_cells_per_sample[r], audit$ccc_readiness[r]), stage=STAGE)

tsel <- ccc == "MPNST-Tumor"
tumor_comp <- do.call(rbind, lapply(sort(unique(det_l2[tsel])), function(a){ s <- tsel & det_l2==a
  ct <- table(factor(grp[s], levels=lv))
  data.frame(detailed_tumor_annotation=a, cells=sum(s),
    percent_of_mpnst_tumor=round(100*sum(s)/sum(tsel),2),
    clusters=paste(sort(unique(pcl[s])), collapse=";"),
    samples=paste(lv[ct>0], collapse=";"), patients=paste(lv[ct>0], collapse=";"),
    confidence=paste(sort(unique(md[[DET_CONF]][s])), collapse=";"),
    sample_counts=paste(sprintf("%s=%d", lv, as.integer(ct)), collapse=";"),
    stringsAsFactors=FALSE) }))
write.table(tumor_comp, file.path(out_dir,"MPNST_TUMOR_COMPOSITION.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info("MPNST-Tumor internal composition:", stage=STAGE)
for (r in seq_len(nrow(tumor_comp)))
  log_info(sprintf("  %-34s %5d cells (%5.2f%% of MPNST-Tumor)  clusters %s  conf %s",
    tumor_comp$detailed_tumor_annotation[r], tumor_comp$cells[r], tumor_comp$percent_of_mpnst_tumor[r],
    tumor_comp$clusters[r], tumor_comp$confidence[r]), stage=STAGE)

ctab <- table(ccc, grp)
cnt <- as.data.frame(ctab); names(cnt) <- c("annotation_ccc","sample","n_cells"); cnt$note <- NOTE
prp <- cnt; prp$pct_of_sample <- round(100*prp$n_cells/colSums(ctab)[as.character(prp$sample)],3)
prp$pct_of_ccc_population <- round(100*prp$n_cells/rowSums(ctab)[as.character(prp$annotation_ccc)],3)
for (nm in c("sample","patient")) {
  c2 <- cnt; p2 <- prp; names(c2)[2] <- nm; names(p2)[2] <- nm
  write.table(c2, file.path(comp_dir, sprintf("ccc_counts_by_%s.tsv", nm)), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(p2, file.path(comp_dir, sprintf("ccc_proportions_by_%s.tsv", nm)), sep="\t", row.names=FALSE, quote=FALSE)
}
writeLines(c("ccc_counts_by_condition.tsv and ccc_proportions_by_condition.tsv were NOT generated.","",
 "No biological-condition field exists anywhere in this dataset. The M10 audit of all 93 metadata",
 "columns found no treatment arm, disease stage, anatomical site, primary/metastatic status, NF1",
 "status, tumour grade, age or sex. Fabricating a condition variable to satisfy a filename would be",
 "scientifically invalid. sample_id is simultaneously the dataset and the patient, so the by-sample",
 "and by-patient tables already cover every grouping the data supports.","",
 "See reports/phase2/HARMONY_VARIABLE_DECISION.md section 2."),
 file.path(comp_dir,"ccc_by_condition_NOT_APPLICABLE.txt"))

# ---------------------------------------------------------------------------
# 5. CCC-level marker summary. Reuses the M14 marker tables (no rediscovery) plus
#    mean SCT expression of canonical genes per CCC identity.
# ---------------------------------------------------------------------------
log_info("Summarising representative markers per CCC identity (reusing M14 tables; no marker rediscovery) ...", stage=STAGE)
mk <- read.delim("results/phase2/markers/filtered_cluster_markers.tsv", sep="\t")
ccc_by_cluster <- setNames(mapping$ccc_compartment, mapping$original_cluster)
ccc_lab_cluster <- setNames(vapply(mapping$original_cluster, function(k) names(sort(table(ccc[pcl==k]), decreasing=TRUE))[1], character(1)),
                            mapping$original_cluster)
mk$annotation_ccc <- ccc_lab_cluster[as.character(mk$cluster)]
msum <- do.call(rbind, lapply(sort(unique(mk$annotation_ccc)), function(a){
  d <- mk[mk$annotation_ccc==a, ]; d <- d[order(-d$avg_log2FC), ]
  g <- unique(d$gene)
  data.frame(annotation_ccc=a, n_contributing_clusters=length(unique(d$cluster)),
    contributing_clusters=paste(sort(unique(d$cluster)), collapse=";"),
    n_significant_markers=nrow(d), top20_representative_markers=paste(head(g,20), collapse=";"),
    max_avg_log2FC=round(max(d$avg_log2FC),3),
    scope="Cluster-characterisation markers reused from M14 (SCT assay, PrepSCTFindMarkers applied). NOT condition-level differential expression.",
    stringsAsFactors=FALSE) }))
write.table(msum, file.path(out_dir,"CCC_MARKER_SUMMARY.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# ---------------------------------------------------------------------------
# 6. FIGURES (15)
# ---------------------------------------------------------------------------
log_info("Generating the CCC figure suite ...", stage=STAGE)
frows <- list()
sf <- function(p,name,w,h,ftype,afield,params,notes=""){ for (ext in c("pdf","png")) {
    fp <- file.path(fig_dir, paste0(name,".",ext))
    ggsave(fp,p,width=w,height=h,units="in",dpi=200,device=ext,limitsize=FALSE)
    frows[[length(frows)+1]] <<- data.frame(figure_path=fp, phase="phase2", milestone="M15A",
      figure_type=ftype, annotation_field=afield, input_object=input_rds,
      input_checksum=in_md5, generating_script="scripts/R/phase2/build_ccc_annotation.R",
      parameters=params, notes=notes, stringsAsFactors=FALSE) }
  log_info(sprintf("  wrote %s", name), stage=STAGE) }
th <- theme_bw(base_size=13)+theme(panel.grid.minor=element_blank(),
      plot.title=element_text(face="bold"), plot.subtitle=element_text(size=9.5,colour="grey30"))
uxy <- Embeddings(obj,"postint_umap_harmony"); set.seed(seed); ov <- sample.int(n_cells)
D <- data.frame(x=uxy[ov,1], y=uxy[ov,2], cl=pcl[ov], det=det_l2[ov], ccc=ccc[ov],
  comp=ccc_comp[ov], conf=md[[DET_CONF]][ov], sample_id=grp[ov], stringsAsFactors=FALSE)
D$conf <- factor(D$conf, levels=c("High","Moderate","Low","Uncertain"))
bigpal <- function(n) grDevices::colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(n)
pal_ccc <- c(`MPNST-Tumor`="#B2182B", `Candidate-Malignant-Unresolved`="#E8A0A0",
  `CD8-T`="#08519C", `CD4-T`="#4292C6", NK="#9ECAE1", `T-cell-other`="#C6DBEF",
  Macrophage="#FD8D3C", Monocyte="#FDBE85", Dendritic="#E6550D", `Plasmacytoid-DC`="#A63603",
  `B-cell`="#807DBA", `Plasma-cell`="#54278F", Fibroblast="#D9A441", `Pericyte-VSMC"`="#8C6D1F",
  `Pericyte-VSMC`="#8C6D1F", Endothelial="#1B7837", Uncertain="#BEBEBE", `Low-quality-excluded`="#6E6E6E")
pal_comp <- c(`MPNST-Tumor`="#B2182B", Immune="#2166AC", `Fibroblast-Stromal`="#D9A441",
  Endothelial="#1B7837", Other="#8C8C8C", Uncertain="#BEBEBE")
pal_conf <- c(High="#1A9850", Moderate="#FEE08B", Low="#F46D43", Uncertain="#9E9E9E")
um <- function(col,title,sub,pal=NULL,lab=FALSE,ncl=1,ptsize=0.32){
  p <- ggplot(D,aes(x,y,colour=.data[[col]]))+geom_point(size=ptsize,alpha=0.82,stroke=0)+
    labs(title=title,subtitle=sub,x="UMAP 1",y="UMAP 2",colour=NULL)+th+
    guides(colour=guide_legend(override.aes=list(size=4,alpha=1),ncol=ncl))
  p <- if (!is.null(pal)) p+scale_colour_manual(values=pal,na.value="grey80") else
       p+scale_colour_manual(values=bigpal(length(unique(D[[col]]))))
  if (lab){ cen <- aggregate(cbind(x,y)~get(col), data=D, FUN=median); names(cen)[1]<-"l"
    p <- p+ggplot2::annotate("text",x=cen$x,y=cen$y,label=cen$l,size=3.4,fontface="bold") }
  p }

sf(um("cl", sprintf("Harmony primary clusters (resolution 1.0, %d clusters)", length(clusters)),
   "Reference figure: the cluster identifiers used throughout M14/M15 and in CCC_ANNOTATION_MAPPING.tsv",
   lab=TRUE)+theme(legend.position="none"),
   "01_harmony_clusters", 10.5, 8.5, "UMAP", "postint_harmony_primary_cluster", "res=1.0;dims=1:30;labelled")
sf(um("det","Detailed literature-supported annotation (Level 2)",
   "The SOURCE OF TRUTH. Retained unchanged in the final object alongside annotation_ccc.", ncl=1),
   "02_detailed_annotation", 12.5, 8.5, "UMAP", "postint_celltype_level2_refined", "res=1.0;detailed level2")
sf(um("ccc","CCC-oriented annotation (annotation_ccc)",
   sprintf("Evidence-supported malignant states collapsed into MPNST-Tumor (%d cells, %.1f%%); immune, stromal and endothelial identities retained. MPNST-Tumor is NOT 'everything non-immune'.",
           sum(ccc=="MPNST-Tumor"), 100*sum(ccc=="MPNST-Tumor")/n_cells), pal=pal_ccc, ncl=1),
   "03_ccc_annotation", 12, 8.5, "UMAP", "annotation_ccc", "res=1.0;annotation_ccc", "Primary CCC figure")
D$hl <- ifelse(D$ccc=="MPNST-Tumor","MPNST-Tumor",
        ifelse(D$ccc=="Candidate-Malignant-Unresolved","Candidate malignant (NOT collapsed)","All other cells"))
sf(ggplot(D[D$hl=="All other cells",],aes(x,y))+geom_point(colour="grey86",size=0.22,stroke=0)+
   geom_point(data=D[D$hl!="All other cells",],aes(colour=hl),size=0.45,alpha=0.9,stroke=0)+
   scale_colour_manual(values=c(`MPNST-Tumor`="#B2182B",`Candidate malignant (NOT collapsed)`="#F4A582"))+
   guides(colour=guide_legend(override.aes=list(size=4.5,alpha=1)))+
   labs(title="MPNST-Tumor compartment highlighted",
        subtitle=sprintf("Red = the %d cells collapsed into MPNST-Tumor. Orange = %d 'candidate malignant' cells deliberately NOT collapsed (Low confidence, patient-private evidence only).",
                         sum(ccc=="MPNST-Tumor"), sum(ccc=="Candidate-Malignant-Unresolved")),
        x="UMAP 1",y="UMAP 2",colour=NULL)+th,
   "04_mpnst_tumor_highlight", 11, 8.5, "UMAP_highlight", "annotation_ccc", "highlight=MPNST-Tumor", "Key figure")
imm <- D$comp=="Immune"
sf(ggplot(D[!imm,],aes(x,y))+geom_point(colour="grey88",size=0.2,stroke=0)+
   geom_point(data=D[imm,],aes(colour=ccc),size=0.42,alpha=0.88,stroke=0)+
   scale_colour_manual(values=pal_ccc)+guides(colour=guide_legend(override.aes=list(size=4.5,alpha=1)))+
   labs(title="Immune compartment cell types",
        subtitle=sprintf("%d immune cells across %d identities; non-immune cells shown in grey. Annotations were not altered for visualisation.",
                         sum(imm), length(unique(D$ccc[imm]))), x="UMAP 1",y="UMAP 2",colour=NULL)+th,
   "05_immune_celltypes", 11.5, 8.5, "UMAP_highlight", "annotation_ccc", "immune compartment emphasised")
sf(um("comp","Tumour-microenvironment compartments",
   "Broad CCC compartments. 'Other' is the technical mitochondrial-high cluster; 'Uncertain' is retained, never absorbed into tumour.",
   pal=pal_comp),
   "06_tme_compartments", 10.5, 8.5, "UMAP", "annotation_ccc_compartment", "broad compartments")

panels <- list(
  `MPNST / Schwann-lineage`=c("SOX10","S100B","PLP1","MPZ","NGFR","PMP22","ERBB3","L1CAM","GFRA3","SHH","APOD","CRYAB"),
  `T cell`=c("CD3D","CD3E","CD3G","CD2","TRAC","CD8A","CD8B","IL7R","CD4"),
  `NK`=c("NKG7","GNLY","KLRD1","KLRF1","NCAM1"),
  `Macrophage / monocyte`=c("LYZ","CD68","AIF1","CSF1R","C1QA","CD163","FOLR2","IL1B","S100A8","S100A9","FCN1"),
  `Dendritic`=c("CD1C","FCER1A","CLEC10A","FLT3","LILRA4","CLEC4C","SPIB"),
  `B / plasma`=c("MS4A1","CD79A","POU2AF1","MZB1","JCHAIN","XBP1","FCRL5"),
  `Fibroblast / stromal`=c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA","PI16","COMP"),
  `Pericyte / VSMC`=c("ACTA2","RGS5","PDGFRB","MYH11","NOTCH3"),
  `Endothelial`=c("PECAM1","VWF","CDH5","CLDN5","ACKR1","DLL4"),
  `Proliferation`=c("MKI67","TOP2A","UBE2C","CDK1"))
panels <- lapply(panels, function(g) intersect(g, rownames(obj)))
panels <- panels[vapply(panels, length, 1L) > 0]
dg <- unique(unlist(panels))
ord_ccc <- c("MPNST-Tumor","Candidate-Malignant-Unresolved","CD8-T","CD4-T","NK","T-cell-other",
  "Macrophage","Monocyte","Dendritic","Plasmacytoid-DC","B-cell","Plasma-cell",
  "Fibroblast","Pericyte-VSMC","Endothelial","Uncertain","Low-quality-excluded")
ord_ccc <- ord_ccc[ord_ccc %in% unique(ccc)]
obj$annotation_ccc <- factor(obj$annotation_ccc, levels=ord_ccc)
Idents(obj) <- obj$annotation_ccc
sf(DotPlot(obj, features=dg, assay="SCT", cluster.idents=FALSE)+
   scale_colour_viridis_c(option="viridis")+
   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7.5), axis.text.y=element_text(size=9))+
   labs(title="Canonical marker evidence for the CCC annotation",
        subtitle="Panels for MPNST/Schwann-lineage, T, NK, macrophage/monocyte, dendritic, B/plasma, fibroblast, pericyte, endothelial and proliferation. Expression from the SCT assay - never from Harmony coordinates.",
        x=NULL,y="annotation_ccc"),
   "07_ccc_annotation_marker_dotplot", max(15,0.24*length(dg)), 8.5, "DotPlot", "annotation_ccc",
   sprintf("assay=SCT;genes=%d", length(dg)), "Evidence figure")
avg <- AverageExpression(obj, features=dg, assays="SCT", layer="data", group.by="ident")$SCT
zz <- t(scale(t(as.matrix(avg))))
hm <- reshape2::melt(zz, varnames=c("gene","ccc"), value.name="z")
hm$gene <- factor(hm$gene, levels=rev(dg)); hm$ccc <- factor(as.character(hm$ccc), levels=ord_ccc)
sf(ggplot(hm,aes(ccc,gene,fill=z))+geom_tile()+
   scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,name="z-score")+
   theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=45,hjust=1,size=9),
     axis.text.y=element_text(size=7.5), panel.grid=element_blank())+
   labs(title="Canonical marker heatmap across annotation_ccc",
        subtitle="z-scored mean SCT expression per CCC identity", x=NULL,y=NULL),
   "08_ccc_annotation_marker_heatmap", 11, 13, "Heatmap", "annotation_ccc", "z-scored AverageExpression")

sdf <- summ; sdf$annotation_ccc <- factor(sdf$annotation_ccc, levels=rev(sdf$annotation_ccc[order(sdf$n_cells)]))
sf(ggplot(sdf,aes(annotation_ccc,n_cells,fill=annotation_ccc))+geom_col()+
   geom_text(aes(label=n_cells),hjust=-0.08,size=3.2)+coord_flip(clip="off")+
   expand_limits(y=max(sdf$n_cells)*1.2)+scale_fill_manual(values=pal_ccc,guide="none")+
   geom_hline(yintercept=100,linetype="dashed",colour="grey30")+
   labs(title="CCC population sizes",
        subtitle="Dashed line = the 100-cell audit threshold below which a population is flagged TOO_RARE for communication inference",
        x=NULL,y="cells")+th,
   "09_ccc_population_sizes", 10, 7, "Barplot", "annotation_ccc", "counts")
sf(ggplot(sdf,aes(annotation_ccc,percent_cells,fill=annotation_ccc))+geom_col()+
   geom_text(aes(label=sprintf("%.2f%%",percent_cells)),hjust=-0.08,size=3.2)+coord_flip(clip="off")+
   expand_limits(y=max(sdf$percent_cells)*1.25)+scale_fill_manual(values=pal_ccc,guide="none")+
   labs(title="CCC population proportions", subtitle=sprintf("%d cells total", n_cells), x=NULL,y="% of all cells")+th,
   "10_ccc_population_proportions", 10, 7, "Barplot", "annotation_ccc", "proportions")
prp$annotation_ccc <- factor(as.character(prp$annotation_ccc), levels=ord_ccc)
sf(ggplot(prp,aes(sample,pct_of_sample,fill=annotation_ccc))+geom_col()+scale_fill_manual(values=pal_ccc)+
   labs(title="CCC composition by sample", subtitle=paste("DESCRIPTIVE ONLY - no differential-abundance test.", "sample_id = dataset = patient."), x=NULL,y="% of sample",fill=NULL)+th,
   "11_ccc_composition_by_sample", 10, 7, "StackedBarplot", "annotation_ccc", "pct_of_sample")
sf(ggplot(prp,aes(sample,pct_of_sample,fill=annotation_ccc))+geom_col()+scale_fill_manual(values=pal_ccc)+
   labs(title="CCC composition by patient", subtitle="One sample = one patient; identical to the by-sample figure by construction. DESCRIPTIVE ONLY.", x=NULL,y="% of patient",fill=NULL)+th,
   "12_ccc_composition_by_patient", 10, 7, "StackedBarplot", "annotation_ccc", "pct_of_patient")
tc <- tumor_comp; tc$detailed_tumor_annotation <- factor(tc$detailed_tumor_annotation,
  levels=tc$detailed_tumor_annotation[order(tc$cells)])
sf((ggplot(tc,aes(detailed_tumor_annotation,cells,fill=detailed_tumor_annotation))+geom_col()+
    geom_text(aes(label=sprintf("%d (%.1f%%)",cells,percent_of_mpnst_tumor)),hjust=-0.06,size=3.2)+
    coord_flip(clip="off")+expand_limits(y=max(tc$cells)*1.3)+
    scale_fill_manual(values=bigpal(nrow(tc)),guide="none")+
    labs(title="MPNST-Tumor internal composition",
         subtitle="The detailed malignant states collapsed into MPNST-Tumor. Retained in the final object so tumour heterogeneity is not permanently hidden.",
         x=NULL,y="cells")+th) /
   (ggplot(reshape2::melt(as.matrix(table(det_l2[tsel], grp[tsel])), varnames=c("state","sample"), value.name="n"),
      aes(sample,n,fill=state))+geom_col(position="fill")+scale_fill_manual(values=bigpal(nrow(tc)))+
      labs(title="Tumour-state contribution per sample", subtitle="Within MPNST-Tumor only", x=NULL,y="fraction",fill=NULL)+th),
   "14_mpnst_tumor_internal_composition", 11, 10, "Barplot", "postint_celltype_level2_refined",
   "internal composition of MPNST-Tumor", "Prevents the collapsed tumour population becoming a black box")
sf(um("conf","Detailed annotation confidence",
   "Confidence attaches to the DETAILED annotation and travels with the CCC label. Low/Uncertain regions must not be read as established cell types.",
   pal=pal_conf),
   "15_annotation_confidence", 10.5, 8.5, "UMAP", "postint_annotation_confidence_refined", "confidence_refined")
writeLines(c("13_ccc_composition_by_condition.pdf was NOT generated.","",
 "No biological-condition field exists anywhere in this dataset (M10 audit of all 93 metadata columns).",
 "Fabricating one to satisfy a filename would be scientifically invalid. Figures 11 and 12 already cover",
 "every grouping the data supports, since sample_id is simultaneously the dataset and the patient.","",
 "See reports/phase2/HARMONY_VARIABLE_DECISION.md section 2."),
 file.path(fig_dir,"13_ccc_composition_by_condition_NOT_APPLICABLE.txt"))
write.table(do.call(rbind, frows), file.path(out_dir,"figure_index_m15a_ccc.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# ---------------------------------------------------------------------------
# 7. Save + round-trip
# ---------------------------------------------------------------------------
obj$annotation_ccc <- as.character(obj$annotation_ccc)
Idents(obj) <- factor(obj$annotation_ccc, levels=ord_ccc)
log_info(sprintf("Saving CCC-annotated object to %s ...", out_rds), stage=STAGE)
saveRDS(obj, out_rds)
o_md5 <- digest(out_rds,file=TRUE,algo="md5"); o_sha <- digest(out_rds,file=TRUE,algo="sha256")
log_info(sprintf("Saved %.2f GB | md5 %s", file.info(out_rds)$size/1024^3, o_md5), stage=STAGE)
ref <- digest(obj$annotation_ccc, algo="md5"); rm(obj); invisible(gc(verbose=FALSE))
o2 <- readRDS(out_rds)
V <- list(cells=ncol(o2)==n_cells,
  ccc_present="annotation_ccc" %in% colnames(o2@meta.data),
  ccc_no_missing=all(!is.na(o2$annotation_ccc)) && all(nzchar(o2$annotation_ccc)),
  ccc_identical=identical(digest(o2$annotation_ccc,algo="md5"), ref),
  ccc_compartment_present="annotation_ccc_compartment" %in% colnames(o2@meta.data),
  ccc_is_tumor_present="annotation_ccc_is_tumor" %in% colnames(o2@meta.data),
  ccc_ready_present="annotation_ccc_ccc_ready" %in% colnames(o2@meta.data),
  detailed_preserved=all(sapply(names(det_digests), function(c)
    is.na(det_digests[[c]]) || identical(digest(o2@meta.data[[c]],algo="md5"), det_digests[[c]]))),
  no_parent_label=!any(o2$annotation_ccc=="T-NK-parent"),
  tumor_only_from_map=all(o2@meta.data[[DET_L1]][o2$annotation_ccc=="MPNST-Tumor"]=="Malignant / tumour"),
  uncertain_not_tumor=!any(o2$annotation_ccc[o2@meta.data[[DET_L1]] %in% c("Uncertain","Other")]=="MPNST-Tumor"),
  fibro_endo_not_tumor=!any(o2$annotation_ccc[o2@meta.data[[DET_L1]] %in% c("Fibroblast/Stromal","Endothelial")]=="MPNST-Tumor"),
  immune_not_tumor=!any(o2$annotation_ccc[o2@meta.data[[DET_L1]] %in% c("T/NK","Myeloid","B/Plasma")]=="MPNST-Tumor"),
  all_resolutions=all(sprintf("postint_harmony_clusters_res_%s", format(seq(0.1,1.0,by=0.1),nsmall=1)) %in% colnames(o2@meta.data)),
  reductions=all(c("pca","umap_preintegration","postint_harmony","postint_umap_harmony") %in% Reductions(o2)))
for (nm in names(V)) log_info(sprintf("  %-26s : %s", nm, V[[nm]]), stage=STAGE)
if (!all(unlist(V))) fail(paste("Round-trip validation failed:", paste(names(V)[!unlist(V)],collapse=", ")))
log_info("PASS - annotation_ccc is present and valid; the detailed annotation is byte-identical; no prohibited population was absorbed into MPNST-Tumor.", stage=STAGE)
rm(o2); invisible(gc(verbose=FALSE))

rec <- list(milestone="M15A", phase="phase2", amendment="CCC-oriented annotation layer",
  timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/R/phase2/build_ccc_annotation.R",
  source_object=input_rds, source_checksum=in_md5,
  source_annotation=c(DET_L1,DET_L2,DET_L3,DET_CONF),
  mapping_table=map_tsv, mapping_table_md5=map_md5,
  output_object=out_rds, output_md5=o_md5, output_sha256=o_sha,
  new_columns=c("annotation_ccc","annotation_ccc_compartment","annotation_ccc_is_tumor","annotation_ccc_ccc_ready"),
  tumor_collapse_criteria="Keyed on the detailed Level 2 label. Collapsed: MPNST-like malignant (SCP-like), MPNST-like malignant (NC-like), Schwann-lineage tumour-like, Cycling tumour-like. NOT collapsed: Candidate malignant (provisional label, Low confidence, patient-private evidence only) and everything outside Level 1 'Malignant / tumour'. The inference 'non-immune => tumour' was never used.",
  tnk_subresolution=list(applied_to="the single T/NK parent cluster only",
    logic="NK = CD3(D/E/G)==0 & (NKG7|GNLY|KLRD1|NCAM1|KLRF1)>0; CD8-T = CD3>0 & (CD8A|CD8B)>0; CD4-T = CD3>0 & CD8==0 & (CD4|IL7R|CCR7)>0; else T-cell-other",
    counts=as.list(setNames(as.integer(tnk_counts), names(tnk_counts)))),
  population_counts=as.list(setNames(as.integer(tb), names(tb))),
  population_percent=as.list(setNames(round(100*as.integer(tb)/n_cells,3), names(tb))),
  mpnst_tumor=list(cells=sum(ccc=="MPNST-Tumor"), percent=round(100*sum(ccc=="MPNST-Tumor")/n_cells,3),
    clusters_collapsed=mapping$original_cluster[mapping$collapsed_to_mpnst_tumor=="TRUE"],
    detailed_states=tumor_comp$detailed_tumor_annotation),
  recomputation_avoided=c("Harmony (M11)","integration assessment (M12)","clustering (M13)","marker discovery (M14)"),
  round_trip_validation=V, warnings=warns,
  slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(),t_start,units="secs")), session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir,"m15a_ccc_annotation_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase2_m15a_ccc_annotation", inputs=list(m16_rds=input_rds, ccc_map=map_tsv),
  outputs=list(ccc_rds=out_rds, mapping=file.path(out_dir,"CCC_ANNOTATION_MAPPING.tsv"),
    summary=file.path(out_dir,"CCC_ANNOTATION_SUMMARY.tsv"),
    audit=file.path(out_dir,"CCC_POPULATION_SIZE_AUDIT.tsv"),
    tumor_composition=file.path(out_dir,"MPNST_TUMOR_COMPOSITION.tsv"),
    marker_summary=file.path(out_dir,"CCC_MARKER_SUMMARY.tsv"),
    record=file.path(out_dir,"m15a_ccc_annotation_record.json")),
  parameters=list(milestone="M15A", seed=seed, map_md5=map_md5, output_md5=o_md5,
    slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, file.path(out_dir,"prov_m15a_ccc.json"))
log_info(sprintf("M15A COMPLETE in %.1f min. Warnings: %d", as.numeric(difftime(Sys.time(),t_start,units="mins")), length(warns)), stage=STAGE)
log_system_usage(stage=STAGE)
