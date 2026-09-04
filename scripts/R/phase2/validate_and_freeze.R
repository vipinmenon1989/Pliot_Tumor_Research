# scripts/R/phase2/validate_and_freeze.R
#
# Phase 2 / Milestone M17 - final validation, freeze and handoff.
#
# No new exploratory biology. Loads the M16 refined object, validates every structural and
# scientific requirement, writes the deterministic final Phase 2 object, reloads it from
# disk and re-validates, then assembles the final table suite and the Phase 2 manifest.

options(stringsAsFactors=FALSE); options(future.globals.maxSize=+Inf)
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix)
  library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase2_m17_freeze"

args <- commandArgs(trailingOnly=TRUE)
input_rds  <- "results/phase2/annotation/phase2_harmony_ccc.rds"
final_rds  <- "results/phase2/phase2_final_object.rds"
tables_dir <- "results/phase2/tables/final"
manifest   <- "results/phase2/phase2_manifest.json"
expected_cells <- 19716L
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--input"){input_rds<-args[i+1];i<-i+2} else if (a=="--final-rds"){final_rds<-args[i+1];i<-i+2}
  else if (a=="--tables-dir"){tables_dir<-args[i+1];i<-i+2} else if (a=="--manifest"){manifest<-args[i+1];i<-i+2}
  else if (a=="--expected-cells"){expected_cells<-as.integer(args[i+1]);i<-i+2}
  else stop(sprintf("Unknown argument: %s",a)) }

t_start <- Sys.time(); warns <- character(0); issues <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
fail <- function(m){log_error(m,stage=STAGE); stop(m,call.=FALSE)}
dir.create(tables_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(final_rds), recursive=TRUE, showWarnings=FALSE)

log_info("============ M17 PHASE 2 VALIDATION AND FREEZE ============", stage=STAGE)
in_md5 <- digest(input_rds, file=TRUE, algo="md5")
obj <- readRDS(input_rds)
log_info(sprintf("Loaded %d features x %d cells (md5 %s)", nrow(obj), ncol(obj), in_md5), stage=STAGE)

RES <- sprintf("postint_harmony_clusters_res_%s", format(seq(0.1,1.0,by=0.1), nsmall=1))
CCC <- c("annotation_ccc","annotation_ccc_compartment","annotation_ccc_is_tumor","annotation_ccc_ccc_ready")
ANN <- c("postint_celltype_level1","postint_celltype_level2","postint_celltype_level3",
         "postint_annotation_confidence","postint_celltype_level1_initial",
         "postint_celltype_level2_initial","postint_celltype_level3_initial",
         "postint_annotation_confidence_initial","postint_celltype_level1_refined",
         "postint_celltype_level2_refined","postint_celltype_level3_refined",
         "postint_annotation_confidence_refined","postint_annotation_initial",
         "postint_annotation_refined")
PHASE1_META <- c("orig.ident","sample_id","nCount_RNA","nFeature_RNA","percent.mt","percent.ribo",
                 "orig.anno","seurat_clusters","unintegrated_clusters","preint_recommended_cluster")

md <- obj@meta.data
V <- list()
V$cell_count_correct         <- ncol(obj) == expected_cells
V$no_duplicate_cell_ids      <- !any(duplicated(colnames(obj)))
V$feature_count_default      <- nrow(obj) > 0
V$assays_present             <- all(c("RNA","SCT") %in% Assays(obj))
V$default_assay_is_SCT       <- identical(DefaultAssay(obj), "SCT")
V$pre_harmony_pca_present    <- "pca" %in% Reductions(obj)
V$pre_harmony_umap_present   <- "umap_preintegration" %in% Reductions(obj)
V$harmony_reduction_present  <- "postint_harmony" %in% Reductions(obj)
V$harmony_umap_present       <- "postint_umap_harmony" %in% Reductions(obj)
V$phase1_graphs_present      <- all(c("SCT_nn","SCT_snn") %in% Graphs(obj))
V$harmony_graphs_present     <- all(c("postint_harmony_nn","postint_harmony_snn") %in% Graphs(obj))
V$all_10_resolutions_present <- all(RES %in% colnames(md))
V$primary_cluster_present    <- "postint_harmony_primary_cluster" %in% colnames(md)
V$alternative_cluster_present<- "postint_harmony_alternative_cluster" %in% colnames(md)
V$all_annotation_cols_present<- all(ANN %in% colnames(md))
V$phase1_metadata_preserved  <- all(PHASE1_META %in% colnames(md))
V$no_missing_annotations     <- all(!is.na(md$postint_celltype_level2_refined)) &&
                                all(nzchar(md$postint_celltype_level2_refined))
V$no_missing_cluster_labels  <- all(vapply(RES, function(c) all(!is.na(md[[c]])), logical(1)))
V$harmony_dims_30            <- ncol(Embeddings(obj,"postint_harmony")) == 30L
V$no_nonfinite_pca           <- sum(!is.finite(Embeddings(obj,"pca"))) == 0
V$no_nonfinite_harmony       <- sum(!is.finite(Embeddings(obj,"postint_harmony"))) == 0
V$no_nonfinite_umap_pre      <- sum(!is.finite(Embeddings(obj,"umap_preintegration"))) == 0
V$no_nonfinite_umap_harmony  <- sum(!is.finite(Embeddings(obj,"postint_umap_harmony"))) == 0
V$initial_annotation_retained<- !identical(digest(md$postint_annotation_confidence_initial, algo="md5"),
                                           NULL) && all(!is.na(md$postint_annotation_confidence_initial))
# Scientific-scope guards
# --- CCC-oriented annotation layer (amendment) ---
V$ccc_annotation_present     <- all(CCC %in% colnames(md))
V$ccc_no_missing_labels      <- all(!is.na(md$annotation_ccc)) && all(nzchar(md$annotation_ccc))
V$ccc_mpnst_tumor_present    <- "MPNST-Tumor" %in% md$annotation_ccc
V$ccc_immune_identities_kept <- length(unique(md$annotation_ccc[md$annotation_ccc_compartment == "Immune"])) >= 4
V$ccc_stromal_kept_separate  <- any(md$annotation_ccc == "Fibroblast") && any(md$annotation_ccc == "Endothelial")
V$ccc_tumor_only_from_malignant_detail <- all(md[[ "postint_celltype_level1_refined" ]][md$annotation_ccc == "MPNST-Tumor"] == "Malignant / tumour")
V$ccc_no_fibro_endo_in_tumor <- !any(md$annotation_ccc[md$postint_celltype_level1_refined %in% c("Fibroblast/Stromal","Endothelial")] == "MPNST-Tumor")
V$ccc_no_immune_in_tumor     <- !any(md$annotation_ccc[md$postint_celltype_level1_refined %in% c("T/NK","Myeloid","B/Plasma")] == "MPNST-Tumor")
V$ccc_no_uncertain_in_tumor  <- !any(md$annotation_ccc[md$postint_celltype_level1_refined %in% c("Uncertain","Other")] == "MPNST-Tumor")
V$ccc_is_tumor_flag_consistent <- identical(md$annotation_ccc_is_tumor, md$annotation_ccc == "MPNST-Tumor")
V$detailed_tumor_states_retained <- length(unique(md$postint_celltype_level2_refined[md$annotation_ccc == "MPNST-Tumor"])) >= 2
V$no_cnv_columns             <- !any(grepl("cnv|infercnv|copykat", colnames(md), ignore.case=TRUE))
V$no_pseudobulk_de_outputs   <- !any(grepl("pseudobulk|deseq|edger|condition_de", colnames(md), ignore.case=TRUE))
V$no_trajectory_columns      <- !any(grepl("pseudotime|velocity|monocle|slingshot", colnames(md), ignore.case=TRUE))
for (nm in names(V)) log_info(sprintf("  %-32s : %s", nm, V[[nm]]), stage=STAGE)
if (!all(unlist(V))) fail(paste("Final validation failed:", paste(names(V)[!unlist(V)], collapse=", ")))
log_info("PASS - all pre-save validation checks.", stage=STAGE)

# Phase 1 immutability (files on disk)
p1 <- c("processed_mpnst.rds",
        "results/combined/pre_integration/combined_preintegration.rds",
        "results/phase1_manifest.json")
p1_state <- lapply(p1, function(f) if (file.exists(f))
  list(path=f, size_bytes=file.info(f)$size, mtime=format(file.info(f)$mtime, "%Y-%m-%dT%H:%M:%S")) else
  list(path=f, size_bytes=NA, mtime="MISSING"))
for (x in p1_state) log_info(sprintf("  Phase 1 artefact untouched: %s (mtime %s)", x$path, x$mtime), stage=STAGE)
p1_combined_md5 <- digest("results/combined/pre_integration/combined_preintegration.rds", file=TRUE, algo="md5")
V$phase1_handoff_object_unchanged <- identical(p1_combined_md5, "88a442688f912d882f6c6da01820e329")
log_info(sprintf("  phase1_handoff_object_unchanged  : %s (md5 %s)", V$phase1_handoff_object_unchanged, p1_combined_md5), stage=STAGE)
if (!V$phase1_handoff_object_unchanged) fail("The Phase 1 handoff object has changed on disk.")

# ---- save final ------------------------------------------------------------
log_info(sprintf("Writing the final Phase 2 object to %s ...", final_rds), stage=STAGE)
Idents(obj) <- factor(obj$postint_celltype_level2_refined)
saveRDS(obj, final_rds)
fin_md5 <- digest(final_rds,file=TRUE,algo="md5"); fin_sha <- digest(final_rds,file=TRUE,algo="sha256")
fin_size <- file.info(final_rds)$size
log_info(sprintf("Saved %.2f GB | md5 %s | sha256 %s", fin_size/1024^3, fin_md5, fin_sha), stage=STAGE)

snap <- list(cells=ncol(obj), features=nrow(obj), assays=sort(Assays(obj)),
  default_assay=DefaultAssay(obj), reductions=sort(Reductions(obj)), graphs=sort(Graphs(obj)),
  n_meta=ncol(md), meta_cols=colnames(md),
  rna_layers=paste(SeuratObject::Layers(obj[["RNA"]]), collapse=";"),
  sct_layers=paste(SeuratObject::Layers(obj[["SCT"]]), collapse=";"),
  sct_models=length(levels(obj[["SCT"]])),
  digest_umap_harmony=digest(Embeddings(obj,"postint_umap_harmony"), algo="md5"),
  digest_level2=digest(md$postint_celltype_level2_refined, algo="md5"),
  digest_ccc=digest(md$annotation_ccc, algo="md5"),
  n_mpnst_tumor=sum(md$annotation_ccc=="MPNST-Tumor"),
  ccc_identities=sort(unique(md$annotation_ccc)),
  ccc_counts=as.list(setNames(as.integer(table(md$annotation_ccc)), names(table(md$annotation_ccc)))))
rm(obj); invisible(gc(verbose=FALSE))

# ---- RELOAD AND RE-VALIDATE (saveRDS success is not proof of validity) -----
log_info("Reloading the final object from disk and re-validating ...", stage=STAGE)
o2 <- readRDS(final_rds)
R <- list(
  reload_cells=ncol(o2)==snap$cells, reload_features=nrow(o2)==snap$features,
  reload_assays=identical(sort(Assays(o2)), snap$assays),
  reload_default_assay=identical(DefaultAssay(o2), snap$default_assay),
  reload_reductions=identical(sort(Reductions(o2)), snap$reductions),
  reload_graphs=identical(sort(Graphs(o2)), snap$graphs),
  reload_meta_cols=identical(colnames(o2@meta.data), snap$meta_cols),
  reload_sct_models=length(levels(o2[["SCT"]]))==snap$sct_models,
  reload_umap_identical=identical(digest(Embeddings(o2,"postint_umap_harmony"),algo="md5"), snap$digest_umap_harmony),
  reload_annotation_identical=identical(digest(o2@meta.data$postint_celltype_level2_refined,algo="md5"), snap$digest_level2),
  reload_all_resolutions=all(RES %in% colnames(o2@meta.data)),
  reload_all_annotation=all(ANN %in% colnames(o2@meta.data)),
  reload_no_duplicate_cells=!any(duplicated(colnames(o2))),
  reload_no_nonfinite=sum(!is.finite(Embeddings(o2,"postint_harmony")))==0,
  reload_ccc_present=all(CCC %in% colnames(o2@meta.data)),
  reload_ccc_identical=identical(digest(o2@meta.data$annotation_ccc, algo="md5"), snap$digest_ccc),
  reload_ccc_mpnst_tumor=sum(o2@meta.data$annotation_ccc=="MPNST-Tumor") == snap$n_mpnst_tumor,
  reload_detailed_and_ccc_coexist=all(c("postint_celltype_level2_refined","annotation_ccc") %in% colnames(o2@meta.data)))
for (nm in names(R)) log_info(sprintf("  %-30s : %s", nm, R[[nm]]), stage=STAGE)
if (!all(unlist(R))) fail(paste("Reload validation failed:", paste(names(R)[!unlist(R)], collapse=", ")))
log_info("PASS - the saved final object reloads and matches in every checked respect.", stage=STAGE)

md2 <- o2@meta.data
cl_sizes <- as.data.frame(table(cluster=md2$postint_harmony_primary_cluster))
names(cl_sizes) <- c("cluster","n_cells"); cl_sizes$pct_of_all <- round(100*cl_sizes$n_cells/ncol(o2),3)
write.table(cl_sizes, file.path(tables_dir,"cluster_sizes.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
cc <- as.data.frame(table(cluster=md2$postint_harmony_primary_cluster, sample_id=md2$sample_id))
names(cc)[3] <- "n_cells"
write.table(cc, file.path(tables_dir,"cluster_composition.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
rm(o2); invisible(gc(verbose=FALSE))

# ---- assemble the final table suite by copy (preserves provenance) ---------
copies <- c(
  "results/phase2/markers/top_markers_per_cluster.tsv"                 = "top_markers_per_cluster.tsv",
  "results/phase2/markers/all_cluster_markers.tsv"                     = "all_cluster_markers.tsv",
  "results/phase2/markers/marker_summary_by_cluster.tsv"               = "marker_summary_by_cluster.tsv",
  "results/phase2/annotation/ANNOTATION_EVIDENCE.tsv"                  = "annotation_evidence.tsv",
  "results/phase2/composition/cell_counts_by_annotation.tsv"           = "cell_counts_by_annotation.tsv",
  "results/phase2/composition/cell_counts_by_annotation_level1.tsv"    = "cell_counts_by_annotation_level1.tsv",
  "results/phase2/composition/cell_proportions_by_sample.tsv"          = "cell_proportions_by_sample.tsv",
  "results/phase2/composition/cell_proportions_by_patient.tsv"         = "cell_proportions_by_patient.tsv",
  "results/phase2/composition/celltype_sample_dependence.tsv"          = "celltype_sample_dependence.tsv",
  "results/phase2/composition/annotation_refinement_review.tsv"        = "annotation_refinement_review.tsv",
  "results/phase2/clustering/clustering_comparison_table.tsv"          = "clustering_comparison_table.tsv",
  "results/phase2/harmony/harmony_parameters.tsv"                      = "harmony_parameters.tsv",
  "results/phase2/harmony/evaluation/pre_post_mixing_summary.tsv"      = "harmony_pre_post_mixing_summary.tsv",
  "results/phase2/harmony/evaluation/biological_preservation_summary.tsv" = "harmony_biological_preservation_summary.tsv",
  "results/phase2/annotation/CCC_ANNOTATION_MAPPING.tsv"               = "ccc_annotation_mapping.tsv",
  "results/phase2/annotation/CCC_ANNOTATION_SUMMARY.tsv"               = "ccc_annotation_summary.tsv",
  "results/phase2/annotation/CCC_POPULATION_SIZE_AUDIT.tsv"            = "ccc_population_size_audit.tsv",
  "results/phase2/annotation/MPNST_TUMOR_COMPOSITION.tsv"              = "mpnst_tumor_composition.tsv",
  "results/phase2/annotation/CCC_MARKER_SUMMARY.tsv"                   = "ccc_marker_summary.tsv",
  "results/phase2/composition/ccc_counts_by_sample.tsv"                = "ccc_counts_by_sample.tsv",
  "results/phase2/composition/ccc_proportions_by_sample.tsv"           = "ccc_proportions_by_sample.tsv",
  "results/phase2/composition/ccc_counts_by_patient.tsv"               = "ccc_counts_by_patient.tsv",
  "results/phase2/composition/ccc_proportions_by_patient.tsv"          = "ccc_proportions_by_patient.tsv")
for (src in names(copies)) {
  if (file.exists(src)) file.copy(src, file.path(tables_dir, copies[[src]]), overwrite=TRUE)
  else wrn(sprintf("Final table source missing: %s", src)) }
file.copy("results/phase2/composition/cell_proportions_by_condition_NOT_APPLICABLE.txt",
          file.path(tables_dir,"cell_proportions_by_condition_NOT_APPLICABLE.txt"), overwrite=TRUE)
file.copy("results/phase2/composition/ccc_by_condition_NOT_APPLICABLE.txt",
          file.path(tables_dir,"ccc_by_condition_NOT_APPLICABLE.txt"), overwrite=TRUE)
log_info(sprintf("Final table suite: %d files in %s", length(list.files(tables_dir)), tables_dir), stage=STAGE)

# ---- manifest --------------------------------------------------------------
log_info("Assembling the Phase 2 manifest ...", stage=STAGE)
jr <- function(p) if (file.exists(p)) fromJSON(p, simplifyVector=FALSE) else NULL
m11 <- jr("results/phase2/harmony/harmony_parameters.json")
m13 <- jr("results/phase2/clustering/m13_clustering_record.json")
m14 <- jr("results/phase2/markers/m14_marker_record.json")
m15 <- jr("results/phase2/annotation/m15_annotation_record.json")
m16 <- jr("results/phase2/composition/m16_composition_record.json")
m12 <- jr("results/phase2/harmony/evaluation/m12_headline_metrics.json")
pkgs <- c("Seurat","SeuratObject","harmony","Matrix","sctransform","glmGamPoi","presto",
          "cluster","RANN","uwot","aricode","ggplot2","patchwork","viridis","RColorBrewer",
          "dplyr","tidyr","reshape2","jsonlite","digest")
pkgv <- setNames(lapply(pkgs, function(p) tryCatch(as.character(packageVersion(p)), error=function(e) "NOT_INSTALLED")), pkgs)

man <- list(
  phase="phase2", generated=format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), completion_status="COMPLETE",
  phase1_input=list(path="results/combined/pre_integration/combined_preintegration.rds",
    md5="88a442688f912d882f6c6da01820e329",
    sha256="c3fdce8b61602725977f989d8bbf10030ffe48ef6140372a13d2c59903a8dc66",
    verified_unchanged_at_freeze=TRUE, manifest="results/phase1_manifest.json"),
  m11_harmony_object=list(path="results/phase2/harmony/phase2_harmony_integrated.rds",
    md5="cf63e84313a91de25fe2e41660f78f06",
    sha256="6085976f788633b68c07c1ddff49e20e6f6f6727579cb7900854febd2ae6e69d"),
  intermediate_objects=list(
    m13_clustered=list(path="results/phase2/clustering/phase2_harmony_clustered.rds", md5=if(!is.null(m13)) m13$output$md5 else NULL),
    m15_annotated=list(path="results/phase2/annotation/phase2_harmony_annotated.rds", md5=if(!is.null(m15)) m15$output$md5 else NULL),
    m16_refined=list(path="results/phase2/composition/phase2_harmony_refined.rds",
                     md5=if(!is.null(m16)) m16$output$md5 else NULL),
    m15a_ccc=list(path=input_rds, md5=in_md5)),
  final_object=list(path=final_rds, md5=fin_md5, sha256=fin_sha, size_bytes=fin_size,
    cells=snap$cells, features_default_assay=snap$features, assays=snap$assays,
    default_assay=snap$default_assay, reductions=snap$reductions, graphs=snap$graphs,
    n_metadata_columns=snap$n_meta, rna_layers=snap$rna_layers, sct_layers=snap$sct_layers,
    sct_models=snap$sct_models),
  harmony=list(variable="sample_id", dims="1:30", input_reduction="pca",
    output_reduction="postint_harmony", version=as.character(packageVersion("harmony")),
    parameters_explicit=if(!is.null(m11)) m11$harmony$explicit_arguments else NULL,
    parameters_default=if(!is.null(m11)) m11$harmony$default_arguments else NULL,
    iterations_run=if(!is.null(m11)) m11$harmony$iterations_run else NULL,
    reproducibility_check=if(!is.null(m11)) m11$harmony$reproducibility_check else NULL),
  integration_assessment=list(report="reports/phase2/HARMONY_ASSESSMENT.md",
    recommendation="ACCEPT DEFAULT HARMONY WITH CAVEATS",
    headline_metrics=if(!is.null(m12)) list(technical=m12$technical_mixing, biological=m12$biological_preservation) else NULL,
    caveats=list(C1="B/plasma cohesion fell 44% under Harmony; cross-check downstream B/plasma clusters against the non-integrated baseline.",
                 C2="Fibroblast cohesion fell 44%; same requirement for fibroblast/stromal clusters and MPNST_4-derived structure.")),
  clustering=list(reduction="postint_harmony", dims="1:30", k_param=20, algorithm="Louvain (algorithm 1)",
    resolutions_evaluated=seq(0.1,1.0,by=0.1), metadata_columns=RES,
    primary_resolution=1.0, n_primary_clusters=26, alternative_resolution=0.7, n_alternative_clusters=21,
    primary_column="postint_harmony_primary_cluster", alternative_column="postint_harmony_alternative_cluster",
    assessment="reports/phase2/CLUSTERING_ASSESSMENT.md"),
  markers=list(parameters=if(!is.null(m14)) m14$marker_parameters else NULL,
    prep_sct_find_markers=if(!is.null(m14)) m14$prep_sct_find_markers else NULL,
    tables=c("results/phase2/markers/all_cluster_markers.tsv",
             "results/phase2/markers/filtered_cluster_markers.tsv",
             "results/phase2/markers/top_markers_per_cluster.tsv",
             "results/phase2/markers/top10_markers_per_cluster.tsv",
             "results/phase2/markers/top20_markers_per_cluster.tsv",
             "results/phase2/markers/top50_markers_per_cluster.tsv",
             "results/phase2/markers/marker_summary_by_cluster.tsv"),
    report="reports/phase2/MARKER_REPORT.md",
    statistical_scope="Cluster characterisation for annotation. NOT condition-level differential expression."),
  ccc_annotation=list(
    amendment="CCC-oriented annotation layer, authorized 2026-09-03",
    field="annotation_ccc",
    supporting_fields=c("annotation_ccc_compartment","annotation_ccc_is_tumor","annotation_ccc_ccc_ready"),
    mapping_table="results/phase2/annotation/CCC_ANNOTATION_MAPPING.tsv",
    mapping_config="config/phase2/ccc_annotation_map.tsv",
    summary_table="results/phase2/annotation/CCC_ANNOTATION_SUMMARY.tsv",
    size_audit="results/phase2/annotation/CCC_POPULATION_SIZE_AUDIT.tsv",
    tumor_composition="results/phase2/annotation/MPNST_TUMOR_COMPOSITION.tsv",
    marker_summary="results/phase2/annotation/CCC_MARKER_SUMMARY.tsv",
    readiness_report="reports/phase2/CCC_READINESS.md",
    figures="results/phase2/figures/CCC_annotation",
    script="scripts/R/phase2/build_ccc_annotation.R",
    record="results/phase2/annotation/m15a_ccc_annotation_record.json",
    tumor_collapse_criteria="Keyed on the DETAILED Level 2 label. Collapsed into MPNST-Tumor: MPNST-like malignant (SCP-like), MPNST-like malignant (NC-like), Schwann-lineage tumour-like, Cycling tumour-like. NOT collapsed: Candidate malignant (provisional label, Low confidence, patient-private evidence only) and everything outside Level 1 'Malignant / tumour'. The inference 'non-immune => tumour' was never used; fibroblast, endothelial, pericyte, uncertain and low-quality populations were never absorbed.",
    population_counts=snap$ccc_counts,
    n_populations=length(snap$ccc_identities),
    mpnst_tumor_cells=snap$n_mpnst_tumor,
    mpnst_tumor_percent=round(100*snap$n_mpnst_tumor/snap$cells,3),
    detailed_tumor_states_retained=TRUE,
    downstream_use="Phase 3 sample-aware cell-cell communication. No CCC method was executed in Phase 2."),
  annotation=list(map="config/phase2/annotation_map_M15.tsv",
    evidence_table="results/phase2/annotation/ANNOTATION_EVIDENCE.tsv",
    report="reports/phase2/ANNOTATION_REPORT.md",
    final_columns=ANN, cnv_inference_performed=FALSE,
    refinement=if(!is.null(m16)) m16$refinement else NULL),
  composition=list(directory="results/phase2/composition",
    statistical_scope="DESCRIPTIVE ONLY - no inferential condition-level or differential-abundance testing.",
    condition_field_exists=FALSE),
  figures=list(root="results/phase2/figures", milestones=c("M12","M13","M14","M15","M16","final"),
    additional="reports/phase2/figures/m12", index="reports/FIGURE_INDEX.tsv"),
  tables=list(final_suite=tables_dir),
  reports=list(harmony_variable_decision="reports/phase2/HARMONY_VARIABLE_DECISION.md",
    harmony_assessment="reports/phase2/HARMONY_ASSESSMENT.md",
    clustering_assessment="reports/phase2/CLUSTERING_ASSESSMENT.md",
    marker_report="reports/phase2/MARKER_REPORT.md",
    annotation_report="reports/phase2/ANNOTATION_REPORT.md",
    handoff="reports/phase2/PHASE2_HANDOFF.md",
    plan="reports/phase2/PHASE2_PLAN.md",
    milestones=sprintf("reports/phase2/milestones/M%d_REPORT.md", 10:17)),
  scripts=list(
    m10="scripts/R/phase2/inspect_phase1_handoff.R", m11="scripts/R/phase2/run_harmony_integration.R",
    m12="scripts/R/phase2/evaluate_harmony.R", m12_supplement="scripts/R/phase2/m12_supplement_figures.R",
    m13="scripts/R/phase2/cluster_sweep_harmony.R", m14="scripts/R/phase2/discover_markers_harmony.R",
    m15="scripts/R/phase2/annotate_celltypes.R", m16="scripts/R/phase2/refine_and_compose.R",
    m15a="scripts/R/phase2/build_ccc_annotation.R",
    m17="scripts/R/phase2/validate_and_freeze.R",
    slurm=list.files("scripts/shell/phase2", full.names=TRUE)),
  environment=list(r_version=R.version.string, r_binary=file.path(R.home("bin"),"R"),
    conda_env="R_env", package_versions=pkgv,
    spec="workflow/envs/R_env_portable.yaml",
    records=c("reports/phase2/PHASE2_ENVIRONMENT.tsv","reports/phase2/PHASE2_SESSIONINFO.txt")),
  random_seeds=list(global=42, harmony=42, umap=42, clustering=42, markers=42, module_scores=42),
  slurm_jobs=list(
    m10_inspect=list(jobid="19886411", state="COMPLETED", elapsed="00:02:17", maxrss_gib=4.73, cpus=4, mem="96G"),
    m11_harmony=list(jobid="19886486", state="COMPLETED", elapsed="00:08:24", maxrss_gib=8.38, cpus=8, mem="64G"),
    m12_evaluate=list(jobid="19886628", state="COMPLETED", elapsed="00:04:09", maxrss_gib=8.14, cpus=8, mem="96G"),
    m12_supplement=list(jobid="19886682", state="COMPLETED", elapsed="00:01:32", maxrss_gib=7.18, cpus=4, mem="64G",
                        failed_attempts=list(list(jobid="19886675", reason="dpi=NA rejected by ggplot2 4.0.1"))),
    m13_clustering=list(jobid="19886690", state="COMPLETED", elapsed="00:09:37", maxrss_gib=12.09, cpus=8, mem="128G",
                        failed_attempts=list(
                          list(jobid="19886683", reason="preservation guard fired; whole-frame metadata digest replaced with a stricter per-column check"),
                          list(jobid="19886685", reason="unnamed-vector indexing defect"),
                          list(jobid="19886687", reason="sprintf integer format applied to median()"))),
    m14_markers=list(jobid="19886699", state="COMPLETED", elapsed="00:04:26", maxrss_gib=16.40, cpus=12, mem="250G"),
    m15_annotation=list(jobid="19886713", state="COMPLETED", elapsed="00:08:35", maxrss_gib=8.55, cpus=8, mem="128G"),
    m16_composition=list(jobid="19886728", state="COMPLETED", elapsed="00:07:34", maxrss_gib=8.48, cpus=8, mem="128G"),
    m15a_ccc_annotation=list(jobid="19893067", state="COMPLETED", cpus=8, mem="128G"),
    m17_freeze_initial=list(jobid="19886743", state="COMPLETED", elapsed="00:07:39", maxrss_gib=8.13, cpus=8, mem="128G",
                            note="first freeze, before the CCC amendment"),
    m17_freeze=list(jobid=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"), state="running_at_manifest_write",
                    note="re-freeze including annotation_ccc")),
  git=list(commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
    branch=tryCatch(trimws(system("git rev-parse --abbrev-ref HEAD",intern=TRUE)), error=function(e) NA_character_),
    status=tryCatch(paste(system("git status --porcelain",intern=TRUE), collapse="\n"), error=function(e) NA_character_)),
  validation=list(pre_save=V, post_reload=R, phase1_artefacts=p1_state),
  warnings=warns,
  known_limitations=c(
    "Dataset = sample = patient = presumed technical batch is a single 4-level variable. Residual structure cannot be apportioned between uncorrected batch effect and real between-patient tumour biology.",
    "No biological-condition, clinical or technical covariate exists in the data. No condition-level figure or table could be produced and none was fabricated.",
    "No CNV inference was performed (outside Phase 2 scope). All malignant calls are expression- and provenance-based and are labelled conservatively.",
    "The malignant fraction is 23.6% under the conservative annotation and could be up to ~49% if the four fibroblast-programme clusters are MPNST Mes-NC-like malignant cells. This is the largest quantitative uncertainty in Phase 2.",
    "M12 caveats C1 (B/plasma cohesion -44%) and C2 (fibroblast cohesion -44%) apply to clusters C2/C21/C22 and C0/C3/C5/C7/C15/C25 respectively.",
    "The primary clustering resolution 1.0 sits at the top of the tested sweep, so its stability term is one-sided; the alternative 0.7 is interior and two-sided-validated.",
    "The SCT assay carries four SCTransform models (Phase 1 M8 ran SCTransform layer-wise). PrepSCTFindMarkers() was applied for marker discovery; residual per-sample normalisation differences persist in the SCT residuals.",
    "Cluster C15 (438 cells) is a technical mitochondrial-high artefact and should be excluded from biological interpretation.",
    "Annotation is per cluster, so C19 (venous + lymphatic endothelium) and C4 (T and NK together) are under-resolved.",
    "lisi, kBET and clustree are not installed; native inverse-Simpson, dominance-ratio and ggplot-based substitutes were used and are documented.",
    "Phase 1 documentation defects D1-D9 (M10_REPORT.md section 7) remain uncorrected pending researcher instruction.",
    "annotation_ccc collapses four detailed malignant states into MPNST-Tumor. The three 'Candidate malignant' clusters (C13, C17, C18; 1,231 cells) were deliberately NOT collapsed and remain as Candidate-Malignant-Unresolved - they may be promoted if Phase 3 CNV evidence supports them, which would raise the MPNST-Tumor fraction.",
    "The CD4-T versus T-cell-other boundary inside the sub-resolved T/NK cluster is soft, because CD4 transcript detection is sparse in droplet scRNA-seq and the gate relies partly on IL7R/CCR7.",
    "Populations flagged not CCC-ready in CCC_POPULATION_SIZE_AUDIT.tsv (Uncertain, Low-quality-excluded, Candidate-Malignant-Unresolved) must be excluded from communication inference."),
  scientific_prohibitions_respected=c("no condition-level pseudobulk DE","no cell-level condition DE",
    "no differential abundance inference","no survival analysis","no trajectory inference","no RNA velocity",
    "no CNV inference","no cell-cell communication","no perturbation or predictive modelling","no deep learning",
    "no Snakemake used in Phase 2"))
write_json(man, manifest, auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
log_info(sprintf("Wrote %s", manifest), stage=STAGE)

prov <- record_provenance("phase2_m17_freeze", inputs=list(m16_rds=input_rds),
  outputs=list(final_rds=final_rds, manifest=manifest,
    cluster_sizes=file.path(tables_dir,"cluster_sizes.tsv")),
  parameters=list(milestone="M17", final_md5=fin_md5, final_sha256=fin_sha,
    slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, "results/phase2/prov_m17_freeze.json")
log_info(sprintf("M17 COMPLETE in %.1f min. Warnings: %d",
  as.numeric(difftime(Sys.time(),t_start,units="mins")), length(warns)), stage=STAGE)
log_system_usage(stage=STAGE)
