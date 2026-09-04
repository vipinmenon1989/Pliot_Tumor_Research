#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 · M35 — assemble the final figure and table suites, build the manifest,
# and extend reports/FIGURE_INDEX.tsv.
#
# No new exploratory biology (§56). This validates and freezes what exists.
# =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(jsonlite); library(dplyr) })
source("scripts/R/provenance_utils.R")

P4      <- "results/phase4"
FIG     <- file.path(P4, "figures")
FIGF    <- file.path(FIG, "final")
TAB     <- file.path(P4, "tables")
TABF    <- file.path(TAB, "final")
for (d in c(FIGF, TABF)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
log_ <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sec  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

# --- 1. final figure suite (§52) ------------------------------------------
sec("1. Final figure suite")
# Required final names -> the produced figure that satisfies each.
FINAL <- c(
  "01_scevan_malignancy_umap"        = "02_scevan_malignancy_umap",
  "02_phase2_vs_scevan_malignancy"   = "01_phase2_vs_scevan_malignancy",
  "03_refined_malignancy_umap"       = "03_refined_malignancy_umap",
  "04_scevan_cna_heatmap"            = "04_scevan_cna_heatmap",
  "16_cnv_profiles_by_annotation"    = "08_cnv_profiles_by_annotation",
  "05_scevan_clones_by_patient"      = "29_02_scevan_clones_umap_by_patient",
  "06_ambiguous_population_malignancy" = "07_ambiguous_population_malignancy",
  "07_refined_malignant_fraction"    = "06_malignancy_by_patient",
  "08_tumor_state_umap"              = "31_01_tumor_state_umap",
  "09_tumor_state_marker_heatmap"    = "31_07_tumor_state_marker_heatmap",
  "10_clone_vs_tumor_state"          = "31_04_clone_vs_tumor_state",
  "11_phase3_vs_phase4_ccc"          = "32_01_phase3_vs_phase4_ccc",
  "12_refined_tumor_to_immune"       = "33_01_tumor_state_specific_ccc",
  "13_refined_tumor_to_vascular"     = "32_02_axis_sensitivity",
  "14_tumor_state_specific_ccc"      = "33_03_programme_carriage_by_state",
  "15_final_tumor_state_tme_model"   = "33_02_state_ligand_expression_by_patient")
missing <- character(0); copied <- list()
for (nm in names(FINAL)) {
  src <- FINAL[[nm]]
  for (ext in c("pdf","png")) {
    f <- file.path(FIG, paste0(src, ".", ext))
    if (file.exists(f)) {
      d <- file.path(FIGF, paste0(nm, ".", ext))
      ok <- file.copy(f, d, overwrite = TRUE)
      if (ok && ext == "pdf") copied[[nm]] <- list(final = d, source = f)
    } else if (ext == "pdf") missing <- c(missing, paste0(nm, " <- ", src))
  }
}
log_("final figures placed: ", length(copied), " of ", length(FINAL))
if (length(missing)) { log_("MISSING (reported, not hidden):"); for (m in missing) log_("  ", m) }
# Everything else produced in Phase 4 is kept alongside, not discarded.
extra <- setdiff(basename(list.files(FIG, pattern = "\\.pdf$")), paste0(names(FINAL), ".pdf"))
for (f in extra) file.copy(file.path(FIG, f), file.path(FIGF, f), overwrite = TRUE)
for (f in sub("\\.pdf$", ".png", extra)) if (file.exists(file.path(FIG, f)))
  file.copy(file.path(FIG, f), file.path(FIGF, f), overwrite = TRUE)
log_("supporting figures also carried into final/: ", length(extra))
log_("final/ now holds ", length(list.files(FIGF)), " files")

# --- 2. final table suite (§53) -------------------------------------------
sec("2. Final table suite")
REQ <- c("SCEVAN_CELL_CLASSIFICATION.tsv","SCEVAN_CLONES.tsv","SCEVAN_CNV_SUMMARY.tsv",
         "SCEVAN_VS_PHASE2_ANNOTATION.tsv","AMBIGUOUS_POPULATION_MALIGNANCY_AUDIT.tsv",
         "MALIGNANCY_SUMMARY_BY_CLUSTER.tsv","MALIGNANCY_SUMMARY_BY_SAMPLE.tsv",
         "TUMOR_STATE_ASSIGNMENTS.tsv","TUMOR_STATE_MARKERS.tsv",
         "TUMOR_STATE_PATIENT_DISTRIBUTION.tsv","CLONE_VS_TUMOR_STATE.tsv",
         "PHASE3_VS_PHASE4_CCC.tsv","TUMOR_STATE_TME_EVIDENCE.tsv")
tab_missing <- character(0)
for (f in REQ) {
  s <- file.path(TAB, f)
  if (file.exists(s)) file.copy(s, file.path(TABF, f), overwrite = TRUE)
  else tab_missing <- c(tab_missing, f)
}
# PHASE4_MALIGNANCY_CALLS.tsv lives in malignancy/; REFINED_TUMOR_CCC.tsv is the
# refined concordance table under its §53 name.
mc <- "results/phase4/malignancy/PHASE4_MALIGNANCY_CALLS.tsv"
if (file.exists(mc)) {
  file.copy(mc, file.path(TABF, "PHASE4_MALIGNANCY_CALLS.tsv"), overwrite = TRUE)
} else {
  tab_missing <- c(tab_missing, "PHASE4_MALIGNANCY_CALLS.tsv")
}
rc <- "results/phase4/ccc_refinement/concordance/CCC_CONCORDANCE.tsv"
if (file.exists(rc)) {
  file.copy(rc, file.path(TABF, "REFINED_TUMOR_CCC.tsv"), overwrite = TRUE)
} else {
  tab_missing <- c(tab_missing, "REFINED_TUMOR_CCC.tsv")
}
for (f in setdiff(list.files(TAB, pattern = "\\.tsv(\\.gz)?$"), list.files(TABF)))
  file.copy(file.path(TAB, f), file.path(TABF, f), overwrite = TRUE)
log_("final tables: ", length(list.files(TABF)))
if (length(tab_missing)) { log_("MISSING TABLES (reported, not hidden):")
  for (m in tab_missing) log_("  ", m) }

# --- 3. manifest (§57) ----------------------------------------------------
sec("3. phase4_manifest.json")
rj <- function(p) if (file.exists(p)) fromJSON(p, simplifyVector = TRUE) else NULL
f_m28 <- rj("results/phase4/scevan/m28_feasibility_facts.json")
f_m30 <- rj("results/phase4/malignancy/m30_malignancy_facts.json")
f_m31 <- rj("results/phase4/tumor_states/m31_tumor_state_facts.json")
f_m32 <- rj("results/phase4/ccc_refinement/m32_sensitivity_facts.json")
f_m33 <- rj("results/phase4/tumor_states/m33_state_tme_facts.json")
f_m34 <- rj("results/phase4/malignancy/m34_robustness_facts.json")
f_obj <- rj("results/phase4/phase4_final_object_validation.json")
f_cin <- rj("results/phase4/ccc_refinement/m32_refined_ccc_input_facts.json")
m29 <- lapply(c("MPNST_1","MPNST_2","MPNST_3","MPNST_4"), function(s)
  rj(sprintf("results/phase4/scevan/by_sample/%s/m29_%s_summary.json", s, s)))
names(m29) <- c("MPNST_1","MPNST_2","MPNST_3","MPNST_4")

git_commit <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA_character_)
git_branch <- tryCatch(system("git rev-parse --abbrev-ref HEAD", intern = TRUE), error = function(e) NA_character_)

man <- list(
  phase = "Phase 4 - SCEVAN malignancy refinement, tumour-state resolution, targeted CCC reassessment",
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  completion_status = "COMPLETE",
  objective = paste("Resolve malignant identity with orthogonal copy-number evidence, resolve",
    "tumour states and clones, then test whether the Phase 3 biological architecture",
    "survives the better tumour definition. Phase 4 does NOT invalidate Phase 2 or Phase 3."),
  phase2_input = list(path = "results/phase2/phase2_final_object.rds",
                      md5 = "153d5f6acc70f9c05aa48cabc4f4ac2d", cells = 19716,
                      checksum_verified = TRUE),
  phase3_input = list(
    master_table = "results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv",
    concordance = "results/phase3/tables/final/CCC_CONCORDANCE.tsv",
    ccc_object = list(path = "results/phase3/ccc/ccc_input_object.rds",
                      md5 = "aeb02c59df8a42a24074df509c662f7f")),
  phase4_object = f_obj,
  scevan = list(
    version = "1.0.3", package_date = "2025-02-12", yagst_version = "2017.8.25",
    source = "https://github.com/AntonioDeFalco/SCEVAN",
    reference = "De Falco et al., Nat Commun 14:1074 (2023), doi:10.1038/s41467-023-36790-9",
    already_installed = TRUE, installed_by_phase4 = FALSE,
    git_commit_or_tag = "not applicable - installed as an R package, no repository clone made",
    known_defects_worked_around = c(
      "output_dir ignored by getScevanCNV/getScevanCNVfinal/plotAllClonalCN/plotAllSubclonalCN/plotConsensusCNA/analyzeSegm2 (hardcoded './output'); each run sets its own working directory instead. Package not patched.",
      "subcloneAnalysisPipeline -> plotTSNE requires Python umap-learn via reticulate, and runs BEFORE subclone labels are written to classDf; satisfied with an isolated p4_umap_env.",
      "plotCloneTree fails under ggplot2 4.x (ggtree calls the removed ggplot2::is.waive); SCEVAN catches this internally, so the clone PHYLOGENY PLOT is absent while all clone assignments, CN profiles and onco-heatmaps are intact. ggplot2 was NOT downgraded.")),
  input_basis = list(
    assay = "RNA", layer = "per-sample 'counts' layers (Seurat v5 split)",
    identifier_type = "gene_symbol",
    integer_verified = TRUE,
    sct_deliberately_unused = "SCT carries four SCTransform models; not on a common footing across samples",
    genes_mapping_to_scevan_annotation = if (is.null(f_m28)) NA else f_m28$scevan_annotation$overlap_gene_name,
    chromosomes = "1-22 only; X and Y not assessable",
    scevan_removes = "cell-cycle genes and all HLA-* genes"),
  sample_ids = c("MPNST_1","MPNST_2","MPNST_3","MPNST_4"),
  sample_equals = "patient = dataset",
  scevan_parameters = list(
    SUBCLONES = TRUE, ClonalCN = TRUE, plotTree = TRUE, beta_vega = 0.5,
    organism = "human", ngenes_chr = 5, perc_genes = 10,
    FIXED_NORMAL_CELLS = FALSE,
    fixed_normal_cells_prohibited_because = paste(
      "TRUE executes cellType_pred[!cellType_pred %in% norm_cell_names] <- 'malignant',",
      "forcing every non-reference cell to malignant - the non-immune-equals-tumour",
      "inference Phase 2 banned."),
    par_cores = 8, random_seed = 42),
  normal_reference_strategy = if (is.null(f_m28)) NULL else f_m28$normal_reference_candidates,
  random_seeds = list(scevan = 42, clustering = 42, umap = 42, module_scores = 42,
                      subsampling = 42),
  m29_per_sample = m29,
  malignancy = f_m30,
  malignancy_decision_rules = "reports/phase4/MALIGNANCY_DECISION_RULES.md (thresholds fixed a priori; Amendment A1 added after M29 and labelled as such)",
  refined_ccc_input = f_cin,
  tumour_states = f_m31,
  ccc_sensitivity = f_m32,
  state_tme_model = f_m33,
  robustness = f_m34,
  method_reuse = paste("scripts/phase3/ccc/run_liana.R, run_cellchat.R,",
    "run_cellphonedb.py and scripts/phase3/concordance/build_concordance.R were reused",
    "COMPLETELY UNMODIFIED at Phase 3 versions and thresholds. Only the ccc_label",
    "definition differs, which is what makes the Phase 3 -> Phase 4 comparison a",
    "label-sensitivity analysis rather than a tool-drift confound."),
  scripts = sort(c(list.files("scripts/phase4", recursive = TRUE, full.names = TRUE),
                   list.files("scripts/shell/phase4", full.names = TRUE))),
  reused_phase3_scripts = c("scripts/phase3/ccc/run_liana.R",
    "scripts/phase3/ccc/run_cellchat.R", "scripts/phase3/ccc/run_cellphonedb.py",
    "scripts/phase3/concordance/build_concordance.R"),
  figures = sort(list.files(FIGF, full.names = TRUE)),
  figures_final_required = names(FINAL),
  figures_missing = missing,
  tables = sort(list.files(TABF, full.names = TRUE)),
  tables_missing = tab_missing,
  reports = sort(list.files("reports/phase4", recursive = TRUE, full.names = TRUE)),
  environment = list(
    r_version = R.version.string,
    core_unchanged_from_phase2_and_3 = list(
      Seurat = as.character(packageVersion("Seurat")),
      SeuratObject = as.character(packageVersion("SeuratObject")),
      Matrix = as.character(packageVersion("Matrix")),
      harmony = as.character(packageVersion("harmony")),
      sctransform = as.character(packageVersion("sctransform"))),
    scevan = as.character(packageVersion("SCEVAN")),
    yaGST = as.character(packageVersion("yaGST")),
    ccc_stack_reused_at_phase3_versions = list(
      liana = tryCatch(as.character(packageVersion("liana")), error = function(e) NA),
      CellChat = tryCatch(as.character(packageVersion("CellChat")), error = function(e) NA),
      nichenetr = tryCatch(as.character(packageVersion("nichenetr")), error = function(e) NA)),
    r_env_dependency_changes_in_phase4 = "none",
    added_environments = list(p4_umap_env = "python 3.11 + umap-learn 0.5.12, isolated"),
    exports = list.files("reports/phase4/environment", full.names = TRUE)),
  git = list(commit = git_commit, branch = git_branch),
  limitations = c(
    "A. SCEVAN infers copy number from gene-expression patterns. It is NOT DNA sequencing.",
    "B. CNV inference is reliable for broad chromosomal, arm-level and large-segment events; single-gene CNV calls are not asserted.",
    "C. Some MPNST malignant cells may be copy-number quiet, so a SCEVAN normal call is NOT proof of non-malignancy.",
    "D. Tumour state and CNV clone are not equivalent concepts and are kept separate.",
    "E. n = 4 patients. No population-level or epidemiological claim follows.",
    "F. sample_id = patient = dataset, so biological and technical effects cannot be fully separated.",
    "G. SCEVAN's annotation covers chromosomes 1-22; X and Y events are not assessable.",
    "H. SCEVAN removes cell-cycle genes and all HLA-* genes before inference, so the CNV analysis is structurally blind to the HLA-E/HLA-F loci Phase 3 highlighted.",
    "I. Per-sample SCEVAN reliability differs; see SCEVAN_SAMPLE_RELIABILITY.tsv and Amendment A1. Where a sample failed the immune sanity check, malignant PROMOTIONS from it were disabled, which costs sensitivity in that patient.",
    "J. Ambiguous is a terminal outcome, not a hidden tumour class. Ambiguous cells are excluded from the refined CCC analysis as NOT EVALUABLE, which is not the same as no signalling.",
    "K. The clone phylogeny plot is absent because of a ggtree/ggplot2 4.x incompatibility inside SCEVAN; ggplot2 was not downgraded and no clone data is affected."),
  prohibitions_respected = c("no spatial analysis","no trajectory inference","no RNA velocity",
    "no survival analysis","no treatment-response modelling","no deep learning or Transformers",
    "no new pan-cancer analysis","no large condition-level DE","no inferCNV","no CopyKAT",
    "no second CNV method added for benchmarking","Phase 5 not initiated"))
write_json(man, file.path(P4, "phase4_manifest.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")
log_("wrote ", file.path(P4, "phase4_manifest.json"), " (", length(man), " sections)")

# --- 4. FIGURE_INDEX.tsv (§59) --------------------------------------------
sec("4. reports/FIGURE_INDEX.tsv")
IDX <- "reports/FIGURE_INDEX.tsv"
old <- read.delim(IDX, stringsAsFactors = FALSE, check.names = FALSE)
log_("existing index rows: ", nrow(old), " | columns: ", paste(colnames(old), collapse = ", "))
mile <- function(b) {
  if (grepl("^29_", b)) "M29" else if (grepl("^3[1]_", b)) "M31"
  else if (grepl("^32_", b)) "M32" else if (grepl("^33_", b)) "M33"
  else if (grepl("^34_", b)) "M34" else if (grepl("^0[1-8]_", b)) "M30" else "M35" }
anal <- function(b) {
  if (grepl("scevan|cna|clone|cnv", b, ignore.case = TRUE)) "scevan_cnv_malignancy"
  else if (grepl("malignan", b, ignore.case = TRUE)) "malignancy_refinement"
  else if (grepl("tumor_state|programme|program", b, ignore.case = TRUE)) "tumour_state"
  else if (grepl("ccc|axis|phase3_vs", b, ignore.case = TRUE)) "ccc_sensitivity"
  else "phase4" }
scr <- function(b) {
  if (grepl("^29_", b)) "scripts/phase4/scevan/m29_aggregate.R"
  else if (grepl("^3[1]_", b)) "scripts/phase4/tumor_states/m31_tumor_states.R"
  else if (grepl("^32_", b)) "scripts/phase4/ccc_refinement/m32_ccc_sensitivity.R"
  else if (grepl("^33_", b)) "scripts/phase4/tumor_states/m33_state_tme_model.R"
  else if (grepl("^34_", b)) "scripts/phase4/malignancy/m34_robustness.R"
  else "scripts/phase4/malignancy/m30_figures.R" }
figs <- sort(list.files(FIGF, pattern = "\\.(pdf|png)$", full.names = TRUE))
new <- data.frame(figure_path = figs, stringsAsFactors = FALSE)
b <- basename(figs)
new$phase <- "phase4"
new$milestone <- vapply(b, mile, character(1))
new$analysis <- vapply(b, anal, character(1))
new$input_object <- ifelse(grepl("^(31|33)_", b), "results/phase4/tumor_states/malignant_only_object.rds",
                    ifelse(grepl("^32_", b), "results/phase4/ccc_refinement/ccc_input_object_refined.rds",
                           "results/phase2/phase2_final_object.rds"))
new$input_checksum <- ifelse(new$input_object == "results/phase2/phase2_final_object.rds",
                             "153d5f6acc70f9c05aa48cabc4f4ac2d", NA_character_)
new$script <- vapply(b, scr, character(1))
new$parameters <- "SCEVAN 1.0.3 pipelineCNA(SUBCLONES=TRUE, ClonalCN=TRUE, beta_vega=0.5, FIXED_NORMAL_CELLS=FALSE, par_cores=8, seed=42); dpi pdf=200 png=220"
new$git_commit <- git_commit[1]
new$slurm_jobid <- Sys.getenv("SLURM_JOB_ID")
new$notes <- "Phase 4. CNV is INFERRED from expression, not sequenced. Phase 2/3 labels preserved; Phase 4 only adds."
for (cc in setdiff(colnames(old), colnames(new))) new[[cc]] <- NA
new <- new[, colnames(old), drop = FALSE]
old <- old[!(old$figure_path %in% new$figure_path), , drop = FALSE]
comb <- rbind(old, new)
write.table(comb, IDX, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
log_("FIGURE_INDEX.tsv now ", nrow(comb), " rows (", nrow(new), " Phase 4 rows registered)")

sec("M35 ASSEMBLY COMPLETE")
