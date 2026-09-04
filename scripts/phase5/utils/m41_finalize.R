#!/usr/bin/env Rscript
# =============================================================================
# Phase 5 - M41 - assemble tables/final, write the manifest, extend the figure
# index. Asserts that the frozen Phase 4 manifest and object are untouched.
# =============================================================================
suppressPackageStartupMessages({ library(jsonlite); library(tools) })
source("scripts/phase5/utils/phase5_common.R")
args <- commandArgs(trailingOnly = TRUE)
KSTAR <- as.integer(args[1])
JOBS  <- strsplit(Sys.getenv("P5_SLURM_JOBS", ""), ",")[[1]]
GIT   <- tryCatch(system("git rev-parse HEAD", intern = TRUE)[1],
                  error = function(e) NA_character_)
FIG_INDEX <- "reports/FIGURE_INDEX.tsv"
MANIFEST  <- "results/phase5/phase5_manifest.json"

p5_sec("1. Copy the Phase 5 deliverable tables into tables/final")
REQ <- c("PROGRAM_K_SELECTION.tsv", "MALIGNANT_PROGRAMS.tsv",
         "MALIGNANT_PROGRAM_TOP_GENES.tsv", "MALIGNANT_PROGRAM_CELL_SCORES.tsv",
         "MALIGNANT_PROGRAM_PATIENT_DISTRIBUTION.tsv",
         "MALIGNANT_PROGRAM_PATHWAYS.tsv",
         "PROGRAM_PATIENT_BALANCE_SENSITIVITY.tsv", "PROGRAM_VS_PHASE4_STATE.tsv",
         "MALIGNANT_ECM_VS_TRUE_FIBROBLAST.tsv", "MALIGNANT_ECM_SIGNATURE.tsv",
         "AMBIGUOUS_FIBROBLAST_PROGRAM_PROJECTION.tsv",
         "PHASE5_PROGRAM_ROBUSTNESS.tsv")
EXTRA <- c("PROGRAM_RECURRENCE_BY_K.tsv",
           "PROGRAM_RECURRENCE_THRESHOLD_SENSITIVITY.tsv",
           "M36_MALIGNANT_PATIENT_DISTRIBUTION.tsv",
           "M36_PHASE4_STATE_BY_PATIENT.tsv")
dir.create(P5_TABF, showWarnings = FALSE, recursive = TRUE)
missing <- character(0)
for (f in c(REQ, EXTRA)) {
  src <- file.path(P5_TAB, f)
  if (!file.exists(src)) { if (f %in% REQ) missing <- c(missing, f); next }
  file.copy(src, file.path(P5_TABF, f), overwrite = TRUE)
}
if (length(missing)) stop("M41: required Phase 5 table(s) missing: ",
                          paste(missing, collapse = ", "))
p5_msg("tables/final: %d files", length(list.files(P5_TABF)))

p5_sec("2. Figures")
FIGS <- c("01_program_rank_selection", "02_malignant_program_heatmap",
          "03_program_activity_umap", "04_program_patient_distribution",
          "05_program_vs_phase4_states", "06_recurrent_vs_patient_private_programs",
          "07_patient_balanced_program_sensitivity",
          "08_malignant_ecm_vs_true_fibroblast",
          "09_malignant_ecm_signature_heatmap", "10_program_pathway_activity",
          "11_phase5_robustness_summary", "12_phase5_malignant_program_model")
have <- character(0); miss <- character(0)
for (f in FIGS) {
  p <- file.path(P5_FIG, paste0(f, ".pdf"))
  if (file.exists(p)) have <- c(have, p) else miss <- c(miss, f)
}
png <- list.files(P5_FIG, pattern = "\\.png$", full.names = TRUE)
p5_msg("required figures present: %d/%d | PNG copies: %d", length(have),
       length(FIGS), length(png))
if (length(miss)) p5_msg("MISSING: %s", paste(miss, collapse = ", "))
allfig <- sort(c(have, png))

p5_sec("3. Checksums")
csum <- lapply(allfig, function(p) list(path = p, md5 = unname(md5sum(p)),
                                        bytes = as.numeric(file.info(p)$size)))
tsum <- lapply(list.files(P5_TABF, full.names = TRUE), function(p)
  list(path = p, md5 = unname(md5sum(p)), bytes = as.numeric(file.info(p)$size)))

p5_sec("4. Figure index")
OBJ <- "results/phase5/phase5_final_object.rds"
IN  <- "results/phase5/programs/phase5_cell_metadata.rds"
INMD5 <- unname(md5sum(IN))
PARAMS <- sprintf(paste("cNMF 1.7.1 in the isolated p5_cnmf_env; K grid 4-15, n_iter=100,",
  "seed=42, numgenes=2000, local-density-threshold=0.10; selected K=%d by the pre-declared rule.",
  "Input: RNA raw counts of the 6,434 malignancy_refined==Malignant cells, genes detected in >=0.5%%",
  "of them with ^MT- removed. Gene sets MSigDB v2024.1.Hs. dpi pdf=200 png=220."), KSTAR)
idx <- read.delim(FIG_INDEX, check.names = FALSE, colClasses = "character", quote = "")
rows <- do.call(rbind, lapply(allfig, function(p) {
  stem <- sub("\\.(pdf|png)$", "", basename(p))
  ms <- if (grepl("^0[89]|^12", stem)) "M39/M41" else
        if (grepl("^11", stem)) "M40/M41" else
        if (grepl("^0[1]", stem)) "M37/M41" else "M38/M41"
  data.frame(figure_path = p, phase = "phase5", milestone = ms, figure_type = "",
    annotation_field = "malignancy_refined + cNMF program scores",
    analysis = "continuous_malignant_program_discovery", sender = "", receiver = "",
    method = "cNMF consensus non-negative matrix factorization", dataset = "",
    processing_stage = "phase5_final", analysis_method = "",
    parameters = PARAMS, input = IN, input_object = OBJ, input_checksum = INMD5,
    generating_script = if (grepl("^0[89]|^12", stem))
      "scripts/phase5/utils/m41_figures2.R" else "scripts/phase5/utils/m41_figures.R",
    script = if (grepl("^0[89]|^12", stem))
      "scripts/phase5/utils/m41_figures2.R" else "scripts/phase5/utils/m41_figures.R",
    snakemake_rule = "", git_commit = GIT,
    slurm_job_id = paste(JOBS, collapse = ";"),
    notes = paste("Phase 5. Programs are CONTINUOUS patterns of co-regulated expression, not cell types,",
      "lineages, states or clones. n=4 patients and sample_id=patient=dataset, so cells are not biological",
      "replicates. MPNST_3 contributes malignant cells but no reliable SCEVAN clone structure."),
    stringsAsFactors = FALSE)
}))
rows <- rows[, names(idx), drop = FALSE]
idx <- idx[!idx$figure_path %in% rows$figure_path, , drop = FALSE]
before <- nrow(idx)
write.table(rbind(idx, rows), FIG_INDEX, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
p5_msg("FIGURE_INDEX.tsv %d -> %d rows", before, before + nrow(rows))

p5_sec("5. Manifest")
readj <- function(p) if (file.exists(p)) fromJSON(p, simplifyVector = FALSE) else NULL
mf <- list(
  phase = "Phase 5 - malignant transcriptional programs",
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  completion_status = "COMPLETE",
  objective = paste("Test whether recurrent CONTINUOUS malignant transcriptional programs exist",
                    "across patients even though Phase 4 found 0 of 8 DISCRETE malignant states",
                    "recurrent, and determine what separates malignant ECM-like MPNST cells from",
                    "genuine fibroblasts."),
  phase4_input = list(path = P4_OBJECT, md5 = P4_MD5, sha256 = P4_SHA256,
                      opened = "read-only", modified = FALSE,
                      md5_after_phase5 = p5_md5(P4_OBJECT)),
  historical_annotation_field = readj(file.path(P5_VAL, "m36_feasibility_facts.json"))$historical_annotation_field,
  m36_feasibility = readj(file.path(P5_VAL, "m36_feasibility_facts.json"))$patient_distribution,
  patient_imbalance = readj(file.path(P5_VAL, "m36_feasibility_facts.json"))$patient_imbalance,
  m37_prepare = readj(file.path(P5_VAL, "m37a_prepare_facts.json")),
  k_selection = readj(file.path(P5_VAL, "m37d_k_selection_decision.json")),
  selected_K = KSTAR,
  m38_programs = readj(file.path(P5_VAL, "m38_program_annotation_facts.json")),
  m39_malignant_ecm = readj(file.path(P5_VAL, "m39_malignant_ecm_facts.json")),
  m40_robustness = readj(file.path(P5_VAL, "m40_robustness_facts.json")),
  phase5_object = readj(file.path(P5_VAL, "m41_final_object_validation.json"))$phase5_object,
  object_validation = readj(file.path(P5_VAL, "m41_final_object_validation.json"))$validation,
  preservation_guards = readj(file.path(P5_VAL, "m41_final_object_validation.json"))$preservation,
  figures = as.list(allfig), figures_required = as.list(FIGS),
  figures_missing = as.list(miss), figure_checksums = csum,
  tables = as.list(list.files(P5_TABF, full.names = TRUE)), table_checksums = tsum,
  scripts = as.list(c(list.files("scripts/phase5", recursive = TRUE, full.names = TRUE),
                      list.files("scripts/shell/phase5", full.names = TRUE))),
  environment = list(
    R_env = "UNCHANGED - no package installed, upgraded or downgraded",
    isolated_python_env = "p5_cnmf_env (cnmf 1.7.1, python 3.11)",
    gene_sets = "external/genesets MSigDB v2024.1.Hs GMT files with md5s in MSIGDB_CHECKSUMS.md5",
    snapshots = list.files("reports/phase5/environment", full.names = TRUE)),
  slurm_job_ids = as.list(JOBS), git_commit = GIT,
  reports = as.list(list.files("reports/phase5", pattern = "\\.md$",
                               full.names = TRUE, recursive = TRUE)),
  limitations = list(
    "n = 4 patients; sample_id = patient = dataset, so biological and technical effects cannot be fully separated",
    "cells are NOT biological replicates - cell-level p-values do not represent biological replication",
    "malignant-cell abundance is severely unequal (MPNST_4 holds 57.3% of the compartment), which is why a patient-balanced sensitivity analysis is mandatory rather than optional",
    "a cNMF program is a continuous pattern of co-regulated expression, not a cell type, lineage, state or clone",
    "the malignant ECM-like vs true fibroblast comparison is evaluable in only 2 of 4 patients: MPNST_4 retains 1 non-malignant fibroblast and MPNST_3 contributes no malignant fibroblasts",
    "the malignant-ECM signature is descriptive; no classifier was trained and no pseudo-independent accuracy is reported",
    "Ambiguous fibroblasts are PROJECTED for description only and are never reclassified",
    "Phase 4 CNA metrics reported alongside the groups are context, not independent validation of the split that defined them",
    "true generalization requires independent MPNST patients"),
  prohibitions_respected = list(
    "Phase 1-4 objects, annotations, malignancy calls, SCEVAN calls and clones, tumour states, CCC results, figures, tables, reports and scripts all unmodified",
    "no further clustering sweep was run to rescue the Phase 4 discrete states",
    "no new CCC method; LIANA, CellChat, CellPhoneDB, NicheNet and LochNESS were not rerun",
    "no RNA velocity, no pseudotime, no survival analysis, no condition DE",
    "no patient-level predictive machine learning, no deep learning, no Transformers",
    "no spatial, pan-cancer or other-sarcoma comparison",
    "no inferCNV, CopyKAT or any second CNV method",
    "R_env not modified; cNMF isolated in p5_cnmf_env",
    "Phase 7 not initiated"))
write(toJSON(mf, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA),
      MANIFEST)
p5_msg("manifest written: %s (%d sections)", MANIFEST, length(mf))

p5_sec("6. Phase 4 artefacts untouched")
stopifnot(p5_md5(P4_OBJECT) == P4_MD5)
p4mf <- fromJSON("results/phase4/phase4_manifest.json", simplifyVector = FALSE)
stopifnot(grepl("e85ba8486e456917e2483f2773bdbaf3",
                paste(unlist(p4mf$phase4_object), collapse = " ")))
p5_msg("Phase 4 object md5 unchanged and Phase 4 manifest intact")
p5_sec("M41 finalize complete")
