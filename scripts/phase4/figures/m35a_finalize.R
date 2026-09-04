#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 - M35A - finalize
#   * checksums every new figure and plotting table
#   * appends rows to reports/FIGURE_INDEX.tsv
#   * adds an m35a_figure_consolidation block to results/phase4/phase4_manifest.json
#     and appends the new figure paths to its `figures` list
#
# The final Phase 4 object is NOT touched. Its manifest entry and checksum are
# asserted unchanged before and after this script writes anything.
# =============================================================================
suppressPackageStartupMessages({ library(jsonlite); library(tools) })
source("scripts/phase4/figures/m35a_common.R")

FIG_INDEX <- "reports/FIGURE_INDEX.tsv"
MANIFEST  <- "results/phase4/phase4_manifest.json"
JOBIDS    <- c("19899301", "19899307", "19899308", "19899311", "19899312")
JOB_FINAL <- Sys.getenv("M35A_FINAL_JOBID", unset = "19899312")

GIT <- tryCatch(system("git rev-parse HEAD", intern = TRUE)[1],
                error = function(e) NA_character_)
m35a_msg("git commit: %s", GIT)

# ---- the figures M35A produced, with their provenance -----------------------
FIGS <- list(
  list(stem = "17_scevan_native_cna_MPNST1", png = FALSE,
       script = "scripts/phase4/figures/m35a_native_figures.R",
       analysis = "scevan_native_cna_heatmap",
       input = paste("results/phase4/scevan/by_sample/MPNST_1/primary/output/*.png",
                     "results/phase4/malignancy/phase4_malignancy_metadata.rds", sep = "; "),
       notes = "Native SCEVAN CNA heatmaps embedded VERBATIM as rasters; no CNA value redrawn. Page 1 subclone-annotated heatmap + companion Phase 2/clone/malignancy panel; page 2 all-cell, tumour-only and consensus profile. Primary run only."),
  list(stem = "18_scevan_native_cna_MPNST2", png = FALSE,
       script = "scripts/phase4/figures/m35a_native_figures.R",
       analysis = "scevan_native_cna_heatmap",
       input = paste("results/phase4/scevan/by_sample/MPNST_2/primary/output/*.png",
                     "results/phase4/malignancy/phase4_malignancy_metadata.rds", sep = "; "),
       notes = "As figure 17, for MPNST_2. All four MPNST_2 clones are fibroblast-dominated."),
  list(stem = "19_scevan_native_cna_MPNST4", png = FALSE,
       script = "scripts/phase4/figures/m35a_native_figures.R",
       analysis = "scevan_native_cna_heatmap",
       input = paste("results/phase4/scevan/by_sample/MPNST_4/primary/output/*.png",
                     "results/phase4/malignancy/phase4_malignancy_metadata.rds", sep = "; "),
       notes = "As figure 17, for MPNST_4. Seven of eight MPNST_4 clones are fibroblast-dominated; clone 8 is endothelial."),
  list(stem = "20_scevan_cna_reliable_patients", png = FALSE,
       script = "scripts/phase4/figures/m35a_native_figures.R",
       analysis = "scevan_native_cna_heatmap",
       input = "results/phase4/scevan/by_sample/{MPNST_1,MPNST_2,MPNST_4}/primary/output/*heatmap_subclones.png",
       notes = "Three reliable patients on one page. MPNST_3 excluded (failed immune sanity gate)."),
  list(stem = "21_scevan_clone_composition_phase2", png = TRUE,
       script = "scripts/phase4/figures/m35a_clone_composition.R",
       analysis = "scevan_clone_composition_by_phase2_annotation",
       input = "results/phase4/malignancy/phase4_malignancy_metadata.rds",
       notes = "PRIMARY M35A evidence figure. SCEVAN builds clones from CNA alone and never saw a Phase 2 label, so clone composition is non-circular. MPNST_3 omitted (figure 25)."),
  list(stem = "22_fibroblast_cna_burden_by_patient", png = TRUE,
       script = "scripts/phase4/figures/m35a_fibroblast_cna.R",
       analysis = "fibroblast_cna_burden_by_patient",
       input = "results/phase4/malignancy/phase4_malignancy_metadata.rds",
       notes = "Patient-stratified, never pooled. Effect size is Cliff's delta, not a pooled cell-level p-value. MPNST_3 marked unreliable; MPNST_4 contrast NOT EVALUABLE (1 non-malignant fibroblast)."),
  list(stem = "23_fibroblast_malignant_vs_nonmalignant_cna_profile", png = FALSE,
       script = "scripts/phase4/figures/m35a_fibroblast_cna.R",
       analysis = "fibroblast_cna_profile_comparison",
       input = "results/phase4/scevan/by_sample/*/primary/output/*_CNAmtx.RData; *_count_mtx_annot.RData; results/phase4/malignancy/phase4_malignancy_metadata.rds",
       notes = "Mean per-gene relative CNA from the STORED native SCEVAN CNA matrices. SCEVAN was not re-run. Shared hg38 genome axis. Groups drawn only at n >= 20."),
  list(stem = "24_phase2_to_phase4_malignancy_transition", png = TRUE,
       script = "scripts/phase4/figures/m35a_transition.R",
       analysis = "phase2_to_phase4_identity_transition",
       input = "results/phase4/malignancy/phase4_malignancy_metadata.rds",
       notes = "Disputed compartments only; canonical immune populations deliberately excluded so they cannot bury the transition."),
  list(stem = "25_MPNST3_scevan_failure_qc", png = TRUE,
       script = "scripts/phase4/figures/m35a_mpnst3_qc.R",
       analysis = "mpnst3_scevan_failure_qc",
       input = "results/phase4/scevan/by_sample/*/m29_*_summary.json; results/phase4/tables/SCEVAN_POPULATION_MALIGNANT_FRACTION.tsv; results/phase4/malignancy/phase4_malignancy_metadata.rds",
       notes = "MPNST_3 shown failing, not omitted. Depth shown as CONTEXT ONLY; low depth is NOT claimed as the cause."),
  list(stem = "26_fibroblast_malignancy_threshold_robustness", png = TRUE,
       script = "scripts/phase4/figures/m35a_threshold.R",
       analysis = "malignancy_threshold_sensitivity",
       input = "results/phase4/tables/MALIGNANCY_THRESHOLD_SENSITIVITY.tsv",
       notes = "Fibroblast -> Malignant is 4,036 at every tested pop_frac_low (0.15-0.40); retained MPNST-Tumor and the 32.63% cohort fraction are threshold-sensitive."),
  list(stem = "27_broad_cna_recurrence", png = FALSE,
       script = "scripts/phase4/figures/m35a_broad_cna.R",
       analysis = "broad_cna_recurrence",
       input = "results/phase4/tables/SCEVAN_CNV_SUMMARY.tsv; results/phase4/tables/SCEVAN_RECURRENT_BROAD_EVENTS.tsv",
       notes = "SECONDARY evidence. Broad clonal segments >= 10 Mb only. NO gene-level CNV claim; NF1/NF2 marked only as 'broad segment containing the locus'."),
  list(stem = "28_phase4_scevan_evidence_summary", png = TRUE,
       script = "scripts/phase4/figures/m35a_summary.R",
       analysis = "phase4_scevan_evidence_summary",
       input = "results/phase4/malignancy/phase4_malignancy_metadata.rds; results/phase4/tables/MALIGNANCY_THRESHOLD_SENSITIVITY.tsv; results/phase4/scevan/by_sample/MPNST_2/primary/output/MPNST_2_primaryheatmap_subclones.png",
       notes = "Integrated evidence summary, panels A-H. SCEVAN does NOT provide DNA-level proof; the 32.63% cohort fraction is not patient-robust.")
)

PARAMS <- paste(
  "M35A visualization only - SCEVAN not re-run, no malignancy call, rule, amendment,",
  "clone, tumour state or threshold recomputed. Native SCEVAN CNA heatmaps embedded verbatim.",
  "dpi pdf=200 png=220. Upstream: SCEVAN 1.0.3 pipelineCNA(SUBCLONES=TRUE, ClonalCN=TRUE,",
  "beta_vega=0.5, FIXED_NORMAL_CELLS=FALSE, par_cores=8, seed=42), pop_frac_low=0.25 a priori.")

INPUT_OBJ <- "results/phase4/malignancy/phase4_malignancy_metadata.rds"
INPUT_MD5 <- unname(tools::md5sum(INPUT_OBJ))
m35a_msg("plotting input %s md5 %s", INPUT_OBJ, INPUT_MD5)

# ---- assert the frozen final object is untouched ----------------------------
MANIFEST_BEFORE <- file.path(tempdir(), "phase4_manifest_before.json")
file.copy(MANIFEST, MANIFEST_BEFORE, overwrite = TRUE)
mf <- jsonlite::fromJSON(MANIFEST, simplifyVector = FALSE)
FINAL_MD5 <- "e85ba8486e456917e2483f2773bdbaf3"
obj_md5 <- unlist(mf$phase4_object)[grepl("md5", names(unlist(mf$phase4_object)))]
stopifnot(any(obj_md5 == FINAL_MD5))
m35a_msg("phase4_final_object md5 in manifest unchanged: %s", FINAL_MD5)

# ---- collect paths + checksums ---------------------------------------------
rows <- list(); newpaths <- character(0); csum <- list()
for (f in FIGS) {
  ps <- file.path(M35A_FIG, paste0(f$stem, ".pdf"))
  if (isTRUE(f$png)) ps <- c(ps, file.path(M35A_FIG, paste0(f$stem, ".png")))
  for (p in ps) {
    if (!file.exists(p)) stop("M35A: expected figure missing: ", p)
    md5 <- unname(tools::md5sum(p))
    newpaths <- c(newpaths, p)
    csum[[p]] <- list(path = p, md5 = md5,
                      bytes = as.numeric(file.info(p)$size),
                      script = f$script, analysis = f$analysis)
    rows[[length(rows) + 1L]] <- data.frame(
      figure_path = p, phase = "phase4", milestone = "M35A", figure_type = "",
      annotation_field = "", analysis = f$analysis, sender = "", receiver = "",
      method = "SCEVAN 1.0.3 (inferred CNA; re-used, not re-run)", dataset = "",
      processing_stage = "phase4_final", analysis_method = "",
      parameters = PARAMS, input = f$input, input_object = INPUT_OBJ,
      input_checksum = INPUT_MD5, generating_script = f$script, script = f$script,
      snakemake_rule = "", git_commit = GIT, slurm_job_id = JOB_FINAL,
      notes = f$notes, stringsAsFactors = FALSE)
  }
}
rows <- do.call(rbind, rows)
m35a_msg("new figure files: %d", nrow(rows))

# ---- FIGURE_INDEX.tsv -------------------------------------------------------
idx <- read.delim(FIG_INDEX, check.names = FALSE, colClasses = "character",
                  quote = "")
stopifnot(all(names(rows) %in% names(idx)))
rows <- rows[, names(idx), drop = FALSE]
idx <- idx[!idx$figure_path %in% rows$figure_path, , drop = FALSE]   # idempotent
before <- nrow(idx)
idx2 <- rbind(idx, rows)
write.table(idx2, FIG_INDEX, sep = "\t", quote = FALSE, row.names = FALSE,
            na = "")
m35a_msg("FIGURE_INDEX.tsv %d -> %d rows", before, nrow(idx2))

# ---- plotting tables produced by M35A --------------------------------------
TABS <- file.path(M35A_TABF, c("SCEVAN_CLONE_COMPOSITION_VISUALIZATION.tsv",
                               "M35A_FIBROBLAST_CNA_EFFECT_SIZES.tsv",
                               "M35A_FIBROBLAST_CNA_PROFILE_CORRELATION.tsv"))
tab_entries <- lapply(TABS, function(p) {
  if (!file.exists(p)) stop("M35A: expected table missing: ", p)
  list(path = p, md5 = unname(tools::md5sum(p)),
       bytes = as.numeric(file.info(p)$size))
})

# ---- manifest ---------------------------------------------------------------
mf$figures <- as.list(unique(c(unlist(mf$figures), newpaths)))
mf$tables  <- as.list(unique(c(unlist(mf$tables), TABS)))
mf$m35a_figure_consolidation <- list(
  milestone = "M35A",
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  objective = paste("Surface, consolidate and improve the visualization of the EXISTING Phase 4",
                    "SCEVAN evidence so the malignant fibroblast-like / Mesenchymal-ECM conclusion",
                    "is visually defensible."),
  reruns_performed = "none - SCEVAN, malignancy rules, amendments A1/A2, clone assignments, tumour states, CCC and thresholds were all re-used unchanged",
  final_object_untouched = list(path = "results/phase4/phase4_final_object.rds",
                                md5 = FINAL_MD5, loaded = FALSE, modified = FALSE),
  plotting_input = list(path = INPUT_OBJ, md5 = INPUT_MD5, cells = 19716L,
                        note = "the 6 GB final object was never loaded; this 5 MB metadata frame carries every field the figures need"),
  native_scevan_sources = list(
    index = "results/phase4/tables/SCEVAN_NATIVE_FIGURE_INDEX.tsv",
    surfaced = c("heatmap_subclones.png", "heatmap.png", "onlytumorheatmap.png", "consensus.png"),
    not_surfaced = c("OncoHeat.png / OncoHeat2.png - cytoband labels unreadable at page scale and invite gene-level over-reading",
                     "umap_CNA.png / umap_scRNA.png - redundant with existing figures 01/02, 05, 31_03",
                     "CloneTree.png - all 8 files are BLANK (ggtree/ggplot2 4.x defect, handoff limitation K); phylogeny plotting was not forced"),
    embedded_verbatim = TRUE, cna_values_redrawn = FALSE),
  scripts = c("scripts/phase4/figures/m35a_common.R",
              "scripts/phase4/figures/m35a_native_figures.R",
              "scripts/phase4/figures/m35a_clone_composition.R",
              "scripts/phase4/figures/m35a_fibroblast_cna.R",
              "scripts/phase4/figures/m35a_transition.R",
              "scripts/phase4/figures/m35a_mpnst3_qc.R",
              "scripts/phase4/figures/m35a_threshold.R",
              "scripts/phase4/figures/m35a_broad_cna.R",
              "scripts/phase4/figures/m35a_summary.R",
              "scripts/phase4/figures/m35a_finalize.R",
              "scripts/shell/phase4/run_m35a_figures.sh"),
  figures = unname(csum),
  tables = tab_entries,
  slurm_job_ids = JOBIDS,
  slurm_final_job_id = JOB_FINAL,
  slurm_note = "all runs COMPLETED exit 0; the reruns are figure-quality iterations found by inspecting rendered output, not failures. Peak MaxRSS 2.11 GiB against a 24 G request.",
  git_commit = GIT,
  reports = c("reports/phase4/M35A_SCEVAN_FIGURE_AUDIT.md",
              "reports/phase4/PHASE4_HANDOFF.md (section 28)"),
  numerical_validation = list(
    performed_before_every_figure = TRUE,
    checks = 18L,
    total_cells = 19716L,
    refined = list(Malignant = 6434L, `Non-malignant` = 9078L, Ambiguous = 3766L,
                   `Excluded-low-quality` = 438L),
    fibroblast = list(total = 5064L, Malignant = 4036L, `Non-malignant` = 908L,
                      Ambiguous = 120L),
    candidate_malignant_unresolved = list(Malignant = 836L, Ambiguous = 395L,
                                          `Non-malignant` = 0L),
    mpnst_tumor = list(Malignant = 1405L, Ambiguous = 2015L),
    clones = list(MPNST_1 = 7L, MPNST_2 = 4L, MPNST_3 = 3L, MPNST_4 = 8L),
    result = "all checks passed in every figure script"),
  discrepancy_found = list(
    where = "PHASE4_HANDOFF.md section 12",
    text = "'Every MPNST_1 clone mixes Candidate-Malignant-Unresolved, Fibroblast and MPNST-Tumor together'",
    finding = "holds for 6 of 7 clones: MPNST_1_clone6 (176 cells) carries 108 Candidate-Malignant-Unresolved and 68 MPNST-Tumor and NO Fibroblast cells; MPNST_1_clone7 carries 1",
    accurate_statement = "all 7 MPNST_1 clones mix Candidate-Malignant-Unresolved with MPNST-Tumor; 6 of 7 additionally carry Fibroblast cells",
    impact = "none on any call, count, clone assignment or conclusion; the two exceptions are 230 of MPNST_1's 1,952 malignant cells",
    action = "recorded in handoff section 28.7 and in M35A_SCEVAN_FIGURE_AUDIT.md; section 12 left as written; figure 21 plots the composition as it actually is"),
  key_results = list(
    clone_composition = "MPNST_2 4/4 fibroblast-dominated; MPNST_4 7/8 (clone 8 endothelial); MPNST_1 all 7 mix CMU with MPNST-Tumor, 6/7 also carry Fibroblast",
    cliffs_delta_fibroblast_mal_vs_nonmal = list(MPNST_1 = "0.834-0.917 across the four CNA metrics",
                                                 MPNST_2 = "0.781-0.890",
                                                 MPNST_3 = "NOT EVALUABLE - promotions disabled",
                                                 MPNST_4 = "NOT EVALUABLE - 1 non-malignant fibroblast"),
    cna_profile_correlation = list(
      MPNST_1 = "r(Fib->Mal, MPNST-Tumor->Mal) = 0.934 vs r(Fib->Mal, Fib->Non-mal) = 0.184, immune -0.055",
      MPNST_2 = "r(Fib->Mal, Fib->Non-mal) = 0.346, immune -0.387; no evaluable MPNST-Tumor->Malignant group (n=7)",
      MPNST_4 = "r(Fib->Mal, MPNST-Tumor->Mal) = 0.972, immune -0.381"),
    threshold_stability = "Fibroblast -> Malignant = 4036 at pop_frac_low 0.15/0.20/0.25/0.30/0.40; retained MPNST-Tumor 3266 -> 1405; refined fraction 42.07% -> 32.63%",
    mpnst3 = "primary/sensitivity agreement 0.0864; clones are T/NK, plasma/pDC/B and plasma-dominated; 6 of 9 immune populations ~100% malignant; median 1,594 genes/cell shown as context only",
    honest_limit = "the malignant vs non-malignant fibroblast contrast exists in 2 of 4 patients only (MPNST_4 has 1 non-malignant fibroblast, MPNST_3 has no malignant fibroblasts)"),
  prohibitions_respected = c(
    "SCEVAN not re-run", "malignancy_refined / malignancy_confidence unchanged",
    "malignancy decision rules unchanged", "amendments A1/A2 unchanged",
    "SCEVAN calls and clones unchanged", "tumour states unchanged",
    "CCC not re-run", "NicheNet not re-run", "LochNESS not re-run",
    "tumour cells not re-clustered", "cNMF not run", "Phase 5 not started",
    "inferCNV not installed", "CopyKAT not installed", "R environment unchanged",
    "ggplot2 / ggtree not downgraded", "clone phylogeny plotting not forced",
    "phase4_final_object.rds neither loaded nor modified",
    "no single-gene CNV claim; NF1/NF2 described only as broad segments containing the locus",
    "no claim that SCEVAN provides DNA-level proof")
)
# digits = NA keeps FULL precision. jsonlite's default (4) would silently round
# frozen values already in the manifest - e.g. run agreement 0.99514117 -> 0.9951.
write(jsonlite::toJSON(mf, auto_unbox = TRUE, pretty = TRUE, null = "null",
                       digits = NA), MANIFEST)

# ---- re-assert after writing ------------------------------------------------
# Every pre-existing manifest section must survive byte-for-byte in value, not
# just in name. M35A adds; it never edits what M28-M35 froze.
mf2 <- jsonlite::fromJSON(MANIFEST, simplifyVector = FALSE)
o2 <- unlist(mf2$phase4_object)
stopifnot(any(o2[grepl("md5", names(o2))] == FINAL_MD5))
mf0 <- jsonlite::fromJSON(MANIFEST_BEFORE, simplifyVector = FALSE)
untouched <- setdiff(names(mf0), c("figures", "tables"))
bad <- untouched[!vapply(untouched, function(k) identical(mf0[[k]], mf2[[k]]),
                         logical(1))]
if (length(bad)) stop("M35A: pre-existing manifest sections were altered: ",
                      paste(bad, collapse = ", "))
stopifnot(all(unlist(mf0$figures) %in% unlist(mf2$figures)),
          all(unlist(mf0$tables)  %in% unlist(mf2$tables)))
m35a_msg("all %d pre-existing manifest sections verified unchanged", length(untouched))
m35a_msg("manifest written; phase4_object md5 still %s", FINAL_MD5)
m35a_msg("manifest sections: %d | figures listed: %d | tables listed: %d",
         length(mf2), length(mf2$figures), length(mf2$tables))
m35a_msg("== M35A finalize done ==")
