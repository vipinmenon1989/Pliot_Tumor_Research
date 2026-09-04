#!/usr/bin/env Rscript
# Phase 6 - M50 - tables/final, manifest, figure index.
suppressPackageStartupMessages({ library(jsonlite); library(tools) })
source("scripts/phase6/utils/phase6_common.R")
JOBS <- strsplit(Sys.getenv("P6_SLURM_JOBS", ""), ",")[[1]]
GIT  <- tryCatch(system("git rev-parse HEAD", intern = TRUE)[1],
                 error = function(e) NA_character_)
FIG_INDEX <- "reports/FIGURE_INDEX.tsv"; MANIFEST <- "results/phase6/phase6_manifest.json"

p6_sec("1. tables/final")
REQ <- c("CLONE_PROGRAM_ACTIVITY.tsv", "CLONE_PROGRAM_ASSOCIATION.tsv",
         "CLONE_PROGRAM_DIVERSITY.tsv", "PLASTICITY_METRICS.tsv",
         "PROGRAM_TF_ACTIVITY.tsv", "PROGRAM_PATHWAY_ACTIVITY.tsv",
         "PROGRAM_REGULATORY_EVIDENCE.tsv", "BROAD_CNA_EXPRESSION_EFFECTS.tsv",
         "CNA_PROGRAM_ASSOCIATION.tsv", "PATIENT_TUMOR_ARCHITECTURE.tsv",
         "PHASE6_INTEGRATED_EVIDENCE.tsv", "PHASE6_ROBUSTNESS.tsv")
dir.create(P6_TABF, showWarnings = FALSE, recursive = TRUE)
miss_t <- character(0)
for (f in REQ) {
  s <- file.path(P6_TAB, f)
  if (file.exists(s)) file.copy(s, file.path(P6_TABF, f), overwrite = TRUE)
  else miss_t <- c(miss_t, f)
}
if (length(miss_t)) stop("M50: required table(s) missing: ", paste(miss_t, collapse = ", "))
p6_msg("tables/final: %d files", length(list.files(P6_TABF)))

p6_sec("2. Figures")
FIGS <- c("01_clone_program_heatmap", "02_clone_program_alluvial",
          "03_program_diversity_by_clone", "04_within_vs_between_clone_variation",
          "05_tf_activity_by_program", "06_pathway_activity_by_program",
          "07_program_tf_pathway_network", "08_broad_cna_expression_effects",
          "09_cna_program_associations", "10_patient_tumor_architectures",
          "11_malignant_ecm_regulatory_profile", "12_phase6_robustness_summary",
          "13_final_mpnst_malignant_architecture")
have <- character(0); miss <- character(0)
for (f in FIGS) {
  p <- file.path(P6_FIG, paste0(f, ".pdf"))
  if (file.exists(p)) have <- c(have, p) else miss <- c(miss, f)
}
png <- list.files(P6_FIG, pattern = "\\.png$", full.names = TRUE)
allfig <- sort(c(have, png))
p6_msg("required figures: %d/%d | PNG copies: %d", length(have), length(FIGS), length(png))
if (length(miss)) p6_msg("MISSING: %s", paste(miss, collapse = ", "))

p6_sec("3. Figure index")
OBJ <- "results/phase6/phase6_final_object.rds"
IN  <- file.path(P6_CLONE, "phase6_cell_metadata.rds")
INMD5 <- unname(md5sum(IN))
PARAMS <- paste("decoupleR 2.12.0 run_ulm over CollecTRI (TF) and run_mlm over PROGENy top-500 (pathway)",
  "on the frozen RNA log-normalised layer; Hallmark = mean z of MSigDB v2024.1.Hs set genes;",
  "clone-program eta^2 with 1,000 within-patient label permutations; clone-size minimum 20 cells;",
  "dominant-program margin 0.10. MPNST_3 excluded from every clone-based analysis. dpi pdf=200 png=220.")
idx <- read.delim(FIG_INDEX, check.names = FALSE, colClasses = "character", quote = "")
rows <- do.call(rbind, lapply(allfig, function(p) {
  stem <- sub("\\.(pdf|png)$", "", basename(p))
  scr <- if (grepl("^1[0-3]_", stem)) "scripts/phase6/utils/m50_figures2.R"
         else "scripts/phase6/utils/m50_figures.R"
  ms <- if (grepl("^0[12]_", stem)) "M43/M50" else if (grepl("^0[34]_", stem)) "M44/M50"
        else if (grepl("^05_", stem)) "M45/M50" else if (grepl("^0[67]_", stem)) "M46/M50"
        else if (grepl("^0[89]_", stem)) "M47/M50" else if (grepl("^1[01]_", stem)) "M48/M50"
        else if (grepl("^12_", stem)) "M49/M50" else "M48/M50"
  data.frame(figure_path = p, phase = "phase6", milestone = ms, figure_type = "",
    annotation_field = "tumor_clone_phase4 + Phase 5 program scores",
    analysis = "clone_program_regulatory_architecture", sender = "", receiver = "",
    method = "decoupleR/CollecTRI + PROGENy + Hallmark; cNMF program usages from Phase 5",
    dataset = "", processing_stage = "phase6_final", analysis_method = "",
    parameters = PARAMS, input = IN, input_object = OBJ, input_checksum = INMD5,
    generating_script = scr, script = scr, snakemake_rule = "", git_commit = GIT,
    slurm_job_id = paste(JOBS, collapse = ";"),
    notes = paste("Phase 6. MPNST_3's SCEVAN clone structure failed the Phase 4 immune sanity gate and is",
      "EXCLUDED from every clone-based inference; clone labels are patient-scoped and never homologous.",
      "SCEVAN infers copy number from expression, so CNA-expression association is an internal consistency",
      "analysis, NOT independent validation. Within-clone program diversity is not observed state switching."),
    stringsAsFactors = FALSE)
}))
rows <- rows[, names(idx), drop = FALSE]
idx <- idx[!idx$figure_path %in% rows$figure_path, , drop = FALSE]
before <- nrow(idx)
write.table(rbind(idx, rows), FIG_INDEX, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
p6_msg("FIGURE_INDEX.tsv %d -> %d rows", before, before + nrow(rows))

p6_sec("4. Manifest")
readj <- function(p) if (file.exists(p)) fromJSON(p, simplifyVector = FALSE) else NULL
mf5 <- fromJSON("results/phase5/phase5_manifest.json", simplifyVector = FALSE)
mf <- list(
  phase = "Phase 6 - clonal, regulatory and transcriptional-plasticity architecture",
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  completion_status = "COMPLETE",
  objective = paste("Determine whether malignant transcriptional heterogeneity is associated with",
    "distinct CNA-defined clones, with substantial within-clone program diversity, or a mixture,",
    "and which TF and pathway activities distinguish the Phase 5 malignant programs."),
  phase4_input = list(path = P4_OBJECT, md5 = P4_MD5, modified = FALSE,
                      md5_after_phase6 = p6_md5(P4_OBJECT)),
  phase5_input = list(path = P5_OBJECT, md5 = mf5$phase5_object$md5,
                      sha256 = mf5$phase5_object$sha256, modified = FALSE,
                      md5_after_phase6 = p6_md5(P5_OBJECT)),
  mpnst3_exclusion = list(patient = CLONE_EXCLUDED, reason = CLONE_EXCLUSION_REASON,
    scope = "every clone-based inference: clone-program coupling, clone diversity, plasticity, CNA-program association and clone-based conclusions",
    still_described = "MPNST_3's Phase 5 transcriptional programs are described; only its clone architecture is excluded"),
  m42_feasibility = readj(file.path(P6_VAL, "m42_feasibility_facts.json")),
  m43_clone_program = readj(file.path(P6_VAL, "m43_clone_program_facts.json")),
  m44_diversity = readj(file.path(P6_VAL, "m44_diversity_facts.json")),
  m45_tf = readj(file.path(P6_VAL, "m45_tf_activity_facts.json")),
  m46_pathway = readj(file.path(P6_VAL, "m46_pathway_facts.json")),
  m47_cna_expression = readj(file.path(P6_VAL, "m47_cna_expression_facts.json")),
  m48_integrated = readj(file.path(P6_VAL, "m48_integrated_facts.json")),
  m49_robustness = readj(file.path(P6_VAL, "m49_robustness_facts.json")),
  phase6_object = readj(file.path(P6_VAL, "m50_final_object_validation.json"))$phase6_object,
  object_validation = readj(file.path(P6_VAL, "m50_final_object_validation.json"))$validation,
  preservation_guards = readj(file.path(P6_VAL, "m50_final_object_validation.json"))$preservation,
  figures = as.list(allfig), figures_required = as.list(FIGS),
  figures_missing = as.list(miss),
  figure_checksums = lapply(allfig, function(p) list(path = p,
    md5 = unname(md5sum(p)), bytes = as.numeric(file.info(p)$size))),
  tables = as.list(list.files(P6_TABF, full.names = TRUE)),
  table_checksums = lapply(list.files(P6_TABF, full.names = TRUE), function(p)
    list(path = p, md5 = unname(md5sum(p)), bytes = as.numeric(file.info(p)$size))),
  scripts = as.list(c(list.files("scripts/phase6", recursive = TRUE, full.names = TRUE),
                      list.files("scripts/shell/phase6", full.names = TRUE))),
  environment = list(R_env = "UNCHANGED - decoupleR 2.12.0 and OmnipathR 3.14.0 were already installed, so nothing was added for Phase 6",
    networks_cached = list.files("external/networks", full.names = TRUE),
    snapshots = list.files("reports/phase6/environment", full.names = TRUE)),
  slurm_job_ids = as.list(JOBS), git_commit = GIT,
  reports = as.list(list.files("reports/phase6", pattern = "\\.md$",
                               full.names = TRUE, recursive = TRUE)),
  limitations = list(
    "n = 4 patients; sample_id = patient = dataset, so biological and technical effects cannot be fully separated",
    "MPNST_3's SCEVAN clone structure is unreliable and is excluded from every clone-based inference, so clone conclusions rest on 3 patients",
    "SCEVAN infers broad copy number from expression - it is not DNA sequencing",
    "broad CNA-expression association is an INTERNAL CONSISTENCY analysis, not orthogonal validation, because the CNA call was itself derived from expression",
    "cNMF programs are continuous transcriptional patterns, not discrete cell states, cell types, lineages or clones",
    "within-clone program diversity is consistent with plasticity but does NOT demonstrate an actual state transition; no transition rate, direction or trajectory is claimed",
    "cell-level p-values do not represent biological replication; patient-level direction and concordance carry the evidence",
    "clone labels are patient-scoped and were never compared across patients as homologous",
    "no single-gene CNV claim is made anywhere; loci are described only as contained within a broad inferred segment",
    "true generalization requires independent MPNST patients"),
  prohibitions_respected = list(
    "Phase 1-5 objects, annotations, malignancy calls, SCEVAN calls and clones, tumour states, CCC results and Phase 5 programs all unmodified",
    "no new CCC method; LIANA, CellChat, CellPhoneDB, NicheNet and LochNESS were not rerun",
    "no RNA velocity, no pseudotime, no survival analysis, no condition DE",
    "no patient-level predictive machine learning, no deep learning, no Transformers",
    "no spatial, pan-cancer or other-sarcoma comparison",
    "no inferCNV, CopyKAT or any second CNV method",
    "pySCENIC was not introduced; decoupleR was already present in the frozen stack",
    "program, TF activity and pathway activity were kept as independent layers - no weighted composite score and no single numeric rank",
    "R_env not modified",
    "Phase 7 NOT initiated - it requires separate authorization and new biological evidence"))
write(toJSON(mf, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA), MANIFEST)
p6_msg("manifest written: %s (%d sections)", MANIFEST, length(mf))

p6_sec("5. Upstream artefacts untouched")
stopifnot(p6_md5(P4_OBJECT) == P4_MD5, p6_md5(P5_OBJECT) == mf5$phase5_object$md5)
p6_msg("Phase 4 and Phase 5 objects both unchanged")
p6_sec("M50 finalize complete")
