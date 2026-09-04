#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 · M29 — per-patient SCEVAN execution
#
# Runs SCEVAN::pipelineCNA independently for ONE sample (Phase 4 §14:
# sample_id = patient = dataset, so pooling patients first would confuse
# patient-specific CNV architecture with intra-tumour subclones).
#
# Two runs per sample:
#
#   PRIMARY      norm_cell = NULL  -> SCEVAN's own confident-normal detection,
#                with NO input from Phase 2 labels. This is the run that carries
#                inferential weight (§16 anti-circularity).
#
#   SENSITIVITY  norm_cell = high-confidence IMMUNE cells only, and
#                FIXED_NORMAL_CELLS = FALSE. Never the disputed stromal
#                populations, which are the question being asked.
#
# FIXED_NORMAL_CELLS = TRUE is prohibited in this project: SCEVAN's source does
#   cellType_pred[!cellType_pred %in% norm_cell_names] <- "malignant"
# which forces every non-reference cell to malignant - the
# non-immune-equals-tumour inference Phase 2 explicitly banned.
#
# Usage: Rscript m29_run_scevan.R --sample MPNST_1 [--cores 8]
# =============================================================================
suppressPackageStartupMessages({ library(Matrix); library(jsonlite) })
source("scripts/R/provenance_utils.R")

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(f, default = NULL) {
  i <- which(args == paste0("--", f)); if (!length(i)) return(default); args[i + 1]
}
SAMPLE <- getopt("sample")
CORES  <- as.integer(getopt("cores", Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
if (is.null(SAMPLE)) stop("--sample is required")

ROOT     <- normalizePath(".")
SDIR     <- file.path(ROOT, "results/phase4/scevan/by_sample", SAMPLE)
COUNTS   <- file.path(SDIR, "counts_raw.rds")
META     <- file.path(ROOT, "results/phase4/malignancy/phase4_cell_metadata.tsv.gz")

log_ <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sec  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

set.seed(42)
sec(paste0("M29 · SCEVAN · sample ", SAMPLE, " · ", CORES, " cores"))
stopifnot(file.exists(COUNTS), file.exists(META))

cm <- readRDS(COUNTS)
log_("counts: ", nrow(cm), " genes x ", ncol(cm), " cells | class ", class(cm)[1])
xv <- cm@x
stopifnot(all(xv == floor(xv)))           # raw integer counts, re-asserted here
log_("integer-count assertion passed (min ", min(xv), ", max ", max(xv), ")")

md <- read.delim(META, stringsAsFactors = FALSE)
rownames(md) <- md$cell_id
md <- md[colnames(cm), , drop = FALSE]
stopifnot(identical(rownames(md), colnames(cm)),
          all(md$sample_id == SAMPLE))
log_("metadata aligned to count columns; all cells belong to ", SAMPLE)

# --- normal-reference sets ---------------------------------------------------
IMMUNE_HIGH_CONF <- c("CD4-T","CD8-T","NK","T-cell-other","B-cell","Plasma-cell",
                      "Macrophage","Monocyte","Dendritic","Plasmacytoid-DC")
DISPUTED <- c("MPNST-Tumor","Candidate-Malignant-Unresolved","Fibroblast",
              "Pericyte-VSMC","Uncertain","Low-quality-excluded")
immune_ref <- colnames(cm)[md$annotation_ccc %in% IMMUNE_HIGH_CONF]
log_("immune reference candidates in this sample: ", length(immune_ref))
stopifnot(length(intersect(
  colnames(cm)[md$annotation_ccc %in% DISPUTED], immune_ref)) == 0L)
log_("verified: no disputed population contributes to the reference set")

# --- run helper --------------------------------------------------------------
# SCEVAN 1.0.3 ignores output_dir in getScevanCNV/getScevanCNVfinal and the
# plotAll*/plotConsensusCNA/analyzeSegm2 helpers (hardcoded path = "./output"),
# and plotCNclonal() does not forward it. So each run gets its own directory as
# the process working directory and uses SCEVAN's default output_dir, making the
# parameterised and hardcoded paths agree. The package is not patched.
run_scevan <- function(tag, norm_cell) {
  rundir <- file.path(SDIR, tag)
  dir.create(rundir, showWarnings = FALSE, recursive = TRUE)
  sname  <- paste0(SAMPLE, "_", tag)
  sec(paste0("RUN '", tag, "' · sample name '", sname, "' · ",
             if (is.null(norm_cell)) "norm_cell = NULL (automatic confident-normal detection)"
             else paste0("norm_cell = ", length(norm_cell), " immune cells, FIXED_NORMAL_CELLS = FALSE")))
  owd <- setwd(rundir)
  res <- tryCatch({
    t0 <- Sys.time()
    set.seed(42)
    cls <- SCEVAN::pipelineCNA(
      count_mtx  = cm,
      sample     = sname,
      par_cores  = CORES,
      norm_cell  = norm_cell,
      SUBCLONES  = TRUE,
      beta_vega  = 0.5,
      ClonalCN   = TRUE,
      plotTree   = TRUE,
      organism   = "human",
      FIXED_NORMAL_CELLS = FALSE,     # never TRUE - see header
      output_dir = "./output")
    list(ok = TRUE, cls = cls,
         elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
  setwd(owd)

  if (!isTRUE(res$ok)) {
    log_("RUN '", tag, "' FAILED: ", res$message)
    return(list(tag = tag, status = "FAILED", message = res$message))
  }

  cls <- res$cls
  cls$cell_id <- rownames(cls)
  cls$sample_id <- SAMPLE
  cls$scevan_run <- tag
  # Normalise SCEVAN's column names into stable Phase 4 names.
  cls$scevan_call <- cls$class
  cls$scevan_confident_normal <- ifelse(
    !is.null(cls$confidentNormal) & !is.na(cls$confidentNormal), "yes", "no")
  cls$scevan_subclone <- if ("subclone" %in% colnames(cls)) cls$subclone else NA_integer_
  ord <- c("cell_id","sample_id","scevan_run","scevan_call",
           "scevan_confident_normal","scevan_subclone")
  cls <- cls[, c(ord, setdiff(colnames(cls), ord)), drop = FALSE]

  outtsv <- file.path(rundir, paste0("scevan_classification_", SAMPLE, "_", tag, ".tsv"))
  write.table(cls, outtsv, sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(cls, file.path(rundir, paste0("scevan_classification_", SAMPLE, "_", tag, ".rds")))

  log_("elapsed ", round(res$elapsed_min, 2), " min")
  log_("class table:"); print(table(cls$scevan_call))
  log_("confident normals: ", sum(cls$scevan_confident_normal == "yes"))
  if (!all(is.na(cls$scevan_subclone))) {
    log_("subclones:"); print(table(cls$scevan_subclone, useNA = "ifany"))
  } else log_("no subclones reported by SCEVAN for this sample")

  files <- list.files(file.path(rundir, "output"))
  log_("SCEVAN wrote ", length(files), " files into ", tag, "/output/")

  list(tag = tag, status = "SUCCESS", elapsed_min = res$elapsed_min,
       rundir = rundir, sample_name_in_scevan = sname,
       classification_tsv = outtsv,
       cells_in = ncol(cm), cells_returned = nrow(cls),
       class_table = as.list(table(cls$scevan_call)),
       confident_normal = sum(cls$scevan_confident_normal == "yes"),
       subclone_table = if (all(is.na(cls$scevan_subclone))) NULL
                        else as.list(table(cls$scevan_subclone)),
       n_subclones = if (all(is.na(cls$scevan_subclone))) 0L
                     else length(unique(na.omit(cls$scevan_subclone))),
       files_written = files,
       norm_cell_n = if (is.null(norm_cell)) 0L else length(norm_cell),
       parameters = list(par_cores = CORES, SUBCLONES = TRUE, beta_vega = 0.5,
                         ClonalCN = TRUE, plotTree = TRUE, organism = "human",
                         FIXED_NORMAL_CELLS = FALSE, ngenes_chr = 5, perc_genes = 10,
                         seed = 42))
}

results <- list()
results$primary     <- run_scevan("primary", NULL)
results$sensitivity <- run_scevan("sensitivity", immune_ref)

# --- primary vs sensitivity concordance -------------------------------------
sec("PRIMARY vs SENSITIVITY concordance")
if (identical(results$primary$status, "SUCCESS") &&
    identical(results$sensitivity$status, "SUCCESS")) {
  a <- readRDS(file.path(SDIR, "primary",
        paste0("scevan_classification_", SAMPLE, "_primary.rds")))
  b <- readRDS(file.path(SDIR, "sensitivity",
        paste0("scevan_classification_", SAMPLE, "_sensitivity.rds")))
  common <- intersect(a$cell_id, b$cell_id)
  ta <- setNames(a$scevan_call, a$cell_id)[common]
  tb <- setNames(b$scevan_call, b$cell_id)[common]
  ct <- table(primary = ta, sensitivity = tb)
  print(ct)
  agree <- sum(ta == tb) / length(common)
  log_("cells classified by both runs: ", length(common))
  log_(sprintf("call agreement: %.4f", agree))
  results$concordance <- list(cells_common = length(common), agreement = agree,
                              crosstab = as.data.frame.matrix(ct))
} else {
  log_("one or both runs failed; concordance not computed")
}

# --- provenance --------------------------------------------------------------
results$meta <- list(
  sample = SAMPLE, cells = ncol(cm), genes_input = nrow(cm),
  counts_input = COUNTS, counts_md5 = calculate_file_checksum(COUNTS),
  metadata_md5 = calculate_file_checksum(META),
  scevan_version = as.character(packageVersion("SCEVAN")),
  yagst_version = as.character(packageVersion("yaGST")),
  r_version = R.version.string,
  immune_reference_populations = IMMUNE_HIGH_CONF,
  immune_reference_cells = length(immune_ref),
  disputed_excluded_from_reference = DISPUTED,
  fixed_normal_cells = FALSE,
  fixed_normal_cells_rationale = paste(
    "FIXED_NORMAL_CELLS = TRUE forces every cell outside the reference set to",
    "'malignant' (SCEVAN source: cellType_pred[!cellType_pred %in%",
    "norm_cell_names] <- 'malignant'). That is the non-immune-equals-tumour",
    "inference Phase 2 banned, so it is prohibited in this project."),
  scevan_output_dir_defect = paste(
    "SCEVAN 1.0.3 hardcodes path='./output' in getScevanCNV, getScevanCNVfinal,",
    "plotAllClonalCN, plotAllSubclonalCN, plotConsensusCNA and analyzeSegm2, and",
    "plotCNclonal does not forward output_dir. Each run therefore sets its own",
    "working directory and uses SCEVAN's default output_dir. Package not patched."),
  slurm = list(job_id = Sys.getenv("SLURM_JOB_ID"),
               array_task = Sys.getenv("SLURM_ARRAY_TASK_ID"),
               cpus = CORES, node = Sys.info()[["nodename"]]),
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

write_json(results, file.path(SDIR, paste0("m29_", SAMPLE, "_summary.json")),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")
log_("wrote ", file.path(SDIR, paste0("m29_", SAMPLE, "_summary.json")))

failed <- vapply(results[c("primary","sensitivity")],
                 function(x) identical(x$status, "FAILED"), logical(1))
if (any(failed)) stop("SCEVAN run(s) failed for ", SAMPLE, ": ",
                      paste(names(failed)[failed], collapse = ", "))
sec(paste0("M29 COMPLETE for ", SAMPLE))
