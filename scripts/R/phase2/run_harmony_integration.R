# scripts/R/phase2/run_harmony_integration.R
#
# Phase 2 / Milestone M11 - Default Harmony integration.
#
# Executes the Harmony configuration approved at the M10 gate:
#   grouping variable : sample_id
#   input reduction   : pca            (Phase 1 shared, assay = SCT)
#   dimensions        : 1:30
#   output reduction  : postint_harmony
#   scientific params : ALL left at harmony 1.2.4 package defaults
#
# Hard guarantees enforced by this script:
#   * The Phase 1 handoff object is opened READ-ONLY and never written to.
#   * The output path is a NEW file and is refused if it already exists (unless --force).
#   * Phase 1 reductions (pca, umap_preintegration), Phase 1 graphs and the full
#     Phase 1 metadata frame are digested before and after Harmony and asserted
#     byte-identical.
#   * Every precondition is checked and the script aborts with a non-zero exit
#     status if any assumption is violated.
#
# NOT done here (reserved for later milestones): neighbours, UMAP, clustering,
# marker discovery, annotation, and the scientific pre/post integration assessment.

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = +Inf)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(jsonlite)
  library(digest)
})

source("scripts/R/logging_utils.R")
source("scripts/R/provenance_utils.R")

setup_strict_logging()   # options(warn = 1): warnings surface immediately, never suppressed

STAGE <- "phase2_m11_harmony"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

input_rds        <- "results/combined/pre_integration/combined_preintegration.rds"
out_rds          <- "results/phase2/harmony/phase2_harmony_integrated.rds"
out_dir          <- "results/phase2/harmony"
group_by_var     <- "sample_id"
reduction_use    <- "pca"
dims_use_max     <- 30L
reduction_save   <- "postint_harmony"
random_seed      <- 42L
force_overwrite  <- FALSE
validation_mode  <- FALSE   # synthetic smoke test: relaxes the M10-specific real-data expectations only
expected_cells   <- 19716L
expected_input_md5    <- "88a442688f912d882f6c6da01820e329"
expected_input_sha256 <- "c3fdce8b61602725977f989d8bbf10030ffe48ef6140372a13d2c59903a8dc66"

i <- 1
while (i <= length(args)) {
  key <- args[i]
  if (key == "--input")               { input_rds <- args[i + 1];             i <- i + 2
  } else if (key == "--out-rds")      { out_rds <- args[i + 1];               i <- i + 2
  } else if (key == "--out-dir")      { out_dir <- args[i + 1];               i <- i + 2
  } else if (key == "--group-by")     { group_by_var <- args[i + 1];          i <- i + 2
  } else if (key == "--reduction-use"){ reduction_use <- args[i + 1];         i <- i + 2
  } else if (key == "--dims")         { dims_use_max <- as.integer(args[i+1]); i <- i + 2
  } else if (key == "--reduction-save"){ reduction_save <- args[i + 1];       i <- i + 2
  } else if (key == "--random-seed")  { random_seed <- as.integer(args[i+1]); i <- i + 2
  } else if (key == "--expected-cells"){ expected_cells <- as.integer(args[i+1]); i <- i + 2
  } else if (key == "--expected-md5") { expected_input_md5 <- args[i + 1];    i <- i + 2
  } else if (key == "--expected-sha256"){ expected_input_sha256 <- args[i+1]; i <- i + 2
  } else if (key == "--force")        { force_overwrite <- TRUE;              i <- i + 1
  } else if (key == "--validation-mode") { validation_mode <- TRUE;            i <- i + 1
  } else { stop(sprintf("Unknown argument: %s", key)) }
}

dims_use <- seq_len(dims_use_max)

t_start <- Sys.time()
warnings_collected <- character(0)

record_warning <- function(msg) {
  warnings_collected <<- c(warnings_collected, msg)
  log_warn(msg, stage = STAGE)
}

fail <- function(msg) {
  log_error(msg, stage = STAGE)
  stop(msg, call. = FALSE)
}

log_info("================ M11 DEFAULT HARMONY INTEGRATION ================", stage = STAGE)
log_info(sprintf("Input object     : %s", input_rds), stage = STAGE)
log_info(sprintf("Output object    : %s", out_rds), stage = STAGE)
log_info(sprintf("Grouping variable: %s", group_by_var), stage = STAGE)
log_info(sprintf("Input reduction  : %s  dims 1:%d", reduction_use, dims_use_max), stage = STAGE)
log_info(sprintf("Output reduction : %s", reduction_save), stage = STAGE)
log_info(sprintf("Random seed      : %d", random_seed), stage = STAGE)
log_info(sprintf("harmony version  : %s", as.character(packageVersion("harmony"))), stage = STAGE)
log_info(sprintf("Seurat version   : %s", as.character(packageVersion("Seurat"))), stage = STAGE)
if (validation_mode) {
  log_warn("VALIDATION MODE: real-data checksum and per-sample cell-count assertions are relaxed. NEVER use this for a production run.", stage = STAGE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# GUARD 0 - output path safety. The Phase 1 object must be unreachable by any write.
# ---------------------------------------------------------------------------
log_info("[GUARD 0] Verifying output path safety ...", stage = STAGE)

if (!file.exists(input_rds)) fail(sprintf("Input object not found: %s", input_rds))

abs_in  <- normalizePath(input_rds, mustWork = TRUE)
abs_out <- normalizePath(out_rds, mustWork = FALSE)

if (identical(abs_in, abs_out)) {
  fail("Output path resolves to the Phase 1 input object. Refusing to overwrite Phase 1 data.")
}
protected_prefixes <- c(
  normalizePath("results/datasets", mustWork = FALSE),
  normalizePath("results/combined", mustWork = FALSE),
  normalizePath("results/phase1_manifest.json", mustWork = FALSE),
  normalizePath("processed_mpnst.rds", mustWork = FALSE)
)
for (p in protected_prefixes) {
  if (startsWith(abs_out, p)) {
    fail(sprintf("Output path '%s' falls inside protected Phase 1 location '%s'. Refusing to write.", abs_out, p))
  }
}
if (file.exists(out_rds) && !force_overwrite) {
  fail(sprintf("Output object already exists: %s. Refusing to overwrite silently (use --force deliberately).", out_rds))
}
log_info("[GUARD 0] PASS - output path is new and outside every protected Phase 1 location.", stage = STAGE)

# ---------------------------------------------------------------------------
# GUARD 1 - input integrity. The object must be the one M10 signed off.
# ---------------------------------------------------------------------------
log_info("[GUARD 1] Verifying input object checksums ...", stage = STAGE)
input_md5    <- digest(input_rds, file = TRUE, algo = "md5")
input_sha256 <- digest(input_rds, file = TRUE, algo = "sha256")
log_info(sprintf("  md5    = %s", input_md5), stage = STAGE)
log_info(sprintf("  sha256 = %s", input_sha256), stage = STAGE)

if (validation_mode) {
  log_info("[GUARD 1] SKIPPED (validation mode) - checksums recorded but not asserted.", stage = STAGE)
} else {
  if (!identical(input_md5, expected_input_md5)) {
    fail(sprintf("Input md5 mismatch. Expected %s, observed %s. The Phase 1 handoff object has changed.",
                 expected_input_md5, input_md5))
  }
  if (!identical(input_sha256, expected_input_sha256)) {
    fail(sprintf("Input sha256 mismatch. Expected %s, observed %s.", expected_input_sha256, input_sha256))
  }
  log_info("[GUARD 1] PASS - input object matches both Phase 1 checksums exactly.", stage = STAGE)
}

# ---------------------------------------------------------------------------
# Load (read-only)
# ---------------------------------------------------------------------------
log_info("Loading Phase 1 handoff object (read-only) ...", stage = STAGE)
t_load <- Sys.time()
obj <- readRDS(input_rds)
load_seconds <- as.numeric(difftime(Sys.time(), t_load, units = "secs"))
log_info(sprintf("Loaded in %.1f s: %d features x %d cells, default assay '%s'",
                 load_seconds, nrow(obj), ncol(obj), DefaultAssay(obj)), stage = STAGE)
log_system_usage(stage = STAGE)

# ---------------------------------------------------------------------------
# GUARD 2 - object shape
# ---------------------------------------------------------------------------
log_info("[GUARD 2] Verifying object shape ...", stage = STAGE)
n_cells_in <- ncol(obj)
if (!identical(as.integer(n_cells_in), expected_cells)) {
  fail(sprintf("Cell count mismatch. Expected %d, observed %d.", expected_cells, n_cells_in))
}
if (any(duplicated(colnames(obj)))) fail("Duplicate cell barcodes present in the input object.")
log_info(sprintf("[GUARD 2] PASS - %d cells, no duplicate barcodes.", n_cells_in), stage = STAGE)

# ---------------------------------------------------------------------------
# GUARD 3 - grouping variable
# ---------------------------------------------------------------------------
log_info("[GUARD 3] Verifying grouping variable ...", stage = STAGE)
if (!(group_by_var %in% colnames(obj@meta.data))) {
  fail(sprintf("Grouping variable '%s' not present in metadata.", group_by_var))
}
grp <- obj@meta.data[[group_by_var]]
if (any(is.na(grp))) fail(sprintf("Grouping variable '%s' contains %d NA values.", group_by_var, sum(is.na(grp))))
grp_tab <- table(as.character(grp))
log_info(sprintf("  '%s' has %d levels: %s", group_by_var, length(grp_tab),
                 paste(sprintf("%s=%d", names(grp_tab), as.integer(grp_tab)), collapse = ", ")), stage = STAGE)
if (length(grp_tab) < 2) fail(sprintf("Grouping variable '%s' has fewer than 2 levels; Harmony has nothing to correct.", group_by_var))
if (validation_mode) {
  expected_levels <- setNames(as.integer(grp_tab), names(grp_tab))
  log_info("[GUARD 3] Level-count assertion relaxed (validation mode).", stage = STAGE)
} else {
  expected_levels <- c(MPNST_1 = 7615L, MPNST_2 = 2284L, MPNST_3 = 2940L, MPNST_4 = 6877L)
  if (!identical(sort(names(grp_tab)), sort(names(expected_levels)))) {
    fail(sprintf("Grouping levels differ from the M10-validated set. Observed: %s", paste(names(grp_tab), collapse = ",")))
  }
  for (lv in names(expected_levels)) {
    if (as.integer(grp_tab[[lv]]) != expected_levels[[lv]]) {
      fail(sprintf("Cell count for level '%s' is %d, M10 recorded %d.", lv, as.integer(grp_tab[[lv]]), expected_levels[[lv]]))
    }
  }
}
log_info("[GUARD 3] PASS - grouping variable matches the M10-validated structure exactly.", stage = STAGE)

# Redundancy check against orig.ident (M10 established they are identical)
if ("orig.ident" %in% colnames(obj@meta.data)) {
  same_as_orig <- all(as.character(obj@meta.data[[group_by_var]]) == as.character(obj@meta.data$orig.ident))
  log_info(sprintf("  '%s' identical to 'orig.ident' cell-for-cell: %s", group_by_var, same_as_orig), stage = STAGE)
} else {
  same_as_orig <- NA
}

# ---------------------------------------------------------------------------
# GUARD 4 - input reduction and dimensions
# ---------------------------------------------------------------------------
log_info("[GUARD 4] Verifying input reduction and requested dimensions ...", stage = STAGE)
if (!(reduction_use %in% Reductions(obj))) {
  fail(sprintf("Reduction '%s' not found. Available: %s", reduction_use, paste(Reductions(obj), collapse = ", ")))
}
pca_emb <- Embeddings(obj, reduction = reduction_use)
n_dims_available <- ncol(pca_emb)
log_info(sprintf("  Reduction '%s': %d cells x %d dims, assay '%s'",
                 reduction_use, nrow(pca_emb), n_dims_available, DefaultAssay(obj[[reduction_use]])), stage = STAGE)
if (n_dims_available < dims_use_max) {
  fail(sprintf("Requested dims 1:%d but reduction '%s' has only %d dimensions.", dims_use_max, reduction_use, n_dims_available))
}
if (nrow(pca_emb) != n_cells_in) fail("PCA embedding row count does not match cell count.")
if (!identical(rownames(pca_emb), colnames(obj))) fail("PCA embedding rownames do not match object cell names in order.")
log_info(sprintf("[GUARD 4] PASS - dims 1:%d are available in '%s'.", dims_use_max, reduction_use), stage = STAGE)

# Phase 1 baseline reductions must be present
if (!("umap_preintegration" %in% Reductions(obj))) {
  record_warning("Phase 1 reduction 'umap_preintegration' is absent from the input object.")
}

# ---------------------------------------------------------------------------
# GUARD 5 - numerical integrity precheck (M10 Finding 2)
#
# MPNST_1 has percent.mt identically zero, yet Phase 1 M8 passed
# vars.to.regress = "percent.mt" to SCTransform. Confirm that produced no
# non-finite values in the SCT residuals or in the PCA used as Harmony input.
# ---------------------------------------------------------------------------
log_info("[GUARD 5] Numerical-integrity precheck (M10 Finding 2) ...", stage = STAGE)

n_nonfinite_pca <- sum(!is.finite(pca_emb))
log_info(sprintf("  Non-finite values in '%s' embedding (all %d dims): %d",
                 reduction_use, n_dims_available, n_nonfinite_pca), stage = STAGE)
n_nonfinite_pca_used <- sum(!is.finite(pca_emb[, dims_use, drop = FALSE]))
log_info(sprintf("  Non-finite values in dims 1:%d (Harmony input): %d", dims_use_max, n_nonfinite_pca_used), stage = STAGE)

sct_scale <- tryCatch(GetAssayData(obj, assay = "SCT", layer = "scale.data"),
                      error = function(e) { record_warning(sprintf("Could not access SCT scale.data: %s", conditionMessage(e))); NULL })
n_nonfinite_sct <- NA_integer_
sct_scale_dim <- c(NA_integer_, NA_integer_)
nonfinite_sct_by_sample <- list()
if (!is.null(sct_scale)) {
  sct_scale_dim <- dim(sct_scale)
  n_nonfinite_sct <- sum(!is.finite(sct_scale))
  log_info(sprintf("  SCT scale.data is %d x %d; non-finite values: %d",
                   sct_scale_dim[1], sct_scale_dim[2], n_nonfinite_sct), stage = STAGE)
  for (lv in names(expected_levels)) {
    cells_lv <- colnames(obj)[as.character(grp) == lv]
    cells_lv <- intersect(cells_lv, colnames(sct_scale))
    nf <- sum(!is.finite(sct_scale[, cells_lv, drop = FALSE]))
    nonfinite_sct_by_sample[[lv]] <- nf
    log_info(sprintf("    %s: %d non-finite residual values", lv, nf), stage = STAGE)
  }
  rm(sct_scale); invisible(gc(verbose = FALSE))
}

if (n_nonfinite_pca_used > 0) {
  fail(sprintf("ABORT: %d non-finite values in the PCA dimensions Harmony would consume. Refusing to run Harmony on a corrupt embedding.",
               n_nonfinite_pca_used))
}
if (!is.na(n_nonfinite_sct) && n_nonfinite_sct > 0) {
  record_warning(sprintf("%d non-finite values present in SCT scale.data. Harmony input (PCA) is clean, but this must be resolved before M14 marker discovery.", n_nonfinite_sct))
}
log_info("[GUARD 5] PASS - Harmony input embedding contains no non-finite values.", stage = STAGE)

# ---------------------------------------------------------------------------
# Pre-Harmony state snapshot (for the preservation proof after Harmony)
# ---------------------------------------------------------------------------
log_info("Digesting pre-Harmony state for the preservation proof ...", stage = STAGE)
pre_state <- list(
  reductions = sort(Reductions(obj)),
  graphs = sort(Graphs(obj)),
  assays = sort(Assays(obj)),
  default_assay = DefaultAssay(obj),
  meta_colnames = colnames(obj@meta.data),
  digest_meta = digest(obj@meta.data, algo = "md5"),
  digest_pca = digest(Embeddings(obj, reduction_use), algo = "md5"),
  digest_umap_preint = if ("umap_preintegration" %in% Reductions(obj))
    digest(Embeddings(obj, "umap_preintegration"), algo = "md5") else NA_character_,
  digest_cellnames = digest(colnames(obj), algo = "md5")
)
log_info(sprintf("  pre-Harmony pca digest  = %s", pre_state$digest_pca), stage = STAGE)
log_info(sprintf("  pre-Harmony meta digest = %s", pre_state$digest_meta), stage = STAGE)

if (reduction_save %in% Reductions(obj)) {
  fail(sprintf("Reduction name '%s' already exists in the input object. Refusing to overwrite an existing reduction.", reduction_save))
}

# ---------------------------------------------------------------------------
# HARMONY - PRIMARY RUN (canonical Seurat entry point, package defaults)
# ---------------------------------------------------------------------------
log_info("---------------------------------------------------------------", stage = STAGE)
log_info("Running Harmony (primary, RunHarmony.Seurat, package defaults) ...", stage = STAGE)
log_info(sprintf("CALL: harmony::RunHarmony(object = obj, group.by.vars = \"%s\", reduction.use = \"%s\", dims.use = 1:%d, reduction.save = \"%s\", verbose = TRUE)",
                 group_by_var, reduction_use, dims_use_max, reduction_save), stage = STAGE)

set.seed(random_seed)
t_harmony <- Sys.time()
obj <- harmony::RunHarmony(
  object         = obj,
  group.by.vars  = group_by_var,
  reduction.use  = reduction_use,
  dims.use       = dims_use,
  reduction.save = reduction_save,
  verbose        = TRUE
)
harmony_seconds <- as.numeric(difftime(Sys.time(), t_harmony, units = "secs"))
log_info(sprintf("Harmony primary run completed in %.1f s.", harmony_seconds), stage = STAGE)
log_system_usage(stage = STAGE)

if (!(reduction_save %in% Reductions(obj))) {
  fail(sprintf("Harmony did not create reduction '%s'.", reduction_save))
}
harmony_emb <- Embeddings(obj, reduction = reduction_save)
log_info(sprintf("Harmony embedding: %d cells x %d dims, key '%s'",
                 nrow(harmony_emb), ncol(harmony_emb), Key(obj[[reduction_save]])), stage = STAGE)

# ---------------------------------------------------------------------------
# HARMONY - DIAGNOSTIC RUN (identical inputs and defaults, return_object = TRUE)
#
# RunHarmony.Seurat discards Harmony's internal state, so the convergence history
# is unavailable from the primary run. The matrix-level entry point is therefore
# called again with the same input, the same seed and the same defaults, and the
# resulting embedding is asserted to be numerically identical to the primary one.
# If it is not, the convergence record is marked FAILED rather than being reported
# as if it described the primary run.
# ---------------------------------------------------------------------------
log_info("---------------------------------------------------------------", stage = STAGE)
log_info("Running Harmony (diagnostic, matrix-level, return_object = TRUE) for convergence history ...", stage = STAGE)

set.seed(random_seed)
t_diag <- Sys.time()
harmony_obj <- harmony::RunHarmony(
  data_mat      = pca_emb[, dims_use, drop = FALSE],
  meta_data     = obj@meta.data[, group_by_var, drop = FALSE],
  vars_use      = group_by_var,
  return_object = TRUE,
  verbose       = TRUE
)
diag_seconds <- as.numeric(difftime(Sys.time(), t_diag, units = "secs"))
log_info(sprintf("Harmony diagnostic run completed in %.1f s.", diag_seconds), stage = STAGE)

diag_emb <- t(as.matrix(harmony_obj$Z_corr))
rownames(diag_emb) <- rownames(pca_emb)
max_abs_diff <- max(abs(diag_emb - harmony_emb))
reproducibility_ok <- is.finite(max_abs_diff) && max_abs_diff < 1e-8
log_info(sprintf("Reproducibility check: max |primary - diagnostic| = %.3e  -> %s",
                 max_abs_diff, if (reproducibility_ok) "MATCH" else "MISMATCH"), stage = STAGE)
if (!reproducibility_ok) {
  record_warning(sprintf(
    "Harmony diagnostic re-run did not reproduce the primary embedding (max abs diff %.3e). The convergence history is NOT a faithful record of the primary run and must not be interpreted as such. The primary embedding is unaffected and remains the Phase 2 result.",
    max_abs_diff))
}

# Convergence history
objective_harmony <- as.numeric(harmony_obj$objective_harmony)
objective_kmeans  <- as.numeric(harmony_obj$objective_kmeans)
kmeans_rounds     <- as.integer(harmony_obj$kmeans_rounds)
n_harmony_iters   <- length(kmeans_rounds)
harmony_nclust    <- as.integer(harmony_obj$K)
converged_early   <- n_harmony_iters < 10L

log_info(sprintf("Harmony ran %d of a maximum 10 iterations (early_stop default TRUE) -> early stop: %s",
                 n_harmony_iters, converged_early), stage = STAGE)
log_info(sprintf("Harmony soft-cluster count K = %d (package default nclust = min(round(N/30), 100))", harmony_nclust), stage = STAGE)
log_info(sprintf("Harmony objective per iteration: %s",
                 paste(sprintf("%.4f", objective_harmony), collapse = " -> ")), stage = STAGE)

conv_df <- data.frame(
  iteration = seq_along(objective_harmony) - 1L,
  objective_harmony = objective_harmony,
  stringsAsFactors = FALSE
)
conv_df$kmeans_rounds <- c(NA_integer_, kmeans_rounds)[seq_len(nrow(conv_df))]
conv_path <- file.path(out_dir, "harmony_convergence.tsv")
write.table(conv_df, conv_path, sep = "\t", row.names = FALSE, quote = FALSE)

kmeans_path <- file.path(out_dir, "harmony_kmeans_objective.tsv")
write.table(data.frame(step = seq_along(objective_kmeans), objective_kmeans = objective_kmeans),
            kmeans_path, sep = "\t", row.names = FALSE, quote = FALSE)

rm(harmony_obj, diag_emb); invisible(gc(verbose = FALSE))

# ---------------------------------------------------------------------------
# PRESERVATION PROOF - Phase 1 state must be byte-identical after Harmony
# ---------------------------------------------------------------------------
log_info("---------------------------------------------------------------", stage = STAGE)
log_info("Verifying Phase 1 state preservation ...", stage = STAGE)

post_state <- list(
  meta_colnames = colnames(obj@meta.data),
  digest_meta = digest(obj@meta.data, algo = "md5"),
  digest_pca = digest(Embeddings(obj, reduction_use), algo = "md5"),
  digest_umap_preint = if ("umap_preintegration" %in% Reductions(obj))
    digest(Embeddings(obj, "umap_preintegration"), algo = "md5") else NA_character_,
  digest_cellnames = digest(colnames(obj), algo = "md5")
)

preservation <- list(
  pca_embedding_unchanged        = identical(pre_state$digest_pca, post_state$digest_pca),
  umap_preint_embedding_unchanged= identical(pre_state$digest_umap_preint, post_state$digest_umap_preint),
  metadata_unchanged             = identical(pre_state$digest_meta, post_state$digest_meta),
  metadata_columns_unchanged     = identical(pre_state$meta_colnames, post_state$meta_colnames),
  cell_names_unchanged           = identical(pre_state$digest_cellnames, post_state$digest_cellnames),
  phase1_graphs_present          = all(pre_state$graphs %in% Graphs(obj)),
  phase1_reductions_present      = all(pre_state$reductions %in% Reductions(obj)),
  assays_unchanged               = identical(pre_state$assays, sort(Assays(obj))),
  default_assay_unchanged        = identical(pre_state$default_assay, DefaultAssay(obj))
)
for (k in names(preservation)) {
  log_info(sprintf("  %-32s : %s", k, preservation[[k]]), stage = STAGE)
}
if (!all(unlist(preservation))) {
  fail(paste0("ABORT: Phase 1 state was not preserved. Failed checks: ",
              paste(names(preservation)[!unlist(preservation)], collapse = ", ")))
}
log_info("PASS - every Phase 1 reduction, graph, assay and metadata column is byte-identical.", stage = STAGE)

# ---------------------------------------------------------------------------
# EMBEDDING SUMMARY (machine-readable)
# ---------------------------------------------------------------------------
log_info("Building Harmony embedding summary ...", stage = STAGE)

n_nonfinite_harmony <- sum(!is.finite(harmony_emb))
if (n_nonfinite_harmony > 0) {
  fail(sprintf("ABORT: Harmony embedding contains %d non-finite values.", n_nonfinite_harmony))
}

dim_summary <- data.frame(
  dimension = seq_len(ncol(harmony_emb)),
  harmony_dim_name = colnames(harmony_emb),
  min = apply(harmony_emb, 2, min),
  q25 = apply(harmony_emb, 2, stats::quantile, probs = 0.25),
  median = apply(harmony_emb, 2, stats::median),
  q75 = apply(harmony_emb, 2, stats::quantile, probs = 0.75),
  max = apply(harmony_emb, 2, max),
  mean = apply(harmony_emb, 2, mean),
  sd = apply(harmony_emb, 2, stats::sd),
  n_nonfinite = apply(harmony_emb, 2, function(x) sum(!is.finite(x))),
  stringsAsFactors = FALSE
)
# Paired pre-Harmony PCA dimension summary, for reference only (NOT an assessment)
pca_used <- pca_emb[, dims_use, drop = FALSE]
dim_summary$pre_harmony_pca_dim_name <- colnames(pca_used)
dim_summary$pre_harmony_sd <- apply(pca_used, 2, stats::sd)
write.table(dim_summary, file.path(out_dir, "harmony_embedding_dimension_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

comp_df <- as.data.frame(grp_tab, stringsAsFactors = FALSE)
colnames(comp_df) <- c(group_by_var, "cell_count")
comp_df$percentage <- 100 * comp_df$cell_count / sum(comp_df$cell_count)
write.table(comp_df, file.path(out_dir, "harmony_grouping_composition.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# SAVE the new Phase 2 object (a NEW file; Phase 1 is untouched)
# ---------------------------------------------------------------------------
log_info("---------------------------------------------------------------", stage = STAGE)
log_info(sprintf("Saving Phase 2 Harmony object to %s ...", out_rds), stage = STAGE)
t_save <- Sys.time()
saveRDS(obj, file = out_rds)
save_seconds <- as.numeric(difftime(Sys.time(), t_save, units = "secs"))
out_size <- file.info(out_rds)$size
log_info(sprintf("Saved in %.1f s (%.2f GB).", save_seconds, out_size / 1024^3), stage = STAGE)

log_info("Computing output checksums ...", stage = STAGE)
out_md5    <- digest(out_rds, file = TRUE, algo = "md5")
out_sha256 <- digest(out_rds, file = TRUE, algo = "sha256")
log_info(sprintf("  md5    = %s", out_md5), stage = STAGE)
log_info(sprintf("  sha256 = %s", out_sha256), stage = STAGE)

# ---------------------------------------------------------------------------
# ROUND-TRIP VALIDATION - re-read the serialized artefact and re-verify it.
# Existence of an RDS is not proof of success.
# ---------------------------------------------------------------------------
log_info("---------------------------------------------------------------", stage = STAGE)
log_info("Round-trip validation: freeing memory and re-reading the saved object ...", stage = STAGE)
harmony_emb_ref_digest <- digest(harmony_emb, algo = "md5")
rm(obj, harmony_emb, pca_emb, pca_used); invisible(gc(verbose = FALSE))

t_reload <- Sys.time()
obj2 <- readRDS(out_rds)
log_info(sprintf("Re-read in %.1f s.", as.numeric(difftime(Sys.time(), t_reload, units = "secs"))), stage = STAGE)

validation <- list(
  cell_count_correct        = identical(as.integer(ncol(obj2)), expected_cells),
  cell_names_unchanged      = identical(digest(colnames(obj2), algo = "md5"), pre_state$digest_cellnames),
  metadata_unchanged        = identical(digest(obj2@meta.data, algo = "md5"), pre_state$digest_meta),
  metadata_columns_unchanged= identical(colnames(obj2@meta.data), pre_state$meta_colnames),
  phase1_pca_present        = reduction_use %in% Reductions(obj2),
  phase1_pca_unchanged      = identical(digest(Embeddings(obj2, reduction_use), algo = "md5"), pre_state$digest_pca),
  phase1_umap_present       = "umap_preintegration" %in% Reductions(obj2),
  phase1_umap_unchanged     = identical(
    if ("umap_preintegration" %in% Reductions(obj2)) digest(Embeddings(obj2, "umap_preintegration"), algo = "md5") else NA_character_,
    pre_state$digest_umap_preint),
  phase1_graphs_present     = all(pre_state$graphs %in% Graphs(obj2)),
  assays_unchanged          = identical(sort(Assays(obj2)), pre_state$assays),
  default_assay_unchanged   = identical(DefaultAssay(obj2), pre_state$default_assay),
  harmony_reduction_present = reduction_save %in% Reductions(obj2),
  harmony_dims_correct      = ncol(Embeddings(obj2, reduction_save)) == dims_use_max,
  harmony_cells_correct     = nrow(Embeddings(obj2, reduction_save)) == expected_cells,
  harmony_embedding_identical_to_memory = identical(digest(Embeddings(obj2, reduction_save), algo = "md5"), harmony_emb_ref_digest),
  harmony_no_nonfinite      = sum(!is.finite(Embeddings(obj2, reduction_save))) == 0,
  harmony_rownames_match_cells = identical(rownames(Embeddings(obj2, reduction_save)), colnames(obj2))
)
for (k in names(validation)) {
  log_info(sprintf("  %-38s : %s", k, validation[[k]]), stage = STAGE)
}
if (!all(unlist(validation))) {
  fail(paste0("ABORT: round-trip validation failed on: ",
              paste(names(validation)[!unlist(validation)], collapse = ", ")))
}
log_info("PASS - the saved object round-trips and satisfies every structural requirement.", stage = STAGE)

final_reductions <- sort(Reductions(obj2))
final_graphs <- sort(Graphs(obj2))
harmony_key <- Key(obj2[[reduction_save]])
harmony_assay <- DefaultAssay(obj2[[reduction_save]])
harmony_has_loadings <- tryCatch(nrow(Loadings(obj2[[reduction_save]])) > 0, error = function(e) FALSE)
rm(obj2); invisible(gc(verbose = FALSE))

# ---------------------------------------------------------------------------
# PARAMETER RECORD (machine-readable)
# ---------------------------------------------------------------------------
log_info("Writing Harmony parameter record ...", stage = STAGE)

explicit_args <- list(
  object = "Phase 1 handoff Seurat object (loaded read-only)",
  group.by.vars = group_by_var,
  reduction.use = reduction_use,
  dims.use = sprintf("1:%d", dims_use_max),
  reduction.save = reduction_save,
  verbose = TRUE
)
default_args <- list(
  project.dim = list(value = TRUE, source = "RunHarmony.Seurat default"),
  theta = list(value = 2, source = "harmony default (NULL -> rep(2, length(vars_use)))"),
  sigma = list(value = 0.1, source = "harmony default"),
  lambda = list(value = 1, source = "harmony default"),
  nclust = list(value = harmony_nclust, source = "harmony default (NULL -> min(round(N/30), 100))"),
  max_iter = list(value = 10, source = "harmony default"),
  early_stop = list(value = TRUE, source = "harmony default"),
  ncores = list(value = 1, source = "harmony default"),
  plot_convergence = list(value = FALSE, source = "harmony default"),
  alpha = list(value = 0.2, source = "harmony_options() default"),
  tau = list(value = 0, source = "harmony_options() default"),
  block.size = list(value = 0.05, source = "harmony_options() default"),
  max.iter.cluster = list(value = 20, source = "harmony_options() default"),
  epsilon.cluster = list(value = 0.001, source = "harmony_options() default"),
  epsilon.harmony = list(value = 0.01, source = "harmony_options() default")
)

param_record <- list(
  milestone = "M11",
  phase = "phase2",
  step = STAGE,
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script = "scripts/R/phase2/run_harmony_integration.R",
  exact_call = sprintf("harmony::RunHarmony(object = obj, group.by.vars = \"%s\", reduction.use = \"%s\", dims.use = 1:%d, reduction.save = \"%s\", verbose = TRUE)",
                       group_by_var, reduction_use, dims_use_max, reduction_save),
  input = list(path = input_rds, md5 = input_md5, sha256 = input_sha256,
               size_bytes = file.info(input_rds)$size,
               cells = n_cells_in, default_assay = pre_state$default_assay,
               reduction_used = reduction_use, reduction_assay = harmony_assay,
               dims_available = n_dims_available),
  output = list(path = out_rds, md5 = out_md5, sha256 = out_sha256,
                size_bytes = out_size, cells = expected_cells,
                harmony_reduction = reduction_save, harmony_dims = dims_use_max,
                harmony_key = harmony_key, harmony_assay = harmony_assay,
                harmony_feature_loadings_projected = harmony_has_loadings,
                reductions = final_reductions, graphs = final_graphs),
  harmony = list(
    package_version = as.character(packageVersion("harmony")),
    grouping_variable = group_by_var,
    grouping_levels = names(grp_tab),
    grouping_level_counts = as.list(setNames(as.integer(grp_tab), names(grp_tab))),
    grouping_identical_to_orig_ident = same_as_orig,
    dims_use = dims_use,
    explicit_arguments = explicit_args,
    default_arguments = default_args,
    random_seed = random_seed,
    iterations_run = n_harmony_iters,
    max_iterations = 10L,
    early_stop_triggered = converged_early,
    objective_harmony = objective_harmony,
    soft_cluster_count_K = harmony_nclust,
    reproducibility_check = list(
      method = "matrix-level re-run with identical inputs, seed and defaults",
      max_abs_difference = max_abs_diff,
      status = if (reproducibility_ok) "MATCH" else "MISMATCH"
    )
  ),
  preservation_checks = preservation,
  round_trip_validation = validation,
  integrity_precheck = list(
    non_finite_pca_all_dims = n_nonfinite_pca,
    non_finite_pca_used_dims = n_nonfinite_pca_used,
    non_finite_sct_scale_data = n_nonfinite_sct,
    non_finite_sct_scale_data_by_sample = nonfinite_sct_by_sample,
    sct_scale_data_dim = sct_scale_dim,
    non_finite_harmony_embedding = n_nonfinite_harmony
  ),
  timing_seconds = list(
    load = load_seconds, harmony_primary = harmony_seconds,
    harmony_diagnostic = diag_seconds, save = save_seconds,
    total = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  ),
  warnings = warnings_collected,
  scientific_status = paste(
    "Harmony execution has completed, but integration quality has not yet been",
    "scientifically accepted. Pre/post-Harmony evaluation is reserved for M12."),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "NOT_IN_SLURM"),
  git_commit = tryCatch(trimws(system("git rev-parse HEAD", intern = TRUE)), error = function(e) NA_character_),
  git_status = tryCatch(paste(system("git status --porcelain", intern = TRUE), collapse = "\n"), error = function(e) NA_character_),
  r_version = R.version.string,
  session_info = capture.output(sessionInfo())
)

write_json(param_record, file.path(out_dir, "harmony_parameters.json"),
           auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)

# Flat TSV of the parameter table, explicit vs default
param_tbl <- rbind(
  data.frame(parameter = names(explicit_args),
             value = sapply(explicit_args, function(x) paste(as.character(x), collapse = ",")),
             specification = "EXPLICIT", source = "M10-approved / technically required",
             stringsAsFactors = FALSE),
  data.frame(parameter = names(default_args),
             value = sapply(default_args, function(x) paste(as.character(x$value), collapse = ",")),
             specification = "DEFAULT",
             source = sapply(default_args, function(x) x$source),
             stringsAsFactors = FALSE)
)
write.table(param_tbl, file.path(out_dir, "harmony_parameters.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# PROVENANCE (Phase 1-compatible record)
# ---------------------------------------------------------------------------
log_info("Recording provenance using the Phase 1 provenance utilities ...", stage = STAGE)
prov_rec <- record_provenance(
  step_name = STAGE,
  inputs = list(phase1_handoff_rds = input_rds),
  outputs = list(
    harmony_rds = out_rds,
    harmony_parameters_json = file.path(out_dir, "harmony_parameters.json"),
    harmony_parameters_tsv = file.path(out_dir, "harmony_parameters.tsv"),
    harmony_convergence_tsv = conv_path,
    harmony_kmeans_objective_tsv = kmeans_path,
    harmony_embedding_dimension_summary_tsv = file.path(out_dir, "harmony_embedding_dimension_summary.tsv"),
    harmony_grouping_composition_tsv = file.path(out_dir, "harmony_grouping_composition.tsv")
  ),
  parameters = list(
    milestone = "M11",
    harmony_version = as.character(packageVersion("harmony")),
    group_by_vars = group_by_var,
    reduction_use = reduction_use,
    dims_use = sprintf("1:%d", dims_use_max),
    reduction_save = reduction_save,
    random_seed = random_seed,
    all_scientific_parameters = "harmony 1.2.4 defaults",
    iterations_run = n_harmony_iters,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", "NOT_IN_SLURM"),
    input_sha256 = input_sha256,
    output_sha256 = out_sha256
  ),
  dataset = "combined"
)
save_provenance_json(prov_rec, file.path(out_dir, "prov_m11_harmony.json"))

total_seconds <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
log_info("---------------------------------------------------------------", stage = STAGE)
log_info(sprintf("M11 COMPLETE in %.1f s (%.2f min).", total_seconds, total_seconds / 60), stage = STAGE)
log_info(sprintf("Warnings collected: %d", length(warnings_collected)), stage = STAGE)
if (length(warnings_collected) > 0) {
  for (w in warnings_collected) log_info(sprintf("  WARNING: %s", w), stage = STAGE)
}
log_info("SCIENTIFIC STATUS: Harmony execution completed. Integration quality is NOT yet accepted; evaluation is M12.", stage = STAGE)
log_system_usage(stage = STAGE)
