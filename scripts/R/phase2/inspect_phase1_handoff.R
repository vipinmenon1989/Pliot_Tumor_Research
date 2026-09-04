# scripts/R/phase2/inspect_phase1_handoff.R
#
# Phase 2 / Milestone M10 - READ-ONLY structural validation of the Phase 1 handoff object.
#
# This script does NOT run Harmony, does NOT cluster, does NOT modify any Phase 1
# artifact. It loads the frozen Phase 1 combined pre-integration object inside a
# SLURM allocation and records the exact structural facts required to design the
# Phase 2 Harmony integration script (M11):
#   - file checksums (md5 to match Phase 1 provenance, sha256 to match phase1_manifest.json)
#   - assays, layers per assay, default assay
#   - number of SCT models (decides whether PrepSCTFindMarkers() is required in M14)
#   - dimensional reductions and their dimensionality
#   - graphs present
#   - metadata columns, types, cardinality
#   - candidate Harmony grouping-variable cross-tabulations
#   - legacy annotation cross-tabulation (PROVENANCE ONLY - never used to pick parameters)
#
# Outputs are small JSON/TSV files that can then be inspected on the login node.

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = +Inf)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(jsonlite)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- ""
out_dir <- ""

i <- 1
while (i <= length(args)) {
  if (args[i] == "--input") {
    input_rds <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--out-dir") {
    out_dir <- args[i + 1]; i <- i + 2
  } else {
    stop(sprintf("Unknown argument: %s", args[i]))
  }
}

if (input_rds == "" || out_dir == "") {
  stop("Usage: inspect_phase1_handoff.R --input <rds> --out-dir <dir>")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...)))
  flush.console()
}

t_start <- Sys.time()

# ---------------------------------------------------------------------------
# 1. Checksums of the frozen input (computed on the compute node, not login)
# ---------------------------------------------------------------------------
log_msg("Computing checksums for %s ...", input_rds)
file_size_bytes <- file.info(input_rds)$size
md5_sum <- digest(input_rds, file = TRUE, algo = "md5")
sha256_sum <- digest(input_rds, file = TRUE, algo = "sha256")
log_msg("md5    = %s", md5_sum)
log_msg("sha256 = %s", sha256_sum)

# ---------------------------------------------------------------------------
# 2. Load object (read-only; never re-saved by this script)
# ---------------------------------------------------------------------------
log_msg("Loading Seurat object ...")
obj <- readRDS(input_rds)
log_msg("Loaded: %d features x %d cells", nrow(obj), ncol(obj))

# ---------------------------------------------------------------------------
# 3. Assay / layer structure
# ---------------------------------------------------------------------------
assay_records <- list()
for (a in Assays(obj)) {
  as_obj <- obj[[a]]
  layers_present <- tryCatch(SeuratObject::Layers(as_obj), error = function(e) NA_character_)
  n_models <- NA_integer_
  model_names <- NA_character_
  if (inherits(as_obj, "SCTAssay")) {
    model_names_vec <- tryCatch(levels(as_obj), error = function(e) character(0))
    n_models <- length(model_names_vec)
    model_names <- paste(model_names_vec, collapse = ";")
  }
  assay_records[[a]] <- list(
    assay = a,
    class = class(as_obj)[1],
    is_default = identical(a, DefaultAssay(obj)),
    n_features = nrow(as_obj),
    n_cells = ncol(as_obj),
    layers = paste(layers_present, collapse = ";"),
    n_variable_features = length(tryCatch(VariableFeatures(obj, assay = a), error = function(e) character(0))),
    n_sct_models = n_models,
    sct_model_names = model_names
  )
}

# ---------------------------------------------------------------------------
# 4. Reductions and graphs
# ---------------------------------------------------------------------------
reduction_records <- list()
for (r in Reductions(obj)) {
  emb <- Embeddings(obj, reduction = r)
  reduction_records[[r]] <- list(
    reduction = r,
    key = Key(obj[[r]]),
    n_cells = nrow(emb),
    n_dims = ncol(emb),
    assay_used = tryCatch(DefaultAssay(obj[[r]]), error = function(e) NA_character_),
    has_feature_loadings = tryCatch(nrow(Loadings(obj[[r]])) > 0, error = function(e) FALSE),
    has_stdev = tryCatch(length(Stdev(obj, reduction = r)) > 0, error = function(e) FALSE)
  )
}
graph_names <- Graphs(obj)

# ---------------------------------------------------------------------------
# 5. Metadata inventory
# ---------------------------------------------------------------------------
md <- obj@meta.data
meta_rows <- lapply(colnames(md), function(f) {
  v <- md[[f]]
  n_uniq <- length(unique(v[!is.na(v)]))
  uniq_str <- if (n_uniq <= 20) paste(sort(unique(as.character(v[!is.na(v)]))), collapse = ";") else "cardinality > 20"
  data.frame(
    metadata_field = f,
    type = class(v)[1],
    n_unique = n_uniq,
    unique_values = uniq_str,
    n_missing = sum(is.na(v)),
    pct_missing = 100 * sum(is.na(v)) / length(v),
    stringsAsFactors = FALSE
  )
})
meta_df <- do.call(rbind, meta_rows)
write.table(meta_df, file.path(out_dir, "handoff_metadata_inventory.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# 6. Candidate Harmony grouping-variable structure
# ---------------------------------------------------------------------------
candidate_fields <- intersect(c("sample_id", "orig.ident"), colnames(md))
xtab_path <- file.path(out_dir, "candidate_batch_crosstab.tsv")
if (length(candidate_fields) == 2) {
  xt <- as.data.frame(table(md[["sample_id"]], md[["orig.ident"]]), stringsAsFactors = FALSE)
  colnames(xt) <- c("sample_id", "orig.ident", "cell_count")
  write.table(xt, xtab_path, sep = "\t", row.names = FALSE, quote = FALSE)
  identical_mapping <- all(as.character(md[["sample_id"]]) == as.character(md[["orig.ident"]]))
} else {
  identical_mapping <- NA
}

# Cells per candidate batch level
cells_per_sample <- as.data.frame(table(md[["sample_id"]]), stringsAsFactors = FALSE)
colnames(cells_per_sample) <- c("sample_id", "cell_count")
write.table(cells_per_sample, file.path(out_dir, "cells_per_sample.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# 7. Legacy annotation cross-tabulation -- PROVENANCE ONLY.
#    Recorded so the researcher can see what the legacy analysis claimed.
#    It is NOT used to choose any Phase 2 parameter, and must not be.
# ---------------------------------------------------------------------------
if ("orig.anno" %in% colnames(md)) {
  legacy <- as.data.frame(table(md[["orig.anno"]], md[["sample_id"]]), stringsAsFactors = FALSE)
  colnames(legacy) <- c("legacy_orig_anno", "sample_id", "cell_count")
  write.table(legacy, file.path(out_dir, "legacy_orig_anno_by_sample_PROVENANCE_ONLY.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# ---------------------------------------------------------------------------
# 8. QC covariate summary per candidate batch level
# ---------------------------------------------------------------------------
qc_fields <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), colnames(md))
qc_rows <- list()
for (s in sort(unique(as.character(md[["sample_id"]])))) {
  sel <- as.character(md[["sample_id"]]) == s
  for (q in qc_fields) {
    qc_rows[[paste(s, q)]] <- data.frame(
      sample_id = s, qc_metric = q,
      n = sum(sel),
      mean = mean(md[[q]][sel]), median = median(md[[q]][sel]),
      sd = sd(md[[q]][sel]),
      min = min(md[[q]][sel]), max = max(md[[q]][sel]),
      stringsAsFactors = FALSE
    )
  }
}
write.table(do.call(rbind, qc_rows), file.path(out_dir, "qc_summary_by_sample.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# 9. Structural summary JSON
# ---------------------------------------------------------------------------
summary_rec <- list(
  script = "scripts/R/phase2/inspect_phase1_handoff.R",
  milestone = "M10",
  phase = "phase2",
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  input = list(
    path = input_rds,
    size_bytes = file_size_bytes,
    md5 = md5_sum,
    sha256 = sha256_sum
  ),
  object = list(
    n_features_default_assay = nrow(obj),
    n_cells = ncol(obj),
    default_assay = DefaultAssay(obj),
    assays = unname(assay_records),
    reductions = unname(reduction_records),
    graphs = graph_names,
    n_metadata_columns = ncol(md),
    metadata_columns = colnames(md)
  ),
  candidate_batch_variables = list(
    fields = candidate_fields,
    sample_id_equals_orig_ident = identical_mapping,
    cells_per_sample = setNames(as.list(cells_per_sample$cell_count), cells_per_sample$sample_id)
  ),
  seurat_v5_marker_prerequisites = list(
    sct_model_count = assay_records[["SCT"]]$n_sct_models,
    rna_layers = assay_records[["RNA"]]$layers,
    note = paste(
      "If SCT model count > 1, PrepSCTFindMarkers() is required before FindAllMarkers()",
      "on the SCT assay. If RNA layers are split (counts.<sample>), JoinLayers() is required",
      "before RNA-assay differential expression."
    )
  ),
  git_commit = tryCatch(trimws(system("git rev-parse HEAD", intern = TRUE)), error = function(e) NA_character_),
  git_status = tryCatch(paste(system("git status --porcelain", intern = TRUE), collapse = "\n"), error = function(e) NA_character_),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "NOT_IN_SLURM"),
  elapsed_seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
  session_info = capture.output(sessionInfo())
)

write_json(summary_rec, file.path(out_dir, "handoff_structure.json"),
           auto_unbox = TRUE, pretty = TRUE, null = "null")

# Flat TSV of assays/reductions for quick login-node inspection
assay_df <- do.call(rbind, lapply(assay_records, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
write.table(assay_df, file.path(out_dir, "handoff_assays.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
red_df <- do.call(rbind, lapply(reduction_records, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
write.table(red_df, file.path(out_dir, "handoff_reductions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

log_msg("Inspection complete in %.1f seconds. Outputs written to %s",
        as.numeric(difftime(Sys.time(), t_start, units = "secs")), out_dir)

# Explicitly confirm nothing was written back over the input
stopifnot(!identical(normalizePath(input_rds), normalizePath(out_dir, mustWork = FALSE)))
log_msg("Phase 1 input object was opened READ-ONLY and was not modified.")
