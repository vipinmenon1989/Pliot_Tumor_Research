# scripts/R/combine_pre_integration.R
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = +Inf)

# Load libraries
library(Seurat)
library(jsonlite)
library(digest)

source("scripts/R/logging_utils.R")
source("scripts/R/provenance_utils.R")

setup_strict_logging()
log_info("Starting pre-integration Seurat object combination pass...", stage = "pre_integration_combination")

# Argument Parsing
args <- commandArgs(trailingOnly = TRUE)
input_rds_files <- NULL
dataset_ids <- NULL
out_rds <- ""
out_dict <- ""
out_inv <- ""
out_prov <- ""
random_seed <- 42

i <- 1
while (i <= length(args)) {
  if (args[i] == "--inputs") {
    # Parse space-separated list of paths
    input_rds_files <- strsplit(args[i+1], " ")[[1]]
    i <- i + 2
  } else if (args[i] == "--datasets") {
    # Parse space-separated list of dataset IDs
    dataset_ids <- strsplit(args[i+1], " ")[[1]]
    i <- i + 2
  } else if (args[i] == "--out-rds") {
    out_rds <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-dict") {
    out_dict <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-inv") {
    out_inv <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-prov") {
    out_prov <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--random-seed") {
    random_seed <- as.integer(args[i+1])
    i <- i + 2
  } else {
    stop(sprintf("Unknown argument: %s", args[i]))
  }
}

if (is.null(input_rds_files) || is.null(dataset_ids) || out_rds == "" || out_dict == "" || out_inv == "" || out_prov == "") {
  stop("Missing required arguments.")
}

if (length(input_rds_files) != length(dataset_ids)) {
  stop("Number of input files does not match number of dataset IDs.")
}

set.seed(random_seed)

# Resolution mapping helper function
get_resolutions <- function(ds) {
  if (ds == "MPNST_1") {
    list(m6_rec = "0.6", m7_rec = "0.6", m7_alt = "0.3")
  } else if (ds == "MPNST_2") {
    list(m6_rec = "0.3", m7_rec = "0.3", m7_alt = "0.5")
  } else if (ds == "MPNST_3") {
    list(m6_rec = "0.6", m7_rec = "0.6", m7_alt = "0.3")
  } else if (ds == "MPNST_4") {
    list(m6_rec = "0.7", m7_rec = "0.7", m7_alt = "0.5")
  } else {
    # Synthetic samples
    list(m6_rec = "0.5", m7_rec = "0.5", m7_alt = "0.3")
  }
}

get_cluster_column <- function(obj, res_str) {
  col_name <- sprintf("cluster_res_%s", res_str)
  if (res_str == "1.0") {
    col_name <- "cluster_res_1"
  }
  if (!(col_name %in% colnames(obj@meta.data))) {
    snn_col <- sprintf("SCT_snn_res.%s", ifelse(res_str == "1.0", "1", res_str))
    if (snn_col %in% colnames(obj@meta.data)) {
      col_name <- snn_col
    } else {
      stop(sprintf("Clustering column for resolution %s not found in metadata.", res_str))
    }
  }
  return(obj[[col_name]][, 1])
}

# 1. Load and Namespace Metadata for Each Object
obj_list <- list()
expected_total_cells <- 0
input_checksums <- list()

for (idx in seq_along(input_rds_files)) {
  path <- input_rds_files[idx]
  ds <- dataset_ids[idx]
  
  log_info(sprintf("Loading clustered object for dataset %s from %s...", ds, path), stage = "pre_integration_combination")
  obj <- readRDS(path)
  expected_total_cells <- expected_total_cells + ncol(obj)
  input_checksums[[ds]] <- calculate_file_checksum(path)
  
  # Ensure dataset_id column exists
  obj$sample_id <- ds
  obj$orig.ident <- ds
  
  # Get resolution info
  res_info <- get_resolutions(ds)
  
  # Detect existing resolutions in the object metadata
  col_names <- colnames(obj@meta.data)
  res_cols <- col_names[grepl("^cluster_res_", col_names)]
  resolutions <- gsub("cluster_res_", "", res_cols)
  resolutions <- sapply(resolutions, function(r) if(r == "1") "1.0" else r)
  
  for (res in resolutions) {
    target_col <- sprintf("preint_%s_res_%s", ds, res)
    obj[[target_col]] <- get_cluster_column(obj, res)
  }
  
  # Preserve explicit milestone recommended columns
  obj$preint_M6_recommended_resolution <- res_info$m6_rec
  obj$preint_M7_recommended_resolution <- res_info$m7_rec
  obj$preint_M7_alternative_resolution <- res_info$m7_alt
  
  # Generate globally unique cluster labels
  rec_cluster <- get_cluster_column(obj, res_info$m7_rec)
  alt_cluster <- get_cluster_column(obj, res_info$m7_alt)
  
  obj$preint_recommended_cluster <- sprintf("%s_C%02d", ds, as.integer(as.character(rec_cluster)))
  obj$preint_alternative_cluster <- sprintf("%s_C%02d", ds, as.integer(as.character(alt_cluster)))
  
  # Strip unneeded reductions and graphs to keep object clean before recalculation
  for (red in Reductions(obj)) {
    obj[[red]] <- NULL
  }
  
  # Set default assay to RNA to ensure compatible starting point
  DefaultAssay(obj) <- "RNA"
  
  # Keep only RNA and SCT assays, drop any legacy integrated assays like mnn.reconstructed
  for (as in Assays(obj)) {
    if (as != "RNA" && as != "SCT") {
      obj[[as]] <- NULL
    }
  }
  
  obj_list[[ds]] <- obj
}

# 2. Merge Objects
log_info("Merging Seurat objects...", stage = "pre_integration_combination")
combined_obj <- merge(
  x = obj_list[[1]],
  y = obj_list[2:length(obj_list)],
  add.cell.ids = NULL
)

# 3. Cell and Feature Integrity Checks
log_info("Executing cell and feature integrity checks...", stage = "pre_integration_combination")
actual_combined_cells <- ncol(combined_obj)

log_info(sprintf("Expected Cell Count: %d", expected_total_cells), stage = "pre_integration_combination")
log_info(sprintf("Actual Combined Cell Count: %d", actual_combined_cells), stage = "pre_integration_combination")

if (expected_total_cells != actual_combined_cells) {
  log_error("CRITICAL ERROR: Cell count mismatch between constituent objects sum and combined object!", stage = "pre_integration_combination")
  stop("Cell count mismatch!")
}

if (any(duplicated(colnames(combined_obj)))) {
  log_error("CRITICAL ERROR: Duplicate cell names detected in combined object!", stage = "pre_integration_combination")
  stop("Duplicate cell names detected!")
}

log_info("Cell counts and cell name integrity verified successfully.", stage = "pre_integration_combination")

# Write machine-readable integrity table
integrity_summary_path <- "reports/combined/pre_integration/integrity_summary.tsv"
dir.create(dirname(integrity_summary_path), recursive = TRUE, showWarnings = FALSE)
integrity_df <- data.frame(
  check = c("total_cells_sum", "combined_cells_actual", "duplicate_cell_names", "retained_provenance"),
  value = c(expected_total_cells, actual_combined_cells, any(duplicated(colnames(combined_obj))), length(unique(combined_obj$sample_id)) == length(dataset_ids)),
  status = c("OK", "OK", "OK", "OK"),
  stringsAsFactors = FALSE
)
write.table(integrity_df, integrity_summary_path, sep = "\t", row.names = FALSE, quote = FALSE)

# 4. Generate Metadata Dictionary
log_info("Creating metadata dictionary...", stage = "pre_integration_combination")
meta_dict <- data.frame(
  metadata_field = character(),
  source_milestone = character(),
  source_dataset = character(),
  resolution_applicable = character(),
  pre_post_integration_status = character(),
  data_type = character(),
  description = character(),
  provenance = character(),
  stringsAsFactors = FALSE
)

# Populate dictionary details
standard_fields <- list(
  list(field = "orig.ident", milestone = "M1", dataset = "global", res = "N/A", status = "pre-integration", type = "character", desc = "Dataset/sample identifier", prov = "Seurat metadata"),
  list(field = "nCount_RNA", milestone = "M2", dataset = "global", res = "N/A", status = "pre-integration", type = "numeric", desc = "Total raw RNA count (library size)", prov = "QC extraction"),
  list(field = "nFeature_RNA", milestone = "M2", dataset = "global", res = "N/A", status = "pre-integration", type = "integer", desc = "Unique genes detected in RNA assay", prov = "QC extraction"),
  list(field = "sample_id", milestone = "M1", dataset = "global", res = "N/A", status = "pre-integration", type = "character", desc = "Dataset/sample identifier", prov = "Seurat metadata"),
  list(field = "percent.mt", milestone = "M2", dataset = "global", res = "N/A", status = "pre-integration", type = "numeric", desc = "Mitochondrial gene transcript percentage", prov = "QC extraction"),
  list(field = "percent.ribo", milestone = "M2", dataset = "global", res = "N/A", status = "pre-integration", type = "numeric", desc = "Ribosomal gene transcript percentage", prov = "QC extraction"),
  list(field = "doublet_class", milestone = "M3", dataset = "global", res = "N/A", status = "pre-integration", type = "character", desc = "Doublet class assignment (singlet/doublet)", prov = "scDblFinder doublet call"),
  list(field = "doublet_score", milestone = "M3", dataset = "global", res = "N/A", status = "pre-integration", type = "numeric", desc = "Doublet detection probability score", prov = "scDblFinder doublet score"),
  list(field = "preint_M6_recommended_resolution", milestone = "M6", dataset = "global", res = "N/A", status = "pre-integration", type = "character", desc = "Dataset-specific computationally recommended resolution", prov = "Clustering sweep"),
  list(field = "preint_M7_recommended_resolution", milestone = "M7", dataset = "global", res = "N/A", status = "pre-integration", type = "character", desc = "Dataset-specific finalized recommended resolution", prov = "Marker discovery review"),
  list(field = "preint_M7_alternative_resolution", milestone = "M7", dataset = "global", res = "N/A", status = "pre-integration", type = "character", desc = "Dataset-specific alternative resolution", prov = "Marker discovery review"),
  list(field = "preint_recommended_cluster", milestone = "M8", dataset = "global", res = "dataset-specific", status = "pre-integration", type = "character", desc = "Globally unique cluster label at recommended resolution", prov = "M8 metadata harmonization"),
  list(field = "preint_alternative_cluster", milestone = "M8", dataset = "global", res = "dataset-specific", status = "pre-integration", type = "character", desc = "Globally unique cluster label at alternative resolution", prov = "M8 metadata harmonization")
)

for (sf in standard_fields) {
  meta_dict[nrow(meta_dict) + 1, ] <- sf
}

# Add namespaced resolution sweep columns
for (ds in dataset_ids) {
  obj <- obj_list[[ds]]
  col_names <- colnames(obj@meta.data)
  res_cols <- col_names[grepl("^cluster_res_", col_names)]
  resolutions <- gsub("cluster_res_", "", res_cols)
  resolutions <- sapply(resolutions, function(r) if(r == "1") "1.0" else r)
  
  for (res in resolutions) {
    meta_dict[nrow(meta_dict) + 1, ] <- list(
      metadata_field = sprintf("preint_%s_res_%s", ds, res),
      source_milestone = "M6",
      source_dataset = ds,
      resolution_applicable = res,
      pre_post_integration_status = "pre-integration",
      data_type = "factor",
      description = sprintf("Clustering partition for %s at resolution %s", ds, res),
      provenance = "Independent M6 resolution sweep"
    )
  }
}

write.table(meta_dict, out_dict, sep = "\t", row.names = FALSE, quote = FALSE)
log_info(sprintf("Saved metadata dictionary to %s", out_dict), stage = "pre_integration_combination")

# 5. Generate Metadata Inventory
log_info("Creating metadata inventory...", stage = "pre_integration_combination")
inv_rows <- list()
all_fields <- colnames(combined_obj@meta.data)

for (f in all_fields) {
  vals <- combined_obj[[f]][, 1]
  type_str <- class(vals)[1]
  num_uniq <- length(unique(vals))
  uniq_str <- ""
  if (num_uniq <= 20) {
    uniq_str <- paste(sort(as.character(unique(vals))), collapse = ";")
  } else {
    uniq_str <- "cardinality > 20"
  }
  
  na_count <- sum(is.na(vals))
  na_pct <- (na_count / length(vals)) * 100
  
  classification <- "technical"
  if (grepl("preint_", f) || f == "orig.ident" || f == "sample_id") {
    classification <- "provenance"
  } else if (f == "orig.anno") {
    classification = "biological"
  }
  
  src <- "M1-M7"
  if (grepl("preint_recommended", f) || grepl("preint_alternative", f)) {
    src <- "M8"
  }
  
  used_vis <- "no"
  reason_vis <- "Not selected for pre-integration baseline visualizations."
  if (f %in% c("sample_id", "orig.ident", "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "preint_recommended_cluster")) {
    used_vis <- "yes"
    reason_vis <- "Required pre-integration baseline metadata visualization."
  }
  
  inv_rows[[f]] <- data.frame(
    metadata_field = f,
    type = type_str,
    num_unique_values = num_uniq,
    unique_values = uniq_str,
    missing_count = na_count,
    missing_percentage = na_pct,
    classification = classification,
    source = src,
    used_in_visualization = used_vis,
    reason_if_not_used = reason_vis,
    stringsAsFactors = FALSE
  )
}

meta_inv_df <- do.call(rbind, inv_rows)
write.table(meta_inv_df, out_inv, sep = "\t", row.names = FALSE, quote = FALSE)
log_info(sprintf("Saved metadata inventory to %s", out_inv), stage = "pre_integration_combination")

# 6. Save Combined Temporary Object
log_info(sprintf("Saving merged pre-integration Seurat object to %s...", out_rds), stage = "pre_integration_combination")
saveRDS(combined_obj, file = out_rds)

# 7. Record Provenance
log_info("Recording execution provenance...", stage = "pre_integration_combination")
prov_inputs <- list()
for (idx in seq_along(input_rds_files)) {
  prov_inputs[[sprintf("input_rds_%d", idx)]] <- input_rds_files[idx]
}
prov_outputs <- list(
  merged_rds = out_rds,
  dict_tsv = out_dict,
  inv_tsv = out_inv,
  integrity_tsv = integrity_summary_path
)
prov_params <- list(
  datasets = dataset_ids,
  expected_cells = expected_total_cells,
  actual_cells = actual_combined_cells,
  checksums = input_checksums,
  random_seed = random_seed
)
prov_rec <- record_provenance(
  step_name = "pre_integration_combination",
  inputs = prov_inputs,
  outputs = prov_outputs,
  parameters = prov_params,
  dataset = "combined"
)
save_provenance_json(prov_rec, out_prov)

log_info("Pre-integration combination completed successfully.", stage = "pre_integration_combination")
