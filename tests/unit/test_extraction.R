# tests/unit/test_extraction.R
library(jsonlite)
library(Seurat)
source("scripts/R/logging_utils.R")

setup_strict_logging()
log_info("Starting dataset extraction validation tests...", stage = "test_extraction")

# Detect which datasets are present in results/datasets
datasets <- list.dirs("results/datasets", full.names = FALSE, recursive = FALSE)
log_info(sprintf("Detected datasets on disk: %s", paste(datasets, collapse = ", ")), stage = "test_extraction")

if (length(datasets) == 0) {
  log_error("No extracted datasets found to validate.", stage = "test_extraction")
  quit(status = 1)
}

for (ds in datasets) {
  rds_path <- sprintf("results/datasets/%s/%s_raw.rds", ds, ds)
  val_md_path <- sprintf("reports/datasets/%s/DATASET_VALIDATION.md", ds)
  manifest_path <- sprintf("reports/datasets/%s/manifest.json", ds)
  
  log_info(sprintf("Validating dataset: %s", ds), stage = "test_extraction")
  
  # Check file presence
  if (!file.exists(rds_path)) {
    log_error(sprintf("Missing RDS: %s", rds_path), stage = "test_extraction")
    quit(status = 1)
  }
  if (!file.exists(val_md_path)) {
    log_error(sprintf("Missing validation MD: %s", val_md_path), stage = "test_extraction")
    quit(status = 1)
  }
  if (!file.exists(manifest_path)) {
    log_error(sprintf("Missing manifest JSON: %s", manifest_path), stage = "test_extraction")
    quit(status = 1)
  }
  
  # Load RDS
  obj <- readRDS(rds_path)
  
  # Class check
  if (!inherits(obj, "Seurat")) {
    log_error(sprintf("Object at %s is not Seurat class.", rds_path), stage = "test_extraction")
    quit(status = 1)
  }
  
  # Reductions check: must be empty
  if (length(obj@reductions) > 0) {
    log_error(sprintf("Object at %s has non-empty reductions slot (should be cleared).", rds_path), stage = "test_extraction")
    quit(status = 1)
  }
  
  # Metadata check
  if (!"percent.ribo" %in% colnames(obj@meta.data)) {
    log_error(sprintf("Object at %s is missing percent.ribo column.", rds_path), stage = "test_extraction")
    quit(status = 1)
  }
  
  # Counts check
  layers <- Layers(obj, assay = "RNA")
  has_counts <- any(grepl("counts", layers))
  if (!has_counts) {
    log_error(sprintf("Object at %s has no counts layer in RNA assay.", rds_path), stage = "test_extraction")
    quit(status = 1)
  }
  
  # Parse manifest
  manifest <- read_json(manifest_path, simplifyVector = TRUE)
  if (manifest$cells != ncol(obj)) {
    log_error(sprintf("Manifest cells count (%d) does not match object cells (%d)", manifest$cells, ncol(obj)), stage = "test_extraction")
    quit(status = 1)
  }
  if (manifest$genes != nrow(obj)) {
    log_error(sprintf("Manifest genes count (%d) does not match object genes (%d)", manifest$genes, nrow(obj)), stage = "test_extraction")
    quit(status = 1)
  }
  
  log_info(sprintf("Dataset %s validation PASSED.", ds), stage = "test_extraction")
}

log_info("All extraction validation tests PASSED.", stage = "test_extraction")
quit(status = 0)
