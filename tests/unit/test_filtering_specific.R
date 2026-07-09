# tests/unit/test_filtering_specific.R
library(jsonlite)
library(Seurat)
library(yaml)
source("scripts/R/logging_utils.R")

setup_strict_logging()
log_info("Starting dataset-specific QC filtering and doublet detection validation tests...", stage = "test_filtering_specific")

# Detect which datasets are present in results/datasets
datasets <- list.dirs("results/datasets", full.names = FALSE, recursive = FALSE)
log_info(sprintf("Detected datasets on disk: %s", paste(datasets, collapse = ", ")), stage = "test_filtering_specific")

if (length(datasets) == 0) {
  log_error("No datasets found to validate.", stage = "test_filtering_specific")
  quit(status = 1)
}

# Determine if we are testing on synthetic data or real data
is_test_mode <- "sample_1" %in% datasets

config_path <- if (is_test_mode) "config/config.test.yaml" else "config/config.yaml"
log_info(sprintf("Loading config file from %s for validation...", config_path), stage = "test_filtering_specific")
config <- read_yaml(config_path)

for (ds in datasets) {
  # Skip based on test vs production mode
  if (is_test_mode && grepl("MPNST", ds)) next
  if (!is_test_mode && grepl("sample", ds)) next
  
  rds_path <- sprintf("results/datasets/%s/%s_filtered_specific.rds", ds, ds)
  report_path <- sprintf("reports/datasets/%s/FILTER_REPORT_SPECIFIC.md", ds)
  manifest_path <- sprintf("reports/datasets/%s/manifest_filtered_specific.json", ds)
  doublet_path <- sprintf("reports/datasets/%s/doublet_report_specific.json", ds)
  
  log_info(sprintf("Validating specific-filtered dataset: %s", ds), stage = "test_filtering_specific")
  
  # Check file presence
  if (!file.exists(rds_path)) {
    log_error(sprintf("Missing specific filtered RDS: %s", rds_path), stage = "test_filtering_specific")
    quit(status = 1)
  }
  if (!file.exists(report_path)) {
    log_error(sprintf("Missing filter report specific MD: %s", report_path), stage = "test_filtering_specific")
    quit(status = 1)
  }
  if (!file.exists(manifest_path)) {
    log_error(sprintf("Missing specific manifest JSON: %s", manifest_path), stage = "test_filtering_specific")
    quit(status = 1)
  }
  if (!file.exists(doublet_path)) {
    log_error(sprintf("Missing specific doublet report JSON: %s", doublet_path), stage = "test_filtering_specific")
    quit(status = 1)
  }
  
  # Load RDS
  obj <- readRDS(rds_path)
  
  # Class check
  if (!inherits(obj, "Seurat")) {
    log_error(sprintf("Object at %s is not Seurat class.", rds_path), stage = "test_filtering_specific")
    quit(status = 1)
  }
  
  # Reductions check: must be empty
  if (length(obj@reductions) > 0) {
    log_error(sprintf("Object at %s has non-empty reductions slot.", rds_path), stage = "test_filtering_specific")
    quit(status = 1)
  }
  
  # Metadata columns check
  required_cols <- c("doublet_class", "doublet_score", "percent.ribo")
  for (col in required_cols) {
    if (!col %in% colnames(obj@meta.data)) {
      log_error(sprintf("Object at %s is missing %s column in metadata.", rds_path, col), stage = "test_filtering_specific")
      quit(status = 1)
    }
  }
  
  # Retrieve dataset specific thresholds from config
  ds_config <- config$dataset_specific_qc[[ds]]
  if (is.null(ds_config)) {
    log_error(sprintf("No dataset-specific config found for %s in %s.", ds, config_path), stage = "test_filtering_specific")
    quit(status = 1)
  }
  
  min_features <- ds_config$min_features
  min_counts <- ds_config$min_counts
  max_mt <- ds_config$max_percent_mt
  max_ribo <- ds_config$max_percent_ribo
  
  # Filter compliance checks
  meta <- obj@meta.data
  
  # 1. Min Features
  if (any(meta$nFeature_RNA < min_features)) {
    log_error(sprintf("Filtered object %s contains cells with nFeature_RNA < %d", ds, min_features), stage = "test_filtering_specific")
    quit(status = 1)
  }
  # 2. Min Counts
  if (any(meta$nCount_RNA < min_counts)) {
    log_error(sprintf("Filtered object %s contains cells with nCount_RNA < %d", ds, min_counts), stage = "test_filtering_specific")
    quit(status = 1)
  }
  # 3. Max MT
  if (any(meta$percent.mt > max_mt)) {
    log_error(sprintf("Filtered object %s contains cells with percent.mt > %.1f%%", ds, max_mt), stage = "test_filtering_specific")
    quit(status = 1)
  }
  # 4. Max Ribo
  if (any(meta$percent.ribo > max_ribo)) {
    log_error(sprintf("Filtered object %s contains cells with percent.ribo > %.1f%%", ds, max_ribo), stage = "test_filtering_specific")
    quit(status = 1)
  }
  # 5. Doublets (no doublets should be retained)
  if (any(meta$doublet_class == "doublet")) {
    log_error(sprintf("Filtered object %s contains predicted doublets.", ds), stage = "test_filtering_specific")
    quit(status = 1)
  }
  
  # Parse manifest
  manifest <- read_json(manifest_path, simplifyVector = TRUE)
  if (manifest$cells_after != ncol(obj)) {
    log_error(sprintf("Manifest cells count (%d) does not match object cells (%d)", manifest$cells_after, ncol(obj)), stage = "test_filtering_specific")
    quit(status = 1)
  }
  if (manifest$genes != nrow(obj)) {
    log_error(sprintf("Manifest genes count (%d) does not match object genes (%d)", manifest$genes, nrow(obj)), stage = "test_filtering_specific")
    quit(status = 1)
  }
  
  log_info(sprintf("Dataset %s validation PASSED.", ds), stage = "test_filtering_specific")
}

log_info("All specific-filtering validation tests PASSED.", stage = "test_filtering_specific")
quit(status = 0)
