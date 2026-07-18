# tests/unit/test_normalization.R
library(jsonlite)
library(Seurat)
library(yaml)
source("scripts/R/logging_utils.R")

setup_strict_logging()
log_info("Starting normalization and variable features validation tests...", stage = "test_normalization")

# Detect which datasets are present in results/datasets
datasets <- list.dirs("results/datasets", full.names = FALSE, recursive = FALSE)
log_info(sprintf("Detected datasets on disk: %s", paste(datasets, collapse = ", ")), stage = "test_normalization")

if (length(datasets) == 0) {
  log_error("No datasets found to validate.", stage = "test_normalization")
  quit(status = 1)
}

# Determine if we are testing on synthetic data or real data
is_test_mode <- "sample_1" %in% datasets

config_path <- if (is_test_mode) "config/config.test.yaml" else "config/config.yaml"
log_info(sprintf("Loading config file from %s for validation...", config_path), stage = "test_normalization")
config <- read_yaml(config_path)

norm_method <- config$normalization$method
n_features <- config$normalization$n_features
active_assay <- if (norm_method == "SCTransform") "SCT" else "RNA"

for (ds in datasets) {
  # Skip based on test vs production mode
  if (is_test_mode && grepl("MPNST", ds)) next
  if (!is_test_mode && grepl("sample", ds)) next
  
  rds_path <- sprintf("results/datasets/%s/%s_normalized.rds", ds, ds)
  report_path <- sprintf("reports/datasets/%s/NORM_REPORT.md", ds)
  prov_path <- sprintf("results/datasets/%s/normalization_provenance.json", ds)
  hvf_tsv_path <- sprintf("reports/datasets/%s/variable_features.tsv", ds)
  fig_index_path <- sprintf("reports/datasets/%s/figure_index_m4.tsv", ds)
  
  log_info(sprintf("Validating normalized dataset: %s", ds), stage = "test_normalization")
  
  # Check file presence
  files_to_check <- c(
    rds_path, report_path, prov_path, hvf_tsv_path, fig_index_path,
    sprintf("reports/datasets/%s/var_features_scatter.pdf", ds),
    sprintf("reports/datasets/%s/var_features_scatter.png", ds),
    sprintf("reports/datasets/%s/var_features_distribution.pdf", ds),
    sprintf("reports/datasets/%s/var_features_distribution.png", ds),
    sprintf("reports/datasets/%s/top_features_violins.pdf", ds),
    sprintf("reports/datasets/%s/top_features_violins.png", ds)
  )
  
  for (f in files_to_check) {
    if (!file.exists(f)) {
      log_error(sprintf("Missing expected output file: %s", f), stage = "test_normalization")
      quit(status = 1)
    }
  }
  
  # Load RDS
  obj <- readRDS(rds_path)
  
  # Class check
  if (!inherits(obj, "Seurat")) {
    log_error(sprintf("Object at %s is not Seurat class.", rds_path), stage = "test_normalization")
    quit(status = 1)
  }
  
  # Default assay check
  if (DefaultAssay(obj) != active_assay) {
    log_error(sprintf("Default assay at %s is '%s' (expected '%s').", rds_path, DefaultAssay(obj), active_assay), stage = "test_normalization")
    quit(status = 1)
  }
  
  # Variable features count check
  v_features <- VariableFeatures(obj, assay = active_assay)
  if (length(v_features) != n_features) {
    # If the total number of genes in the dataset is smaller than n_features, it might select fewer.
    # Check if total genes in active assay is smaller than n_features
    total_active_genes <- nrow(obj[[active_assay]])
    expected_count <- min(n_features, total_active_genes)
    if (length(v_features) != expected_count) {
      log_error(sprintf("Object at %s contains %d variable features (expected %d).", rds_path, length(v_features), expected_count), stage = "test_normalization")
      quit(status = 1)
    }
  }
  
  # Scale data check: must be populated in active assay
  scaled_data <- GetAssayData(obj, assay = active_assay, layer = "scale.data")
  if (length(scaled_data) == 0 || is.null(scaled_data)) {
    log_error(sprintf("Object at %s has empty 'scale.data' slot in active assay '%s'.", rds_path, active_assay), stage = "test_normalization")
    quit(status = 1)
  }
  
  # Reductions check: must be empty
  if (length(obj@reductions) > 0) {
    log_error(sprintf("Object at %s has non-empty reductions slot (should be cleared/empty before PCA).", rds_path), stage = "test_normalization")
    quit(status = 1)
  }
  
  # Verify metadata was preserved from M3
  required_cols <- c("doublet_class", "doublet_score", "percent.ribo", "percent.mt")
  for (col in required_cols) {
    if (!col %in% colnames(obj@meta.data)) {
      log_error(sprintf("Object at %s is missing %s column in metadata (metadata regression).", rds_path, col), stage = "test_normalization")
      quit(status = 1)
    }
  }
  
  log_info(sprintf("Dataset %s normalization validation PASSED.", ds), stage = "test_normalization")
}

log_info("All normalization validation tests PASSED.", stage = "test_normalization")
quit(status = 0)
