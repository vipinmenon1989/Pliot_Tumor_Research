# Smoke Test for Configuration Parsing and Package Loading
# Verifies R_env environment integrity.

library(yaml)
library(Seurat)
source("scripts/R/logging_utils.R")

setup_strict_logging()
log_info("Starting smoke tests...", stage = "smoke_test")

test_config_loading <- function(config_path) {
  log_info(sprintf("Testing config load: %s", config_path), stage = "smoke_test")
  if (!file.exists(config_path)) {
    log_error(sprintf("Config file not found: %s", config_path), stage = "smoke_test")
    return(FALSE)
  }

  cfg <- tryCatch({
    read_yaml(config_path)
  }, error = function(e) {
    log_error(sprintf("Failed to parse config file: %s", conditionMessage(e)), stage = "smoke_test")
    return(NULL)
  })

  if (is.null(cfg)) return(FALSE)

  # Check key fields
  required_keys <- c("input_rds", "random_seed", "dataset_id_column", "qc", "normalization", "pca", "clustering", "output_dirs")
  for (k in required_keys) {
    if (!k %in% names(cfg)) {
      log_error(sprintf("Missing required configuration key: %s", k), stage = "smoke_test")
      return(FALSE)
    }
  }

  log_info(sprintf("Config file %s loaded and validated successfully.", config_path), stage = "smoke_test")
  return(TRUE)
}

test_seurat_synthetic <- function() {
  synthetic_path <- "data/synthetic/synthetic_mpnst.rds"
  log_info(sprintf("Testing synthetic Seurat object load: %s", synthetic_path), stage = "smoke_test")

  if (!file.exists(synthetic_path)) {
    log_error(sprintf("Synthetic dataset not found: %s", synthetic_path), stage = "smoke_test")
    return(FALSE)
  }

  seurat_obj <- tryCatch({
    readRDS(synthetic_path)
  }, error = function(e) {
    log_error(sprintf("Failed to load synthetic RDS: %s", conditionMessage(e)), stage = "smoke_test")
    return(NULL)
  })

  if (is.null(seurat_obj)) return(FALSE)

  # Check class and metadata
  if (!inherits(seurat_obj, "Seurat")) {
    log_error("Object is not a Seurat object class", stage = "smoke_test")
    return(FALSE)
  }

  log_info(sprintf("Successfully verified Seurat object of class: %s", as.character(class(seurat_obj))), stage = "smoke_test")
  log_info(sprintf("Cells: %d, Genes: %d", ncol(seurat_obj), nrow(seurat_obj)), stage = "smoke_test")

  # Check presence of key metadata columns
  expected_cols <- c("sample_id", "patient_id", "condition", "percent.mt", "percent.ribo")
  for (col in expected_cols) {
    if (!col %in% colnames(seurat_obj@meta.data)) {
      log_error(sprintf("Missing metadata column: %s", col), stage = "smoke_test")
      return(FALSE)
    }
  }

  log_info("Metadata columns verified successfully.", stage = "smoke_test")
  return(TRUE)
}

# Run tests
res_prod <- test_config_loading("config/config.yaml")
res_test <- test_config_loading("config/config.test.yaml")
res_synth <- test_seurat_synthetic()

if (res_prod && res_test && res_synth) {
  log_info("All smoke tests PASSED.", stage = "smoke_test")
  quit(status = 0)
} else {
  log_error("Some smoke tests FAILED.", stage = "smoke_test")
  quit(status = 1)
}
