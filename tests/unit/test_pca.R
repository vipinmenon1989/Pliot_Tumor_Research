# tests/unit/test_pca.R
library(jsonlite)
library(Seurat)
library(yaml)
source("scripts/R/logging_utils.R")

setup_strict_logging()
log_info("Starting PCA and PC evaluation validation tests...", stage = "test_pca")

# Argument Parsing
args <- commandArgs(trailingOnly = TRUE)
datasets <- NULL
i <- 1
while (i <= length(args)) {
  if (args[i] == "--datasets") {
    datasets <- strsplit(args[i+1], " ")[[1]]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(datasets)) {
  log_info("No datasets argument provided. Loading config/config.yaml...", stage = "test_pca")
  config <- read_yaml("config/config.yaml")
  datasets <- config$datasets
}

log_info(sprintf("Active datasets to validate: %s", paste(datasets, collapse = ", ")), stage = "test_pca")

# Determine config file path based on dataset naming
is_test_mode <- any(grepl("sample", datasets))
config_path <- if (is_test_mode) "config/config.test.yaml" else "config/config.yaml"
log_info(sprintf("Loading config file from %s for validation...", config_path), stage = "test_pca")
config <- read_yaml(config_path)

norm_method <- config$normalization$method
active_assay <- if (norm_method == "SCTransform") "SCT" else "RNA"
npcs <- config$pca$npcs

for (ds in datasets) {
  # Skip based on test vs production mode
  if (is_test_mode && grepl("MPNST", ds)) next
  if (!is_test_mode && grepl("sample", ds)) next
  
  rds_path <- sprintf("results/datasets/%s/%s_pca.rds", ds, ds)
  report_path <- sprintf("reports/datasets/%s/PCA_REPORT.md", ds)
  prov_path <- sprintf("results/datasets/%s/pca_provenance.json", ds)
  top_genes_path <- sprintf("reports/datasets/%s/top_loading_genes.tsv", ds)
  cor_path <- sprintf("reports/datasets/%s/pc_technical_correlations.tsv", ds)
  var_path <- sprintf("reports/datasets/%s/pca_variance_explained.tsv", ds)
  fig_index_path <- sprintf("reports/datasets/%s/figure_index_m5.tsv", ds)
  
  log_info(sprintf("Validating PCA dataset: %s", ds), stage = "test_pca")
  
  # Check file presence
  files_to_check <- c(
    rds_path, report_path, prov_path, top_genes_path, cor_path, var_path, fig_index_path,
    sprintf("reports/datasets/%s/pca_elbow.pdf", ds),
    sprintf("reports/datasets/%s/pca_elbow.png", ds),
    sprintf("reports/datasets/%s/pca_cumulative_variance.pdf", ds),
    sprintf("reports/datasets/%s/pca_cumulative_variance.png", ds),
    sprintf("reports/datasets/%s/pca_loadings.pdf", ds),
    sprintf("reports/datasets/%s/pca_loadings.png", ds),
    sprintf("reports/datasets/%s/pca_heatmaps.pdf", ds),
    sprintf("reports/datasets/%s/pca_heatmaps.png", ds),
    sprintf("reports/datasets/%s/pca_correlations.pdf", ds),
    sprintf("reports/datasets/%s/pca_correlations.png", ds)
  )
  
  for (f in files_to_check) {
    if (!file.exists(f)) {
      log_error(sprintf("Missing expected output file: %s", f), stage = "test_pca")
      quit(status = 1)
    }
  }
  
  # Load RDS
  obj <- readRDS(rds_path)
  
  # Class check
  if (!inherits(obj, "Seurat")) {
    log_error(sprintf("Object at %s is not Seurat class.", rds_path), stage = "test_pca")
    quit(status = 1)
  }
  
  # Default assay check
  if (DefaultAssay(obj) != active_assay) {
    log_error(sprintf("Default assay at %s is '%s' (expected '%s').", rds_path, DefaultAssay(obj), active_assay), stage = "test_pca")
    quit(status = 1)
  }
  
  # Reductions check: must have pca
  if (!"pca" %in% names(obj@reductions)) {
    log_error(sprintf("Object at %s is missing 'pca' reduction.", rds_path), stage = "test_pca")
    quit(status = 1)
  }
  
  # PC dimensions check
  pca_reduc <- obj@reductions$pca
  embeddings <- Embeddings(pca_reduc)
  loadings <- Loadings(pca_reduc)
  
  if (ncol(embeddings) != npcs) {
    log_error(sprintf("Object PCA embeddings in %s has %d PCs (expected %d).", rds_path, ncol(embeddings), npcs), stage = "test_pca")
    quit(status = 1)
  }
  if (ncol(loadings) != npcs) {
    log_error(sprintf("Object PCA loadings in %s has %d PCs (expected %d).", rds_path, ncol(loadings), npcs), stage = "test_pca")
    quit(status = 1)
  }
  
  # Verify metadata was preserved from M4
  required_cols <- c("doublet_class", "doublet_score", "percent.ribo", "percent.mt")
  for (col in required_cols) {
    if (!col %in% colnames(obj@meta.data)) {
      log_error(sprintf("Object at %s is missing %s column in metadata.", rds_path, col), stage = "test_pca")
      quit(status = 1)
    }
  }
  
  # Load and validate provenance
  prov <- fromJSON(prov_path)
  if (prov$step != "pca_and_pc_evaluation") {
    log_error(sprintf("Provenance step name in %s is '%s' (expected 'pca_and_pc_evaluation').", prov_path, prov$step), stage = "test_pca")
    quit(status = 1)
  }
  
  log_info(sprintf("Dataset %s PCA validation PASSED.", ds), stage = "test_pca")
}

log_info("All PCA validation tests PASSED.", stage = "test_pca")
quit(status = 0)
