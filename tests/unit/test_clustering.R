# tests/unit/test_clustering.R
library(jsonlite)
library(Seurat)
library(yaml)
source("scripts/R/logging_utils.R")

setup_strict_logging()
log_info("Starting Clustering and Resolution Sweep validation tests...", stage = "test_clustering")

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
  log_info("No datasets argument provided. Loading config/config.yaml...", stage = "test_clustering")
  config <- read_yaml("config/config.yaml")
  datasets <- config$datasets
}

log_info(sprintf("Active datasets to validate: %s", paste(datasets, collapse = ", ")), stage = "test_clustering")

# Determine config file path based on dataset naming
is_test_mode <- any(grepl("sample", datasets))
config_path <- if (is_test_mode) "config/config.test.yaml" else "config/config.yaml"
log_info(sprintf("Loading config file from %s for validation...", config_path), stage = "test_clustering")
config <- read_yaml(config_path)

norm_method <- config$normalization$method
active_assay <- if (norm_method == "SCTransform") "SCT" else "RNA"
resolutions <- config$clustering$resolutions

for (ds in datasets) {
  # Skip based on test vs production mode
  if (is_test_mode && grepl("MPNST", ds)) next
  if (!is_test_mode && grepl("sample", ds)) next
  
  m5_rds_path <- sprintf("results/datasets/%s/%s_pca.rds", ds, ds)
  m6_rds_path <- sprintf("results/datasets/%s/%s_clustered.rds", ds, ds)
  prov_path <- sprintf("results/datasets/%s/clustering_provenance.json", ds)
  manifest_path <- sprintf("reports/datasets/%s/manifest_clustered.json", ds)
  
  sweep_stats_path <- sprintf("reports/datasets/%s/clustering_sweep_stats.tsv", ds)
  stability_metrics_path <- sprintf("reports/datasets/%s/clustering_stability_metrics.tsv", ds)
  rec_summary_path <- sprintf("reports/datasets/%s/clustering_recommendation_summary.tsv", ds)
  
  log_info(sprintf("Validating Clustered dataset: %s", ds), stage = "test_clustering")
  
  # 1. Check file presence
  files_to_check <- c(
    m5_rds_path, m6_rds_path, prov_path, manifest_path,
    sweep_stats_path, stability_metrics_path, rec_summary_path,
    sprintf("reports/datasets/%s/pca_umap_grid.pdf", ds),
    sprintf("reports/datasets/%s/pca_umap_grid.png", ds),
    sprintf("reports/datasets/%s/umap_recommended.pdf", ds),
    sprintf("reports/datasets/%s/umap_recommended.png", ds),
    sprintf("reports/datasets/%s/clustering_metrics.pdf", ds),
    sprintf("reports/datasets/%s/clustering_metrics.png", ds),
    sprintf("reports/datasets/%s/clustering_stability.pdf", ds),
    sprintf("reports/datasets/%s/clustering_stability.png", ds),
    sprintf("reports/datasets/%s/clustering_tree.pdf", ds),
    sprintf("reports/datasets/%s/clustering_tree.png", ds)
  )
  
  for (f in files_to_check) {
    if (!file.exists(f)) {
      log_error(sprintf("Missing expected output file: %s", f), stage = "test_clustering")
      quit(status = 1)
    }
  }
  
  # 2. Load Seurat objects
  log_info("Loading Seurat objects for comparison...", stage = "test_clustering")
  obj_m5 <- readRDS(m5_rds_path)
  obj_m6 <- readRDS(m6_rds_path)
  
  # 3. Class checks
  if (!inherits(obj_m5, "Seurat") || !inherits(obj_m6, "Seurat")) {
    log_error("Seurat objects fail class checks.", stage = "test_clustering")
    quit(status = 1)
  }
  
  # 4. Cell count check (M6 must not remove cells)
  if (ncol(obj_m6) != ncol(obj_m5)) {
    log_error(sprintf("Cell count mismatch. M6 count: %d. M5 count: %d. (M6 must not remove cells).", ncol(obj_m6), ncol(obj_m5)), stage = "test_clustering")
    quit(status = 1)
  }
  
  # 5. Assay check
  if (DefaultAssay(obj_m6) != active_assay) {
    log_error(sprintf("Default assay at %s is '%s' (expected '%s').", m6_rds_path, DefaultAssay(obj_m6), active_assay), stage = "test_clustering")
    quit(status = 1)
  }
  
  # 6. PCA reduction check (PCA must be preserved and identical)
  if (!"pca" %in% names(obj_m6@reductions)) {
    log_error("Clustered object is missing 'pca' reduction.", stage = "test_clustering")
    quit(status = 1)
  }
  embed_m5 <- Embeddings(obj_m5@reductions$pca)
  embed_m6 <- Embeddings(obj_m6@reductions$pca)
  if (!all(dim(embed_m5) == dim(embed_m6)) || !all(embed_m5 == embed_m6)) {
    log_error("PCA embeddings mismatch between M5 and M6 objects.", stage = "test_clustering")
    quit(status = 1)
  }
  
  # 7. SNN Neighbor graph check
  graph_name <- paste0(active_assay, "_snn")
  if (!graph_name %in% names(obj_m6@graphs)) {
    log_error(sprintf("Neighbor graph '%s' not found in M6 object.", graph_name), stage = "test_clustering")
    quit(status = 1)
  }
  
  # 8. UMAP check
  if (!"umap" %in% names(obj_m6@reductions)) {
    log_error("Clustered object is missing 'umap' reduction.", stage = "test_clustering")
    quit(status = 1)
  }
  
  # 9. Metadata Columns check
  expected_cols <- c("recommended_resolution", paste0("cluster_res_", resolutions))
  for (col in expected_cols) {
    if (!col %in% colnames(obj_m6@meta.data)) {
      log_error(sprintf("Clustered object is missing expected metadata column: %s", col), stage = "test_clustering")
      quit(status = 1)
    }
    # No NA values check
    if (any(is.na(obj_m6[[col]]))) {
      log_error(sprintf("Clustered object column '%s' contains NA values.", col), stage = "test_clustering")
      quit(status = 1)
    }
  }
  
  # 10. Recommended resolution sanity check
  rec_res <- obj_m6$recommended_resolution[1]
  if (!rec_res %in% resolutions) {
    log_error(sprintf("Recommended resolution %s is not among evaluated resolutions.", rec_res), stage = "test_clustering")
    quit(status = 1)
  }
  
  # 11. Provenance verify
  prov <- fromJSON(prov_path)
  if (prov$step != "clustering_resolution_sweep") {
    log_error("Clustering provenance step field is incorrect.", stage = "test_clustering")
    quit(status = 1)
  }
  if (prov$output_checksum != tools::md5sum(m6_rds_path)) {
    log_error("Provenance output checksum mismatch.", stage = "test_clustering")
    quit(status = 1)
  }
  
  log_info(sprintf("Dataset %s Clustering validation PASSED.", ds), stage = "test_clustering")
}

log_info("All Clustering and Resolution Sweep validation tests PASSED.", stage = "test_clustering")
quit(status = 0)
