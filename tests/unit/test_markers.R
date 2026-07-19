# tests/unit/test_markers.R
library(jsonlite)
library(yaml)
source("scripts/R/logging_utils.R")

setup_strict_logging()
log_info("Starting Marker Discovery validation tests...", stage = "test_markers")

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
  log_info("No datasets argument provided. Loading config/config.yaml...", stage = "test_markers")
  config <- read_yaml("config/config.yaml")
  datasets <- config$datasets
}

log_info(sprintf("Active datasets to validate: %s", paste(datasets, collapse = ", ")), stage = "test_markers")

is_test_mode <- any(grepl("sample", datasets))
config_path <- if (is_test_mode) "config/config.test.yaml" else "config/config.yaml"
log_info(sprintf("Loading config file from %s for validation...", config_path), stage = "test_markers")
config <- read_yaml(config_path)

resolutions <- config$clustering$resolutions
test_use <- config$markers$test_use

for (ds in datasets) {
  if (is_test_mode && grepl("MPNST", ds)) next
  if (!is_test_mode && grepl("sample", ds)) next
  
  log_info(sprintf("Validating Marker Discovery for dataset: %s", ds), stage = "test_markers")
  
  # 1. Verify files exist for every resolution
  for (res in resolutions) {
    res_formatted <- sprintf("%.1f", as.numeric(res))
    res_dir <- sprintf("reports/datasets/%s/markers/resolution_%s", ds, res_formatted)
    
    files_to_check <- c(
      file.path(res_dir, "markers_all.tsv"),
      file.path(res_dir, "markers_filtered.tsv"),
      file.path(res_dir, "top_markers.tsv"),
      file.path(res_dir, "marker_summary.tsv"),
      file.path(res_dir, "provenance.json")
    )
    
    for (f in files_to_check) {
      if (!file.exists(f)) {
        log_error(sprintf("Missing expected marker file: %s", f), stage = "test_markers")
        quit(status = 1)
      }
    }
    
    # 2. Check table formats and columns
    # markers_all
    df_all <- read.table(file.path(res_dir, "markers_all.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    expected_cols <- c("gene", "cluster", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
    if (nrow(df_all) > 0) {
      if (!all(expected_cols %in% colnames(df_all))) {
        log_error(sprintf("markers_all.tsv missing expected columns in %s", res_dir), stage = "test_markers")
        quit(status = 1)
      }
    }
    
    # markers_filtered
    df_filtered <- read.table(file.path(res_dir, "markers_filtered.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (nrow(df_filtered) > 0) {
      if (any(df_filtered$p_val_adj >= 0.05) || any(df_filtered$avg_log2FC <= 0.25)) {
        log_error(sprintf("markers_filtered.tsv contains elements violating filters (p_val_adj < 0.05, avg_log2FC > 0.25) in %s", res_dir), stage = "test_markers")
        quit(status = 1)
      }
    }
    
    # summary check
    summary_df <- read.table(file.path(res_dir, "marker_summary.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (nrow(summary_df) < 1 || summary_df$cluster[1] != "global") {
      log_error(sprintf("marker_summary.tsv has incorrect format or lacks global row in %s", res_dir), stage = "test_markers")
      quit(status = 1)
    }
    
    # provenance check
    prov <- fromJSON(file.path(res_dir, "provenance.json"))
    if (prov$step != "discover_markers" || prov$parameters$dataset_id != ds) {
      log_error(sprintf("Incorrect provenance values in %s", res_dir), stage = "test_markers")
      quit(status = 1)
    }
  }
  
  # 3. Check visualizations exist for the resolution that was actually visualized
  # Find which resolution directory under reports/datasets/ds/markers/ has a "figures" subdirectory
  markers_dir <- sprintf("reports/datasets/%s/markers", ds)
  res_dirs <- list.dirs(markers_dir, recursive = FALSE)
  vis_res_dir <- NULL
  for (rd in res_dirs) {
    if (file.exists(file.path(rd, "figures"))) {
      vis_res_dir <- rd
      break
    }
  }
  
  if (!is.null(vis_res_dir)) {
    # Extract resolution string from path
    res_formatted_name <- basename(vis_res_dir)
    res_val_str <- sub("resolution_", "", res_formatted_name)
    
    vis_dir <- file.path(vis_res_dir, "figures")
    
    heatmap_png <- file.path(vis_dir, sprintf("%s_res%s_top5_heatmap.png", ds, res_val_str))
    heatmap_pdf <- file.path(vis_dir, sprintf("%s_res%s_top5_heatmap.pdf", ds, res_val_str))
    dot_png <- file.path(vis_dir, sprintf("%s_res%s_top5_dotplot.png", ds, res_val_str))
    dot_pdf <- file.path(vis_dir, sprintf("%s_res%s_top5_dotplot.pdf", ds, res_val_str))
    vis_tsv <- file.path(vis_res_dir, "visualized_top5_markers.tsv")
      
      expected_plots <- c(heatmap_png, heatmap_pdf, dot_png, dot_pdf, vis_tsv)
      for (f in expected_plots) {
        if (!file.exists(f)) {
          log_error(sprintf("Missing recommended resolution visualization file: %s", f), stage = "test_markers")
          quit(status = 1)
        }
      }
      log_info(sprintf("Recommended resolution %s visualization files verified.", res_val_str), stage = "test_markers")
    }
  
  log_info(sprintf("Dataset %s Marker Discovery validation PASSED.", ds), stage = "test_markers")
}

log_info("All Marker Discovery validation tests PASSED.", stage = "test_markers")
quit(status = 0)
