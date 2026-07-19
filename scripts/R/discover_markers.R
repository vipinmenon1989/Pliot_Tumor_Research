# scripts/R/discover_markers.R
library(Seurat)
library(jsonlite)
library(dplyr)
source("scripts/R/logging_utils.R")
source("scripts/R/provenance_utils.R")

setup_strict_logging()
log_info("Initializing Milestone 7 Marker Discovery...", stage = "discover_markers")

# Argument Parsing
args <- commandArgs(trailingOnly = TRUE)
input_file <- NULL
dataset_id <- NULL
resolution_str <- NULL
test_use <- "wilcox"
min_pct <- 0.25
logfc_threshold <- 0.25
only_pos <- TRUE
out_dir <- NULL
out_prov <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--input") { input_file <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--dataset-id") { dataset_id <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--resolution") { resolution_str <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--test-use") { test_use <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--min-pct") { min_pct <- as.numeric(args[i+1]); i <- i + 2 }
  else if (args[i] == "--logfc-threshold") { logfc_threshold <- as.numeric(args[i+1]); i <- i + 2 }
  else if (args[i] == "--only-pos") { only_pos <- as.logical(args[i+1]); i <- i + 2 }
  else if (args[i] == "--out-dir") { out_dir <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--out-prov") { out_prov <- args[i+1]; i <- i + 2 }
  else { i <- i + 1 }
}

# Validate inputs
if (is.null(input_file) || is.null(dataset_id) || is.null(resolution_str) || is.null(out_dir) || is.null(out_prov)) {
  log_error("Missing required arguments.", stage = "discover_markers")
  quit(status = 1)
}

resolution <- as.numeric(resolution_str)
log_info(sprintf("Dataset ID  : %s", dataset_id), stage = "discover_markers")
log_info(sprintf("Resolution  : %f", resolution), stage = "discover_markers")
log_info(sprintf("Input File  : %s", input_file), stage = "discover_markers")
log_info(sprintf("Output Dir  : %s", out_dir), stage = "discover_markers")

# Calculate input checksum
input_checksum <- tools::md5sum(input_file)

# Create output directory
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load Seurat object
log_info("Loading Seurat clustered object...", stage = "discover_markers")
obj <- readRDS(input_file)

# Set active identity to the requested resolution
col_name <- paste0("cluster_res_", resolution)
if (!col_name %in% colnames(obj@meta.data)) {
  log_error(sprintf("Resolution column %s not found in metadata.", col_name), stage = "discover_markers")
  quit(status = 1)
}
Idents(obj) <- obj@meta.data[[col_name]]
log_info(sprintf("Set active identity to %s. Total cells: %d. Clusters: %d", col_name, ncol(obj), length(unique(Idents(obj)))), stage = "discover_markers")

# Audit cluster sizes
cluster_sizes <- table(Idents(obj))
valid_clusters <- names(cluster_sizes[cluster_sizes >= 3])

# Output paths
markers_all_path <- file.path(out_dir, "markers_all.tsv")
markers_filtered_path <- file.path(out_dir, "markers_filtered.tsv")
top_markers_path <- file.path(out_dir, "top_markers.tsv")
summary_path <- file.path(out_dir, "marker_summary.tsv")

# Helper to write empty files in case of failure/no markers
write_empty_outputs <- function(msg) {
  log_warn(msg, stage = "discover_markers")
  
  # Empty markers dataframe with correct columns
  empty_df <- data.frame(
    p_val = numeric(), avg_log2FC = numeric(), pct.1 = numeric(), pct.2 = numeric(),
    p_val_adj = numeric(), cluster = character(), gene = character(), stringsAsFactors = FALSE
  )
  write.table(empty_df, markers_all_path, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(empty_df, markers_filtered_path, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(empty_df, top_markers_path, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Empty summary
  empty_summary <- data.frame(
    cluster = "global", cell_count = "0", marker_count = "0",
    specific_marker_count = "0", weak_support = "TRUE", small_cluster = "TRUE",
    stringsAsFactors = FALSE
  )
  write.table(empty_summary, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
}

if (length(valid_clusters) < 2) {
  write_empty_outputs("Fewer than 2 valid clusters with >= 3 cells. Skipping marker discovery.")
} else {
  # Determine active assay dynamically
  active_assay <- DefaultAssay(obj)
  log_info(sprintf("Active assay: %s", active_assay), stage = "discover_markers")
  
  if (active_assay == "SCT") {
    # PrepSCTFindMarkers is required for Seurat v5 SCT FindMarkers
    log_info("Running PrepSCTFindMarkers...", stage = "discover_markers")
    obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = FALSE)
  } else {
    log_info(sprintf("Skipping PrepSCTFindMarkers since active assay is %s", active_assay), stage = "discover_markers")
  }
  
  # Run FindAllMarkers
  log_info("Running FindAllMarkers...", stage = "discover_markers")
  markers_all <- tryCatch({
    FindAllMarkers(
      object = obj,
      assay = active_assay,
      slot = "data",
      test.use = test_use,
      min.pct = min_pct,
      logfc.threshold = logfc_threshold,
      only.pos = only_pos,
      verbose = FALSE
    )
  }, error = function(e) {
    log_error(sprintf("FindAllMarkers failed: %s", e$message), stage = "discover_markers")
    NULL
  })
  
  if (is.null(markers_all) || nrow(markers_all) == 0) {
    write_empty_outputs("No markers identified or FindAllMarkers failed.")
  } else {
    # Ensure column 'gene' is present (sometimes it is in row.names or as a column)
    if (!"gene" %in% colnames(markers_all)) {
      markers_all$gene <- rownames(markers_all)
    }
    
    # Reorder columns to match standard layout
    markers_all <- markers_all[, c("gene", "cluster", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
    
    # Save complete markers table
    write.table(markers_all, markers_all_path, sep = "\t", row.names = FALSE, quote = FALSE)
    log_info(sprintf("Saved complete markers table (%d rows) to %s", nrow(markers_all), markers_all_path), stage = "discover_markers")
    
    # Filter markers: p_val_adj < 0.05 and avg_log2FC > 0.25 (already filtered by logfc.threshold, but let's be strict on p_val_adj)
    markers_filtered <- markers_all %>%
      filter(p_val_adj < 0.05 & avg_log2FC > 0.25)
    
    write.table(markers_filtered, markers_filtered_path, sep = "\t", row.names = FALSE, quote = FALSE)
    log_info(sprintf("Saved filtered markers table (%d rows) to %s", nrow(markers_filtered), markers_filtered_path), stage = "discover_markers")
    
    # Rank and generate Top Markers Table (top 20 markers per cluster)
    # Ranking rule: sort by p_val_adj (ascending), then avg_log2FC (descending), then specificity (pct.1 - pct.2) descending
    top_markers <- markers_filtered %>%
      mutate(specificity = pct.1 - pct.2) %>%
      group_by(cluster) %>%
      arrange(p_val_adj, desc(avg_log2FC), desc(specificity), .by_group = TRUE) %>%
      slice_head(n = 20) %>%
      ungroup()
    
    write.table(top_markers, top_markers_path, sep = "\t", row.names = FALSE, quote = FALSE)
    log_info(sprintf("Saved top markers table (%d rows) to %s", nrow(top_markers), top_markers_path), stage = "discover_markers")
    
    # Calculate specificity and build summary table
    # Specific marker: a marker unique to that cluster in the filtered set
    gene_cluster_counts <- table(markers_filtered$gene)
    unique_genes <- names(gene_cluster_counts[gene_cluster_counts == 1])
    
    summary_rows <- list()
    
    # Tally marker counts per cluster
    all_clusters <- levels(Idents(obj))
    for (clust in all_clusters) {
      clust_markers <- markers_filtered %>% filter(cluster == clust)
      m_count <- nrow(clust_markers)
      spec_count <- sum(clust_markers$gene %in% unique_genes)
      c_size <- as.numeric(cluster_sizes[clust])
      if (is.na(c_size)) c_size <- 0
      
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        cluster = clust,
        cell_count = as.character(c_size),
        marker_count = as.character(m_count),
        specific_marker_count = as.character(spec_count),
        weak_support = ifelse(m_count < 5, "TRUE", "FALSE"),
        small_cluster = ifelse(c_size < 10, "TRUE", "FALSE"),
        stringsAsFactors = FALSE
      )
    }
    summary_df <- do.call(rbind, summary_rows)
    
    # Tally global metrics
    total_cells <- ncol(obj)
    total_unique_markers <- length(unique(markers_filtered$gene))
    median_markers <- median(as.numeric(summary_df$marker_count))
    weak_clusts <- summary_df$cluster[summary_df$weak_support == "TRUE"]
    small_clusts <- summary_df$cluster[summary_df$small_cluster == "TRUE"]
    
    global_row <- data.frame(
      cluster = "global",
      cell_count = as.character(total_cells),
      marker_count = as.character(total_unique_markers),
      specific_marker_count = as.character(round(median_markers, 2)),
      weak_support = ifelse(length(weak_clusts) > 0, paste(weak_clusts, collapse = ","), "None"),
      small_cluster = ifelse(length(small_clusts) > 0, paste(small_clusts, collapse = ","), "None"),
      stringsAsFactors = FALSE
    )
    
    final_summary <- rbind(global_row, summary_df)
    write.table(final_summary, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
    log_info(sprintf("Saved marker summary to %s", summary_path), stage = "discover_markers")
  }
}

# Record and save execution provenance
inputs <- list(clustered_rds = input_file)
outputs <- list(
  markers_all = markers_all_path,
  markers_filtered = markers_filtered_path,
  top_markers = top_markers_path,
  summary = summary_path
)
parameters <- list(
  dataset_id = dataset_id,
  resolution = resolution,
  test_use = test_use,
  min_pct = min_pct,
  logfc_threshold = logfc_threshold,
  only_pos = only_pos,
  input_md5 = input_checksum
)

prov_record <- record_provenance("discover_markers", inputs, outputs, parameters, dataset = dataset_id)
save_provenance_json(prov_record, out_prov)
log_info("Marker Discovery completed successfully.", stage = "discover_markers")
