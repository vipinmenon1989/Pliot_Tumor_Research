# scripts/R/visualize_markers.R
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
source("scripts/R/logging_utils.R")
source("scripts/R/provenance_utils.R")

setup_strict_logging()
log_info("Initializing Milestone 7 Marker Visualization...", stage = "visualize_markers")

# Argument Parsing
args <- commandArgs(trailingOnly = TRUE)
input_file <- NULL
dataset_id <- NULL
resolution_str <- NULL
markers_file <- NULL
out_dir <- NULL
out_prov <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--input") { input_file <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--dataset-id") { dataset_id <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--resolution") { resolution_str <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--markers") { markers_file <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--out-dir") { out_dir <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--out-prov") { out_prov <- args[i+1]; i <- i + 2 }
  else { i <- i + 1 }
}

# Validate inputs
if (is.null(input_file) || is.null(dataset_id) || is.null(resolution_str) || is.null(markers_file) || is.null(out_dir) || is.null(out_prov)) {
  log_error("Missing required arguments.", stage = "visualize_markers")
  quit(status = 1)
}

resolution <- as.numeric(resolution_str)
log_info(sprintf("Dataset ID  : %s", dataset_id), stage = "visualize_markers")
log_info(sprintf("Resolution  : %f", resolution), stage = "visualize_markers")
log_info(sprintf("Input File  : %s", input_file), stage = "visualize_markers")
log_info(sprintf("Markers File: %s", markers_file), stage = "visualize_markers")
log_info(sprintf("Output Dir  : %s", out_dir), stage = "visualize_markers")

# Calculate input checksums
input_checksum <- tools::md5sum(input_file)
markers_checksum <- tools::md5sum(markers_file)

# Create output directories
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Load Seurat object
log_info("Loading Seurat clustered object...", stage = "visualize_markers")
obj <- readRDS(input_file)

# Set active identity to the requested resolution
col_name <- paste0("cluster_res_", resolution)
if (!col_name %in% colnames(obj@meta.data)) {
  log_error(sprintf("Resolution column %s not found in metadata.", col_name), stage = "visualize_markers")
  quit(status = 1)
}
Idents(obj) <- obj@meta.data[[col_name]]

# Load markers and select top 5 per cluster
if (!file.exists(markers_file)) {
  log_error(sprintf("Markers file %s does not exist.", markers_file), stage = "visualize_markers")
  quit(status = 1)
}
markers_all <- read.table(markers_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

if (nrow(markers_all) == 0) {
  log_warn("No markers to visualize. Creating dummy output.", stage = "visualize_markers")
  # Write empty visualized tsv
  empty_df <- data.frame(gene = character(), cluster = character(), rank = numeric(), stringsAsFactors = FALSE)
  write.table(empty_df, file.path(out_dir, "visualized_top5_markers.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save empty provenance
  inputs <- list(clustered_rds = input_file, markers_tsv = markers_file)
  outputs <- list(visualized_tsv = file.path(out_dir, "visualized_top5_markers.tsv"))
  parameters <- list(dataset_id = dataset_id, resolution = resolution)
  prov_record <- record_provenance("visualize_markers", inputs, outputs, parameters, dataset = dataset_id)
  save_provenance_json(prov_record, out_prov)
  quit(status = 0)
}

# Apply ranking rule to identify top 5 markers per cluster
# Sort: p_val_adj (asc), then avg_log2FC (desc), then specificity (pct.1 - pct.2) desc
top5_markers <- markers_all %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
  mutate(specificity = pct.1 - pct.2) %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), desc(specificity), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Save visualized genes list
visualized_tsv_path <- file.path(out_dir, "visualized_top5_markers.tsv")
write.table(top5_markers[, c("gene", "cluster", "rank", "p_val_adj", "avg_log2FC", "pct.1", "pct.2")], 
            visualized_tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)
log_info(sprintf("Saved visualized genes table to %s", visualized_tsv_path), stage = "visualize_markers")

# Retrieve ordered unique genes
unique_genes <- unique(top5_markers$gene)
log_info(sprintf("Top 5 markers count: %d genes across %d clusters", length(unique_genes), length(unique(top5_markers$cluster))), stage = "visualize_markers")

# Detect active assay dynamically
active_assay <- DefaultAssay(obj)
log_info(sprintf("Active assay: %s", active_assay), stage = "visualize_markers")

if (length(unique_genes) > 0) {
  # 1. MARKER HEATMAP
  # Downsample Seurat object to max 100 cells per cluster to keep heatmap readable using base R
  log_info("Generating Marker Heatmap...", stage = "visualize_markers")
  set.seed(42)
  cell_groups <- split(colnames(obj), obj@meta.data[[col_name]])
  cells_to_keep <- unlist(lapply(cell_groups, function(x) {
    sample(x, min(100, length(x)))
  }))
  
  # Clear graphs in a copy of the object to prevent Seurat v5 subset graph errors
  obj_sub <- obj
  obj_sub@graphs <- list()
  obj_sub <- subset(obj_sub, cells = cells_to_keep)
  
  # Scale data for the top genes to ensure DoHeatmap works correctly
  obj_sub <- ScaleData(obj_sub, features = unique_genes, assay = active_assay, verbose = FALSE)
  
  p_heatmap <- DoHeatmap(obj_sub, features = unique_genes, assay = active_assay, slot = "scale.data", size = 3) +
    theme(axis.text.y = element_text(size = 6))
  
  heatmap_prefix <- file.path(fig_dir, sprintf("%s_res%s_top5_heatmap", dataset_id, resolution_str))
  ggsave(paste0(heatmap_prefix, ".png"), plot = p_heatmap, width = 12, height = 10, dpi = 150)
  ggsave(paste0(heatmap_prefix, ".pdf"), plot = p_heatmap, width = 12, height = 10)
  
  # 2. MARKER DOT PLOT
  log_info("Generating Marker Dot Plot...", stage = "visualize_markers")
  p_dot <- DotPlot(obj, features = unique_genes, assay = active_assay) +
    RotatedAxis() +
    theme(axis.text.x = element_text(size = 7)) +
    labs(title = sprintf("Top 5 Cluster Markers for %s at Resolution %s", dataset_id, resolution_str))
  
  dot_prefix <- file.path(fig_dir, sprintf("%s_res%s_top5_dotplot", dataset_id, resolution_str))
  ggsave(paste0(dot_prefix, ".png"), plot = p_dot, width = 12, height = 8, dpi = 150)
  ggsave(paste0(dot_prefix, ".pdf"), plot = p_dot, width = 12, height = 8)
  
  # 3. FEATURE PLOTS PER CLUSTER
  # Loop over clusters and generate a 2x3 multi-panel figure per cluster containing its top 5 markers
  log_info("Generating Feature Plots per cluster...", stage = "visualize_markers")
  DefaultAssay(obj) <- active_assay
  all_clusters <- sort(unique(top5_markers$cluster))
  for (clust in all_clusters) {
    clust_genes <- top5_markers %>% filter(cluster == clust) %>% arrange(rank) %>% pull(gene)
    
    plots <- list()
    for (gene in clust_genes) {
      plots[[gene]] <- FeaturePlot(obj, features = gene, reduction = "umap") +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 10, face = "bold"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right"
        )
    }
    
    # Fill remaining panels of the 2x3 grid (up to 6 panels) with an empty plot or UMAP overview
    if (length(plots) < 6) {
      # Use basic DimPlot as panel 6 to show the cluster's location on the UMAP
      # Highlight the current cluster in red, others in grey
      cells_highlight <- colnames(obj)[Idents(obj) == clust]
      plots[["Overview"]] <- DimPlot(obj, cells.highlight = cells_highlight, cols.highlight = "#d62728", cols = "grey90", sizes.highlight = 1) +
        labs(title = sprintf("Cluster %s Location", clust)) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 10, face = "bold"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none"
        )
    }
    
    cluster_grid <- wrap_plots(plots, ncol = 3) +
      plot_annotation(
        title = sprintf("%s Resolution %s - Cluster %s Top Markers", dataset_id, resolution_str, clust),
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      )
    
    # Cluster ID formatted as 2-digit matching user example
    cluster_num <- as.integer(clust)
    cluster_formatted <- sprintf("%02d", cluster_num)
    
    # Name format: MPNST_1_res0.6_cluster03_top5_featureplots.pdf
    feat_prefix <- file.path(fig_dir, sprintf("%s_res%s_cluster%s_top5_featureplots", dataset_id, resolution_str, cluster_formatted))
    ggsave(paste0(feat_prefix, ".png"), plot = cluster_grid, width = 11, height = 7, dpi = 150)
    ggsave(paste0(feat_prefix, ".pdf"), plot = cluster_grid, width = 11, height = 7)
  }
}

# Record and save execution provenance
inputs <- list(clustered_rds = input_file, markers_tsv = markers_file)
outputs <- list(
  visualized_tsv = visualized_tsv_path,
  heatmap_png = paste0(heatmap_prefix, ".png"),
  dotplot_png = paste0(dot_prefix, ".png")
)
parameters <- list(
  dataset_id = dataset_id,
  resolution = resolution,
  input_md5 = input_checksum,
  markers_md5 = markers_checksum
)

prov_record <- record_provenance("visualize_markers", inputs, outputs, parameters, dataset = dataset_id)
save_provenance_json(prov_record, out_prov)
log_info("Marker Visualization completed successfully.", stage = "visualize_markers")
