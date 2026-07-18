# scripts/R/run_clustering_sweep.R
library(Seurat)
library(yaml)
library(jsonlite)
library(ggplot2)
library(patchwork)
library(mclust)
source("scripts/R/logging_utils.R")
source("scripts/R/provenance_utils.R")

setup_strict_logging()
log_info("Initializing Milestone 6 Clustering Sweep...", stage = "run_clustering_sweep")

# Argument Parsing
args <- commandArgs(trailingOnly = TRUE)
input_file <- NULL
dataset_id <- NULL
out_rds <- NULL
out_plot_dir <- NULL
resolutions_str <- NULL
random_seed <- 42
algorithm <- 1 # Louvain
k_param <- 20
n_resamples <- 5
subsample_prop <- 0.8
pcs_rec_file <- "reports/PCA_RECOMMENDATIONS.tsv"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--input") { input_file <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--dataset-id") { dataset_id <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--out-rds") { out_rds <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--out-plot-dir") { out_plot_dir <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--resolutions") { resolutions_str <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--random-seed") { random_seed <- as.integer(args[i+1]); i <- i + 2 }
  else if (args[i] == "--algorithm") { algorithm <- as.integer(args[i+1]); i <- i + 2 }
  else if (args[i] == "--k-param") { k_param <- as.integer(args[i+1]); i <- i + 2 }
  else if (args[i] == "--n-resamples") { n_resamples <- as.integer(args[i+1]); i <- i + 2 }
  else if (args[i] == "--subsample-prop") { subsample_prop <- as.numeric(args[i+1]); i <- i + 2 }
  else if (args[i] == "--pcs-recommendations-file") { pcs_rec_file <- args[i+1]; i <- i + 2 }
  else { i <- i + 1 }
}

# Validate inputs
if (is.null(input_file) || is.null(dataset_id) || is.null(out_rds) || is.null(out_plot_dir) || is.null(resolutions_str)) {
  log_error("Missing required arguments.", stage = "run_clustering_sweep")
  quit(status = 1)
}

# Parse resolutions
resolutions <- as.numeric(strsplit(resolutions_str, ",")[[1]])
log_info(sprintf("Dataset ID  : %s", dataset_id), stage = "run_clustering_sweep")
log_info(sprintf("Input File  : %s", input_file), stage = "run_clustering_sweep")
log_info(sprintf("Output RDS  : %s", out_rds), stage = "run_clustering_sweep")
log_info(sprintf("Resolutions : %s", paste(resolutions, collapse = ", ")), stage = "run_clustering_sweep")

# Load recommendations to retrieve PCs
pcs_used <- NULL
if (file.exists(pcs_rec_file)) {
  rec_df <- read.table(pcs_rec_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  row_idx <- which(rec_df$Dataset == dataset_id)
  if (length(row_idx) > 0) {
    pcs_used <- as.integer(rec_df$Recommended_PCs[row_idx])
    log_info(sprintf("Retrieved recommended PCs from %s: 1:%d", pcs_rec_file, pcs_used), stage = "run_clustering_sweep")
  }
}

if (is.null(pcs_used)) {
  if (grepl("sample", dataset_id)) {
    pcs_used <- 10 # default for synthetic test
    log_warn("Recommended PCs file not found. Using default test PC count: 10.", stage = "run_clustering_sweep")
  } else {
    pcs_used <- 8 # default fallback
    log_warn("Recommended PCs file not found. Using default production PC count: 8.", stage = "run_clustering_sweep")
  }
}

# Calculate input checksum
input_checksum <- tools::md5sum(input_file)

# Load immutable Seurat object
log_info("Loading Seurat PCA object...", stage = "run_clustering_sweep")
obj <- readRDS(input_file)
active_assay <- DefaultAssay(obj)
log_info(sprintf("Loaded Seurat object. Assay: %s. Cells: %d. Features: %d", active_assay, ncol(obj), nrow(obj)), stage = "run_clustering_sweep")

# Graph Construction
log_info("Constructing SNN graph...", stage = "run_clustering_sweep")
# Ensure any legacy graphs are cleared
obj@graphs <- list()
obj <- FindNeighbors(
  obj,
  dims = 1:pcs_used,
  reduction = "pca",
  k.param = k_param,
  annoy.metric = "euclidean",
  force.recalc = TRUE,
  verbose = FALSE
)
graph_name <- paste0(active_assay, "_snn")

# Compute UMAP for Visual Diagnostics (once, and reuse)
log_info("Computing UMAP embedding...", stage = "run_clustering_sweep")
obj <- RunUMAP(
  obj,
  dims = 1:pcs_used,
  reduction = "pca",
  seed.use = random_seed,
  verbose = FALSE
)

# NMI helper function in pure R
calculate_nmi <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  px <- rowSums(tab) / n
  py <- colSums(tab) / n
  pxy <- tab / n
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  mi <- 0
  for (i in 1:nrow(tab)) {
    for (j in 1:ncol(tab)) {
      if (pxy[i, j] > 0) {
        mi <- mi + pxy[i, j] * log(pxy[i, j] / (px[i] * py[j]))
      }
    }
  }
  if ((hx + hy) == 0) return(0)
  return(2 * mi / (hx + hy))
}

# Loop over resolutions
sweep_data <- list()
stability_data <- list()
covariate_data <- list()

log_info("Beginning resolution sweep and stability checks...", stage = "run_clustering_sweep")

for (res in resolutions) {
  log_info(sprintf("Running clustering at resolution %s...", res), stage = "run_clustering_sweep")
  
  # Run clustering on full dataset
  obj <- FindClusters(
    obj,
    resolution = res,
    algorithm = algorithm,
    random.seed = random_seed,
    verbose = FALSE
  )
  
  # Rename metadata column to standardized name
  col_name <- paste0("cluster_res_", res)
  obj[[col_name]] <- Idents(obj)
  
  # Gather sizes
  sizes <- table(Idents(obj))
  n_clusters <- length(sizes)
  min_size <- min(sizes)
  max_size <- max(sizes)
  median_size <- median(sizes)
  singletons <- sum(sizes == 1)
  small_cells <- sum(sizes[sizes <= 10])
  prop_small <- small_cells / ncol(obj)
  
  sweep_data[[length(sweep_data) + 1]] <- data.frame(
    Resolution = res,
    Clusters = n_clusters,
    Min_Size = min_size,
    Median_Size = median_size,
    Max_Size = max_size,
    Singletons = singletons,
    Prop_Small_Cells = prop_small,
    Algorithm = if (algorithm == 1) "Louvain" else "Leiden",
    Seed = random_seed
  )
  
  # Subsampling-based stability (ARI relative to original)
  stability_scores <- c()
  for (r in 1:n_resamples) {
    set.seed(random_seed + r)
    sub_cells <- sample(colnames(obj), size = floor(subsample_prop * ncol(obj)))
    sub_obj <- subset(obj, cells = sub_cells)
    
    # Run neighbor search and clustering on subset
    sub_obj@graphs <- list()
    sub_obj <- FindNeighbors(
      sub_obj,
      dims = 1:pcs_used,
      reduction = "pca",
      k.param = k_param,
      verbose = FALSE
    )
    sub_obj <- FindClusters(
      sub_obj,
      resolution = res,
      algorithm = algorithm,
      random.seed = random_seed,
      verbose = FALSE
    )
    
    orig_clusts <- as.character(obj[[col_name]][colnames(sub_obj), 1])
    sub_clusts <- as.character(Idents(sub_obj))
    
    ari <- adjustedRandIndex(orig_clusts, sub_clusts)
    stability_scores <- c(stability_scores, ari)
  }
  
  stability_data[[length(stability_data) + 1]] <- data.frame(
    Resolution = res,
    Mean_Stability_ARI = mean(stability_scores),
    Min_Stability_ARI = min(stability_scores),
    Max_Stability_ARI = max(stability_scores)
  )
  
  # Technical covariate R2 calculations
  metadata <- obj@meta.data
  clusters_factor <- factor(Idents(obj))
  
  if (n_clusters >= 2) {
    r2_counts <- if (sd(metadata$nCount_RNA, na.rm = TRUE) > 0) summary(lm(nCount_RNA ~ clusters_factor, data = metadata))$r.squared else 0
    r2_features <- if (sd(metadata$nFeature_RNA, na.rm = TRUE) > 0) summary(lm(nFeature_RNA ~ clusters_factor, data = metadata))$r.squared else 0
    r2_mt <- if (sd(metadata$percent.mt, na.rm = TRUE) > 0) summary(lm(percent.mt ~ clusters_factor, data = metadata))$r.squared else 0
    r2_ribo <- if (sd(metadata$percent.ribo, na.rm = TRUE) > 0) summary(lm(percent.ribo ~ clusters_factor, data = metadata))$r.squared else 0
  } else {
    r2_counts <- r2_features <- r2_mt <- r2_ribo <- 0
  }
  
  covariate_data[[length(covariate_data) + 1]] <- data.frame(
    Resolution = res,
    R2_nCount_RNA = r2_counts,
    R2_nFeature_RNA = r2_features,
    R2_percent_mt = r2_mt,
    R2_percent_ribo = r2_ribo
  )
}

sweep_df <- do.call(rbind, sweep_data)
stability_df <- do.call(rbind, stability_data)
covariate_df <- do.call(rbind, covariate_data)

# Calculate ARI/NMI between adjacent resolutions
adjacent_metrics <- list()
for (i in 1:(length(resolutions)-1)) {
  r1 <- resolutions[i]
  r2 <- resolutions[i+1]
  col1 <- paste0("cluster_res_", r1)
  col2 <- paste0("cluster_res_", r2)
  
  ari_adj <- adjustedRandIndex(as.character(obj[[col1]][,1]), as.character(obj[[col2]][,1]))
  nmi_adj <- calculate_nmi(as.character(obj[[col1]][,1]), as.character(obj[[col2]][,1]))
  
  adjacent_metrics[[length(adjacent_metrics) + 1]] <- data.frame(
    Transition = paste0(r1, "->", r2),
    ARI = ari_adj,
    NMI = nmi_adj
  )
}
adjacent_df <- do.call(rbind, adjacent_metrics)

# Run Recommendation Heuristic
# We look for a resolution in the standard biologically relevant range (0.3 to 0.8) that:
# 1. Has no singletons (minimum cluster size >= 5)
# 2. Reaches high stability (Mean_Stability_ARI >= 0.60)
# 3. Minimizes technical covariate correlation where possible
# 4. Maximizes stability within this range
standard_range <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
standard_indices <- which(resolutions %in% standard_range)
valid_indices <- standard_indices[sweep_df$Min_Size[standard_indices] >= 5 & stability_df$Mean_Stability_ARI[standard_indices] >= 0.60]

if (length(valid_indices) == 0) {
  # Fallback: expand to full resolution range with relaxed stability
  valid_indices <- which(sweep_df$Min_Size >= 3 & stability_df$Mean_Stability_ARI >= 0.50)
}

if (length(valid_indices) == 0) {
  # Direct fallback: take the one with highest stability overall
  rec_idx <- which.max(stability_df$Mean_Stability_ARI)
} else {
  # Filter for low MT covariate correlation if available
  low_mt_indices <- valid_indices[!is.na(covariate_df$R2_percent_mt[valid_indices]) & covariate_df$R2_percent_mt[valid_indices] < 0.35]
  if (length(low_mt_indices) > 0) {
    # Select the resolution with highest stability among those with low MT correlation
    rec_idx <- low_mt_indices[which.max(stability_df$Mean_Stability_ARI[low_mt_indices])]
  } else {
    # Otherwise, select the resolution with highest stability overall in the valid set
    rec_idx <- valid_indices[which.max(stability_df$Mean_Stability_ARI[valid_indices])]
  }
}

recommended_res <- resolutions[rec_idx]
log_info(sprintf("Programmatic Recommended Resolution: %s", recommended_res), stage = "run_clustering_sweep")

# Determine alternatives
cons_idx <- max(1, rec_idx - 2)
high_idx <- min(length(resolutions), rec_idx + 2)
conservative_res <- resolutions[cons_idx]
high_res <- resolutions[high_idx]

# Set recommended resolution as active identity of the final object
col_rec <- paste0("cluster_res_", recommended_res)
obj$recommended_resolution <- recommended_res
Idents(obj) <- obj[[col_rec]][, 1]

# Save Clustered Seurat Object
log_info("Saving clustered Seurat object...", stage = "run_clustering_sweep")
saveRDS(obj, out_rds)
output_checksum <- tools::md5sum(out_rds)

# Ensure plot output directory exists
dir.create(out_plot_dir, recursive = TRUE, showWarnings = FALSE)

# Figure 1: UMAP Grid colored by Resolutions
log_info("Generating diagnostic figures...", stage = "run_clustering_sweep")
plot_list <- list()
for (res in resolutions) {
  col <- paste0("cluster_res_", res)
  # Basic DimPlot
  plot_list[[as.character(res)]] <- DimPlot(obj, group.by = col, label = TRUE, label.size = 3) +
    labs(title = sprintf("Res %s", res)) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10, face = "bold"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}
umap_grid <- wrap_plots(plot_list, ncol = 3) + 
  plot_annotation(title = sprintf("UMAP Resolution Sweep for %s", dataset_id))

ggsave(file.path(out_plot_dir, "pca_umap_grid.png"), plot = umap_grid, width = 12, height = 10, dpi = 150)
ggsave(file.path(out_plot_dir, "pca_umap_grid.pdf"), plot = umap_grid, width = 12, height = 10)

# Figure 2: Recommended Resolution UMAP
recommended_umap <- DimPlot(obj, group.by = col_rec, label = TRUE, label.size = 5, repel = TRUE) +
  labs(title = sprintf("%s UMAP - Recommended Resolution %s", dataset_id, recommended_res)) +
  theme_minimal()
ggsave(file.path(out_plot_dir, "umap_recommended.png"), plot = recommended_umap, width = 8, height = 6, dpi = 150)
ggsave(file.path(out_plot_dir, "umap_recommended.pdf"), plot = recommended_umap, width = 8, height = 6)

# Figure 3: Clustering Metrics Plot
p_clus <- ggplot(sweep_df, aes(x = Resolution, y = Clusters)) +
  geom_line(color = "#1f77b4", linewidth = 1) + geom_point(color = "#1f77b4", size = 2.5) +
  theme_minimal() + labs(title = "Number of Clusters vs Resolution", x = "Resolution", y = "Clusters")

# Prepare a data frame for size distributions
dist_list <- list()
for (res in resolutions) {
  col <- paste0("cluster_res_", res)
  tbl <- table(obj[[col]][, 1])
  dist_list[[as.character(res)]] <- data.frame(
    Resolution = factor(res),
    Size = as.numeric(tbl)
  )
}
dist_df <- do.call(rbind, dist_list)

p_sizes <- ggplot(dist_df, aes(x = Resolution, y = Size)) +
  geom_boxplot(fill = "#aec7e8", color = "grey30") +
  theme_minimal() + labs(title = "Cluster Size Distribution vs Resolution", x = "Resolution", y = "Cells in Cluster")

p_stab <- ggplot(stability_df, aes(x = Resolution, y = Mean_Stability_ARI)) +
  geom_line(color = "#2ca02c", linewidth = 1) + geom_point(color = "#2ca02c", size = 2.5) +
  geom_ribbon(aes(ymin = Min_Stability_ARI, ymax = Max_Stability_ARI), fill = "#2ca02c", alpha = 0.15) +
  theme_minimal() + labs(title = "Subsampling Stability (ARI) vs Resolution", x = "Resolution", y = "Mean ARI (80% Bootstrap)")

metrics_panel <- (p_clus + p_sizes) / p_stab
ggsave(file.path(out_plot_dir, "clustering_metrics.png"), plot = metrics_panel, width = 10, height = 8, dpi = 150)
ggsave(file.path(out_plot_dir, "clustering_metrics.pdf"), plot = metrics_panel, width = 10, height = 8)

# Figure 4: Stability & Technical Covariates R2
cov_long <- data.frame(
  Resolution = rep(covariate_df$Resolution, 4),
  R2 = c(covariate_df$R2_nCount_RNA, covariate_df$R2_nFeature_RNA, covariate_df$R2_percent_mt, covariate_df$R2_percent_ribo),
  Variable = rep(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), each = nrow(covariate_df))
)
p_cov <- ggplot(cov_long, aes(x = Resolution, y = R2, color = Variable)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  theme_minimal() + labs(title = "Technical Covariate Variance Explained (R2)", x = "Resolution", y = "R2 (Covariate ~ Cluster)")

p_adj <- ggplot(adjacent_df, aes(x = Transition, y = ARI, group = 1)) +
  geom_line(color = "#d62728", linewidth = 1) + geom_point(color = "#d62728", size = 2.5) +
  theme_minimal() + labs(title = "Adjacent Resolution Similarity (ARI)", x = "Transition", y = "Adjusted Rand Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

stability_panel <- p_cov / p_adj
ggsave(file.path(out_plot_dir, "clustering_stability.png"), plot = stability_panel, width = 10, height = 8, dpi = 150)
ggsave(file.path(out_plot_dir, "clustering_stability.pdf"), plot = stability_panel, width = 10, height = 8)

# Figure 5: Custom Clustree Transition Plot
generate_transition_plot <- function(obj, resolutions) {
  nodes <- list()
  edges <- list()
  
  for (r in resolutions) {
    col <- paste0("cluster_res_", r)
    clust_sizes <- table(obj[[col]][, 1])
    for (c in names(clust_sizes)) {
      nodes[[length(nodes) + 1]] <- data.frame(
        resolution = r,
        cluster = c,
        size = as.numeric(clust_sizes[c])
      )
    }
  }
  nodes_df <- do.call(rbind, nodes)
  
  for (i in 1:(length(resolutions)-1)) {
    r1 <- resolutions[i]
    r2 <- resolutions[i+1]
    col1 <- paste0("cluster_res_", r1)
    col2 <- paste0("cluster_res_", r2)
    
    tab <- table(obj[[col1]][, 1], obj[[col2]][, 1])
    for (c1 in rownames(tab)) {
      for (c2 in colnames(tab)) {
        count <- tab[c1, c2]
        if (count > 0) {
          frac1 <- count / sum(tab[c1, ])
          if (frac1 > 0.05) { # transition threshold
            edges[[length(edges) + 1]] <- data.frame(
              r1 = r1,
              c1 = as.numeric(c1),
              r2 = r2,
              c2 = as.numeric(c2),
              weight = frac1,
              count = count
            )
          }
        }
      }
    }
  }
  edges_df <- do.call(rbind, edges)
  
  p <- ggplot()
  if (nrow(edges_df) > 0) {
    p <- p + geom_segment(
      data = edges_df,
      aes(x = r1, y = c1, xend = r2, yend = c2, alpha = weight, linewidth = count),
      color = "grey50"
    ) + scale_linewidth_continuous(range = c(0.2, 1.8), name = "Cell Count")
  }
  
  nodes_df$y <- as.numeric(as.character(nodes_df$cluster))
  p <- p + geom_point(
    data = nodes_df,
    aes(x = resolution, y = y, size = size, color = factor(cluster))
  ) +
    geom_text(
      data = nodes_df,
      aes(x = resolution, y = y, label = cluster),
      vjust = 0.5, size = 3.2, color = "black", fontface = "bold"
    ) +
    scale_size_continuous(range = c(4, 13), name = "Cluster Size") +
    theme_minimal() +
    labs(
      x = "Clustering Resolution",
      y = "Cluster ID",
      title = "Cluster Transition Tree Across Resolutions"
    ) +
    scale_x_continuous(breaks = resolutions) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    guides(color = "none")
  
  return(p)
}

transition_tree <- generate_transition_plot(obj, resolutions)
ggsave(file.path(out_plot_dir, "clustering_tree.png"), plot = transition_tree, width = 10, height = 8, dpi = 150)
ggsave(file.path(out_plot_dir, "clustering_tree.pdf"), plot = transition_tree, width = 10, height = 8)

# Write out reports and recommendations
log_info("Writing statistics and recommendations tables...", stage = "run_clustering_sweep")

# Compile sweep statistics
stats_out_df <- data.frame(
  Resolution = sweep_df$Resolution,
  Clusters = sweep_df$Clusters,
  Min_Cluster_Size = sweep_df$Min_Size,
  Median_Cluster_Size = sweep_df$Median_Size,
  Max_Cluster_Size = sweep_df$Max_Size,
  Singletons = sweep_df$Singletons,
  Prop_Small_Cells = sweep_df$Prop_Small_Cells,
  Mean_Stability_ARI = stability_df$Mean_Stability_ARI,
  R2_nCount_RNA = covariate_df$R2_nCount_RNA,
  R2_percent_mt = covariate_df$R2_percent_mt,
  R2_percent_ribo = covariate_df$R2_percent_ribo
)
write.table(stats_out_df, file = file.path(out_plot_dir, "clustering_sweep_stats.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(adjacent_df, file = file.path(out_plot_dir, "clustering_stability_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Generate dataset clustering recommendation details
rec_summary <- data.frame(
  dataset_id = dataset_id,
  cells = ncol(obj),
  pcs_used = pcs_used,
  algorithm = if (algorithm == 1) "Louvain" else "Leiden",
  recommended_resolution = recommended_res,
  n_clusters = sweep_df$Clusters[rec_idx],
  min_cluster_size = sweep_df$Min_Size[rec_idx],
  median_cluster_size = sweep_df$Median_Size[rec_idx],
  max_cluster_size = sweep_df$Max_Size[rec_idx],
  stability_metric = stability_df$Mean_Stability_ARI[rec_idx],
  technical_concern = if (!is.na(covariate_df$R2_percent_mt[rec_idx]) && covariate_df$R2_percent_mt[rec_idx] > 0.35) "MT_Bias" else "None",
  conservative_resolution = conservative_res,
  high_resolution = high_res
)
write.table(rec_summary, file = file.path(out_plot_dir, "clustering_recommendation_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Figure Index
log_info("Writing local figure index...", stage = "run_clustering_sweep")
git_commit <- tryCatch({
  trimws(system("git rev-parse HEAD", intern = TRUE))
}, error = function(e) "NO_GIT_COMMIT")

fig_index <- data.frame(
  figure_path = c(
    file.path("reports/datasets", dataset_id, "pca_umap_grid.png"),
    file.path("reports/datasets", dataset_id, "umap_recommended.png"),
    file.path("reports/datasets", dataset_id, "clustering_metrics.png"),
    file.path("reports/datasets", dataset_id, "clustering_stability.png"),
    file.path("reports/datasets", dataset_id, "clustering_tree.png")
  ),
  dataset = dataset_id,
  processing_stage = "clustering_resolution_sweep",
  analysis_method = c(
    "umap_grid",
    "umap_recommended",
    "clustering_metrics",
    "clustering_stability",
    "clustering_tree"
  ),
  parameters = sprintf("res=%s;k=%d;seed=%d", resolutions_str, k_param, random_seed),
  input_object_checksum = input_checksum,
  generating_script = "scripts/R/run_clustering_sweep.R",
  snakemake_rule = "run_clustering_sweep",
  git_commit = git_commit,
  stringsAsFactors = FALSE
)
write.table(fig_index, file = file.path(out_plot_dir, "figure_index_m6.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Provenance Log
log_info("Writing clustering provenance...", stage = "run_clustering_sweep")
prov <- list(
  dataset_id = dataset_id,
  step = "clustering_resolution_sweep",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  session_info = sessionInfo()$R.version$version.string,
  seurat_version = as.character(packageVersion("Seurat")),
  mclust_version = as.character(packageVersion("mclust")),
  random_seed = random_seed,
  pcs_used = pcs_used,
  k_param = k_param,
  algorithm = if (algorithm == 1) "Louvain" else "Leiden",
  resolutions = resolutions,
  n_resamples = n_resamples,
  subsample_prop = subsample_prop,
  input_file = input_file,
  input_checksum = input_checksum,
  output_file = out_rds,
  output_checksum = output_checksum
)
write(toJSON(prov, auto_unbox = TRUE, pretty = TRUE), file = file.path(dirname(out_rds), "clustering_provenance.json"))

# Manifest
manifest <- list(
  dataset_id = dataset_id,
  seurat_object = out_rds,
  checksum = output_checksum,
  figures = list(
    umap_grid_png = file.path(out_plot_dir, "pca_umap_grid.png"),
    umap_grid_pdf = file.path(out_plot_dir, "pca_umap_grid.pdf"),
    umap_recommended_png = file.path(out_plot_dir, "umap_recommended.png"),
    umap_recommended_pdf = file.path(out_plot_dir, "umap_recommended.pdf"),
    metrics_png = file.path(out_plot_dir, "clustering_metrics.png"),
    metrics_pdf = file.path(out_plot_dir, "clustering_metrics.pdf"),
    stability_png = file.path(out_plot_dir, "clustering_stability.png"),
    stability_pdf = file.path(out_plot_dir, "clustering_stability.pdf"),
    tree_png = file.path(out_plot_dir, "clustering_tree.png"),
    tree_pdf = file.path(out_plot_dir, "clustering_tree.pdf")
  ),
  tables = list(
    sweep_stats = file.path(out_plot_dir, "clustering_sweep_stats.tsv"),
    stability_metrics = file.path(out_plot_dir, "clustering_stability_metrics.tsv"),
    recommendation = file.path(out_plot_dir, "clustering_recommendation_summary.tsv"),
    figure_index = file.path(out_plot_dir, "figure_index_m6.tsv")
  )
)
write(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE), file = file.path(out_plot_dir, "manifest_clustered.json"))

log_info("Milestone 6 Clustering Sweep Completed successfully.", stage = "run_clustering_sweep")
quit(status = 0)
