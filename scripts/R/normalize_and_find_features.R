# scripts/R/normalize_and_find_features.R
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = +Inf)

# Load libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(jsonlite)
library(digest)

source("scripts/R/logging_utils.R")
source("scripts/R/provenance_utils.R")

setup_strict_logging()
log_info("Starting normalization and variable feature selection...", stage = "normalization")

# Argument Parsing
args <- commandArgs(trailingOnly = TRUE)
input_rds <- ""
dataset_id <- ""
norm_method <- "SCTransform"
n_features <- 3000
out_rds <- ""
out_report <- ""
out_prov <- ""
out_plot_dir <- ""
snakemake_rule <- "normalize_and_find_features"
random_seed <- 42

i <- 1
while (i <= length(args)) {
  if (args[i] == "--input") {
    input_rds <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--dataset-id") {
    dataset_id <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--method") {
    norm_method <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--n-features") {
    n_features <- as.integer(args[i+1])
    i <- i + 2
  } else if (args[i] == "--out-rds") {
    out_rds <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-report") {
    out_report <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-prov") {
    out_prov <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-plot-dir") {
    out_plot_dir <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--rule") {
    snakemake_rule <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--random-seed") {
    random_seed <- as.integer(args[i+1])
    i <- i + 2
  } else {
    stop(sprintf("Unknown argument: %s", args[i]))
  }
}

if (input_rds == "" || dataset_id == "" || out_rds == "" || out_report == "" || out_prov == "" || out_plot_dir == "") {
  stop("Missing required arguments.")
}

set.seed(random_seed)

# Ensure plot output directory exists
dir.create(out_plot_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load Seurat Object
log_info(sprintf("Loading filtered specific Seurat object from %s...", input_rds), dataset = dataset_id, stage = "normalization")
seurat_obj <- readRDS(input_rds)
total_cells <- ncol(seurat_obj)
total_genes <- nrow(seurat_obj)

# 2. Determine Regression Variables
# Dynamic MT check: if percent.mt has zero variance or is all 0, we must skip regression to prevent SCT/ScaleData failures
has_mt_variance <- FALSE
if ("percent.mt" %in% colnames(seurat_obj@meta.data)) {
  mt_vals <- seurat_obj$percent.mt
  if (length(unique(mt_vals)) > 1 && var(mt_vals) > 0) {
    has_mt_variance <- TRUE
  }
}

vars_to_regress <- NULL
if (has_mt_variance) {
  vars_to_regress <- "percent.mt"
  log_info("Mitochondrial variance detected. Regression will include: percent.mt", dataset = dataset_id, stage = "normalization")
} else {
  log_info("No mitochondrial variance detected (e.g. mapping difference). Regression will skip: percent.mt", dataset = dataset_id, stage = "normalization")
}

# 3. Perform Normalization
log_info(sprintf("Performing normalization using method: %s...", norm_method), dataset = dataset_id, stage = "normalization")
start_time <- Sys.time()

if (norm_method == "SCTransform") {
  # In Seurat v5, SCTransform uses vst.flavor = "v2" by default
  # We explicitly set seed.use for reproducibility
  seurat_obj <- SCTransform(
    seurat_obj,
    vars.to.regress = vars_to_regress,
    variable.features.n = n_features,
    verbose = FALSE,
    seed.use = random_seed
  )
  active_assay <- "SCT"
} else if (norm_method == "LogNormalize") {
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = n_features,
    verbose = FALSE
  )
  seurat_obj <- ScaleData(
    seurat_obj,
    vars.to.regress = vars_to_regress,
    verbose = FALSE
  )
  active_assay <- "RNA"
} else {
  stop(sprintf("Unsupported normalization method: %s", norm_method))
}

end_time <- Sys.time()
elapsed_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
log_info(sprintf("Normalization completed in %.2f seconds.", elapsed_sec), dataset = dataset_id, stage = "normalization")

# 4. Extract Variable Features and Statistics
log_info("Extracting highly variable feature statistics...", dataset = dataset_id, stage = "normalization")
hvf_info <- HVFInfo(seurat_obj, assay = active_assay)
hvf_info$gene <- rownames(hvf_info)

# Identify sorting metric based on assay method
sort_col <- NULL
if ("residual_variance" %in% colnames(hvf_info)) {
  sort_col <- "residual_variance"
} else if ("variance.standardized" %in% colnames(hvf_info)) {
  sort_col <- "variance.standardized"
} else if ("variance" %in% colnames(hvf_info)) {
  sort_col <- "variance"
}

if (!is.null(sort_col)) {
  hvf_info <- hvf_info[order(hvf_info[[sort_col]], decreasing = TRUE), ]
}

# Write Variable Features TSV Table
hvf_table_path <- file.path(out_plot_dir, "variable_features.tsv")
write.table(hvf_info, hvf_table_path, sep = "\t", row.names = FALSE, quote = FALSE)
log_info(sprintf("Highly variable features table saved to %s", hvf_table_path), dataset = dataset_id, stage = "normalization")

top20_genes <- head(hvf_info$gene, 20)

# 5. Generate Diagnostic Plots (Phase G)
log_info("Generating variable-feature diagnostics plots...", dataset = dataset_id, stage = "normalization")

# Premium themes
theme_premium <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 6)),
      axis.title = element_text(face = "bold", size = 9),
      axis.text = element_text(color = "black", size = 8),
      panel.grid.major = element_line(color = "gray92"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

# Plot 1: HVF Scatter Plot
p_scatter <- VariableFeaturePlot(seurat_obj, assay = active_assay)
p_scatter <- LabelPoints(plot = p_scatter, points = top20_genes, repel = TRUE, max.overlaps = 50, fontface = "bold", size = 3)
p_scatter <- p_scatter + theme_premium() + labs(title = sprintf("Highly Variable Features (%s) - %s", norm_method, dataset_id))

# Plot 2: HVF Metric Distribution (Histogram)
if (!is.null(sort_col)) {
  p_dist <- ggplot(hvf_info, aes(x = .data[[sort_col]])) +
    geom_histogram(bins = 100, fill = "#1f77b4", color = "black", alpha = 0.7) +
    geom_vline(xintercept = hvf_info[[sort_col]][min(n_features, nrow(hvf_info))], linetype = "dashed", color = "red", size = 0.8) +
    theme_premium() +
    labs(
      title = sprintf("Feature Variance Distribution - %s", dataset_id),
      x = sprintf("Metric: %s", sort_col),
      y = "Feature Count"
    )
} else {
  p_dist <- ggplot() + theme_void() + labs(title = "No sorting metric available for variance distribution")
}

# Plot 3: Top Features Violins
# Draw expression levels of top 6 variable features
top6_genes <- head(hvf_info$gene, 6)
p_violins <- VlnPlot(seurat_obj, features = top6_genes, pt.size = 0.05, combine = TRUE, assay = active_assay)
p_violins <- p_violins + plot_annotation(
  title = sprintf("Top 6 Highly Variable Features Expression - %s", dataset_id),
  theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
)

# Save plots
plot_saves <- list(
  var_features_scatter = list(plot = p_scatter, w = 8, h = 6),
  var_features_distribution = list(plot = p_dist, w = 7, h = 5),
  top_features_violins = list(plot = p_violins, w = 10, h = 7)
)

fig_index_rows <- list()
input_checksum <- calculate_file_checksum(input_rds)
git_commit <- tryCatch({
  trimws(system("git rev-parse HEAD", intern = TRUE))
}, error = function(e) "NO_GIT_COMMIT")

for (name in names(plot_saves)) {
  info <- plot_saves[[name]]
  pdf_path <- file.path(out_plot_dir, sprintf("%s.pdf", name))
  png_path <- file.path(out_plot_dir, sprintf("%s.png", name))
  
  log_info(sprintf("Saving diagnostic plot to %s and %s...", pdf_path, png_path), dataset = dataset_id, stage = "normalization")
  ggsave(pdf_path, plot = info$plot, width = info$w, height = info$h, device = "pdf")
  ggsave(png_path, plot = info$plot, width = info$w, height = info$h, dpi = 300, device = "png")
  
  fig_index_rows[[paste0(name, "_pdf")]] <- list(
    figure_path = pdf_path,
    dataset = dataset_id,
    processing_stage = "normalization_and_variable_features",
    analysis_method = name,
    parameters = sprintf("method=%s;n_features=%d;w=%d;h=%d", norm_method, n_features, info$w, info$h),
    input_object_checksum = input_checksum,
    generating_script = "scripts/R/normalize_and_find_features.R",
    snakemake_rule = snakemake_rule,
    git_commit = git_commit
  )
  fig_index_rows[[paste0(name, "_png")]] <- list(
    figure_path = png_path,
    dataset = dataset_id,
    processing_stage = "normalization_and_variable_features",
    analysis_method = name,
    parameters = sprintf("method=%s;n_features=%d;w=%d;h=%d", norm_method, n_features, info$w, info$h),
    input_object_checksum = input_checksum,
    generating_script = "scripts/R/normalize_and_find_features.R",
    snakemake_rule = snakemake_rule,
    git_commit = git_commit
  )
}

# Generate local figure index TSV
fig_index_df <- do.call(rbind, lapply(fig_index_rows, as.data.frame))
fig_index_tsv_path <- file.path(out_plot_dir, "figure_index_m4.tsv")
write.table(fig_index_df, fig_index_tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)

# 6. Save Seurat Object
log_info(sprintf("Saving normalized Seurat object to %s...", out_rds), dataset = dataset_id, stage = "normalization")
saveRDS(seurat_obj, file = out_rds)

# 7. Generate NORM_REPORT.md (Phase H)
report_content <- c(
  sprintf("# Milestone 4 Normalization and Variable Feature Report - %s", dataset_id),
  sprintf("*Generated on: %s*", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. Parameters & Configuration",
  "",
  sprintf("- **Input Filtered Object**: `%s`", input_rds),
  sprintf("- **Output Normalized Object**: `%s`", out_rds),
  sprintf("- **Normalization Method**: `%s`", norm_method),
  sprintf("- **Requested Variable Features**: `%d`", n_features),
  sprintf("- **Mitochondrial Regression**: `%s`", ifelse(has_mt_variance, "YES (percent.mt)", "NO (Zero MT variance)")),
  sprintf("- **Seed Used**: `%d`", random_seed),
  sprintf("- **Execution Time**: `%.2f seconds`", elapsed_sec),
  "",
  "## 2. Dataset Characteristics",
  "",
  sprintf("- **Total Cells**: `%d`", total_cells),
  sprintf("- **Total Raw Genes**: `%d`", total_genes),
  "",
  "## 3. High Variance Features",
  "",
  sprintf("Top 20 highly variable features selected under `%s` (sorted by %s):", norm_method, sort_col),
  "",
  paste(sprintf("%d. **%s** (%s: %.4f, mean: %.4f, variance: %.4f)", 
                1:20, top20_genes, sort_col, 
                hvf_info[[sort_col]][1:min(20, nrow(hvf_info))], hvf_info$mean[1:min(20, nrow(hvf_info))], hvf_info$variance[1:min(20, nrow(hvf_info))]), collapse = "\n"),
  "",
  "## 4. Diagnostics & Visualizations",
  "The following diagnostic plots were generated to assess normalization quality:",
  sprintf("- [Scatter Plot (HVF)](file:///%s)", file.path(out_plot_dir, "var_features_scatter.png")),
  sprintf("- [Distribution Plot (HVF Metric)](file:///%s)", file.path(out_plot_dir, "var_features_distribution.png")),
  sprintf("- [Expression Violins (Top 6)](file:///%s)", file.path(out_plot_dir, "top_features_violins.png")),
  "",
  "---",
  "## 5. Software & Environment Provenance",
  sprintf("- **Seurat Version**: `%s`", as.character(packageVersion("Seurat"))),
  sprintf("- **sctransform Version**: `%s`", as.character(packageVersion("sctransform"))),
  sprintf("- **Git Commit Hash**: `%s`", git_commit)
)

writeLines(report_content, out_report)
log_info(sprintf("Normalization report saved to %s", out_report), dataset = dataset_id, stage = "normalization")

# 8. Record Provenance JSON
prov_inputs <- list(input_rds = input_rds)
prov_outputs <- list(
  out_rds = out_rds,
  out_report = out_report,
  hvf_table = hvf_table_path,
  fig_index = fig_index_tsv_path
)
prov_params <- list(
  method = norm_method,
  n_features = n_features,
  vars_to_regress = vars_to_regress,
  random_seed = random_seed,
  elapsed_sec = elapsed_sec
)
prov_rec <- record_provenance(
  step_name = "normalization_and_variable_features",
  inputs = prov_inputs,
  outputs = prov_outputs,
  parameters = prov_params,
  dataset = dataset_id
)
save_provenance_json(prov_rec, out_prov)
log_info(sprintf("Provenance json saved to %s", out_prov), dataset = dataset_id, stage = "normalization")

log_system_usage(dataset = dataset_id, stage = "normalization")
log_info("Normalization process completed successfully.", dataset = dataset_id, stage = "normalization")
