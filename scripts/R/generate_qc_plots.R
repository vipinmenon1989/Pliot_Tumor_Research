# scripts/R/generate_qc_plots.R
options(stringsAsFactors = FALSE)
set.seed(42)

library(Seurat)
library(ggplot2)
library(patchwork)
library(jsonlite)
library(yaml)

source("scripts/R/logging_utils.R")
source("scripts/R/provenance_utils.R")

setup_strict_logging()
log_info("Starting QC plot generation...", stage = "qc_plots")

# Argument Parsing
args <- commandArgs(trailingOnly = TRUE)
input_rds <- ""
dataset_id <- ""
out_dir <- ""
out_prov <- ""
snakemake_rule <- "generate_qc_plots"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--input") {
    input_rds <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--dataset-id") {
    dataset_id <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-dir") {
    out_dir <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-prov") {
    out_prov <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--rule") {
    snakemake_rule <- args[i+1]
    i <- i + 2
  } else {
    stop(sprintf("Unknown argument: %s", args[i]))
  }
}

if (input_rds == "" || dataset_id == "" || out_dir == "" || out_prov == "") {
  stop("Missing required arguments.")
}

log_info(sprintf("Loading subsetted Seurat object from %s...", input_rds), dataset = dataset_id, stage = "qc_plots")
subset_obj <- readRDS(input_rds)

# Ensure output directory exists
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Retrieve metadata
meta <- subset_obj@meta.data

# Verify that standard columns exist
required_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")
for (col in required_cols) {
  if (!col %in% colnames(meta)) {
    stop(sprintf("Required QC column %s not found in metadata.", col))
  }
}

# Premium ggplot2 theme
theme_premium <- function() {
  theme_minimal(base_size = 12, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 8)),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30", margin = margin(b = 10)),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(color = "black"),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray92"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 15, r = 15, b = 15, l = 15)
    )
}

# Custom color palette (Sleek slate / teal HSL-like theme)
fill_color <- "#2c7fb8"
accent_color <- "#7fcdbb"
scatter_color <- "#3182bd"

# 1. Generate Violin Plots
log_info("Generating QC violin plots...", dataset = dataset_id, stage = "qc_plots")
v1 <- ggplot(meta, aes(x = factor(dataset_id), y = nCount_RNA)) +
  geom_violin(fill = fill_color, alpha = 0.7, color = "gray40") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, coef = 0) +
  theme_premium() +
  labs(x = NULL, y = "nCount_RNA (UMI)", title = "Library Depth")

v2 <- ggplot(meta, aes(x = factor(dataset_id), y = nFeature_RNA)) +
  geom_violin(fill = accent_color, alpha = 0.7, color = "gray40") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, coef = 0) +
  theme_premium() +
  labs(x = NULL, y = "nFeature_RNA (Genes)", title = "Genes Detected")

v3 <- ggplot(meta, aes(x = factor(dataset_id), y = percent.mt)) +
  geom_violin(fill = "#e34a33", alpha = 0.7, color = "gray40") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, coef = 0) +
  theme_premium() +
  labs(x = NULL, y = "percent.mt (%)", title = "Mitochondrial %")

v4 <- ggplot(meta, aes(x = factor(dataset_id), y = percent.ribo)) +
  geom_violin(fill = "#fdbb84", alpha = 0.7, color = "gray40") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, coef = 0) +
  theme_premium() +
  labs(x = NULL, y = "percent.ribo (%)", title = "Ribosomal %")

violins_plot <- (v1 | v2) / (v3 | v4) + 
  plot_annotation(
    title = sprintf("Pre-Filter QC Violin Plots - %s", dataset_id),
    subtitle = sprintf("N = %d cells", nrow(meta)),
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
  )

# 2. Generate Scatter Plot
log_info("Generating QC scatter plots...", dataset = dataset_id, stage = "qc_plots")
scatter_plot <- ggplot(meta, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_viridis_c(option = "plasma", name = "MT %") +
  theme_premium() +
  labs(
    x = "nCount_RNA (Library Size / Depth)",
    y = "nFeature_RNA (Genes Detected)",
    title = "nFeature_RNA vs nCount_RNA",
    subtitle = sprintf("Dataset: %s (N = %d cells)", dataset_id, nrow(meta))
  )

# 3. Generate Histogram Plots
log_info("Generating QC histograms...", dataset = dataset_id, stage = "qc_plots")
h1 <- ggplot(meta, aes(x = nCount_RNA)) +
  geom_histogram(bins = 50, fill = fill_color, alpha = 0.7, color = "gray40") +
  theme_premium() +
  labs(x = "nCount_RNA", y = "Cell Count", title = "Sequencing Depth")

h2 <- ggplot(meta, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 50, fill = accent_color, alpha = 0.7, color = "gray40") +
  theme_premium() +
  labs(x = "nFeature_RNA", y = "Cell Count", title = "Genes Detected")

h3 <- ggplot(meta, aes(x = percent.mt)) +
  geom_histogram(bins = 50, fill = "#e34a33", alpha = 0.7, color = "gray40") +
  theme_premium() +
  labs(x = "percent.mt (%)", y = "Cell Count", title = "Mitochondrial Fraction")

h4 <- ggplot(meta, aes(x = percent.ribo)) +
  geom_histogram(bins = 50, fill = "#fdbb84", alpha = 0.7, color = "gray40") +
  theme_premium() +
  labs(x = "percent.ribo (%)", y = "Cell Count", title = "Ribosomal Fraction")

histograms_plot <- (h1 | h2) / (h3 | h4) +
  plot_annotation(
    title = sprintf("Pre-Filter QC Histograms - %s", dataset_id),
    subtitle = sprintf("N = %d cells", nrow(meta)),
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
  )

# 4. Cell Metrics & Sample Barplot
log_info("Generating cell metrics barplot...", dataset = dataset_id, stage = "qc_plots")
metrics_df <- data.frame(
  Metric = c("Total Cells", "Median UMI", "Median Genes", "Mean MT %", "Mean Ribo %"),
  Value = c(
    nrow(meta),
    median(meta$nCount_RNA),
    median(meta$nFeature_RNA),
    mean(meta$percent.mt),
    mean(meta$percent.ribo)
  )
)

metrics_plot <- ggplot(metrics_df, aes(x = reorder(Metric, -Value), y = Value, fill = Metric)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black", show.legend = FALSE) +
  scale_fill_brewer(palette = "Set2") +
  geom_text(aes(label = ifelse(Metric == "Total Cells" | Metric == "Median UMI" | Metric == "Median Genes", 
                               format(Value, big.mark = ","), 
                               sprintf("%.2f%%", Value))), 
            vjust = -0.5, fontface = "bold") +
  theme_premium() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(
    x = NULL,
    y = "Value",
    title = "Dataset Summary Metrics",
    subtitle = sprintf("Dataset: %s", dataset_id)
  )

# Save plots as PDF and PNG
plot_saves <- list(
  violins = list(plot = violins_plot, w = 10, h = 8),
  scatter = list(plot = scatter_plot, w = 7, h = 6),
  histograms = list(plot = histograms_plot, w = 10, h = 8),
  metrics = list(plot = metrics_plot, w = 7, h = 6)
)

fig_index_rows <- list()
for (name in names(plot_saves)) {
  info <- plot_saves[[name]]
  
  pdf_path <- file.path(out_dir, sprintf("qc_%s.pdf", name))
  png_path <- file.path(out_dir, sprintf("qc_%s.png", name))
  
  log_info(sprintf("Saving plot to %s and %s...", pdf_path, png_path), dataset = dataset_id, stage = "qc_plots")
  
  ggsave(pdf_path, plot = info$plot, width = info$w, height = info$h, device = "pdf")
  ggsave(png_path, plot = info$plot, width = info$w, height = info$h, dpi = 300, device = "png")
  
  # Calculate input checksum
  input_checksum <- calculate_file_checksum(input_rds)
  git_commit <- tryCatch({
    trimws(system("git rev-parse HEAD", intern = TRUE))
  }, error = function(e) "NO_GIT_COMMIT")
  
  # Collect figure index information
  fig_index_rows[[paste0(name, "_pdf")]] <- list(
    figure_path = pdf_path,
    dataset = dataset_id,
    processing_stage = "pre-filter_QC",
    analysis_method = name,
    parameters = sprintf("w=%d;h=%d", info$w, info$h),
    input_object_checksum = input_checksum,
    generating_script = "scripts/R/generate_qc_plots.R",
    snakemake_rule = snakemake_rule,
    git_commit = git_commit
  )
  fig_index_rows[[paste0(name, "_png")]] <- list(
    figure_path = png_path,
    dataset = dataset_id,
    processing_stage = "pre-filter_QC",
    analysis_method = name,
    parameters = sprintf("w=%d;h=%d", info$w, info$h),
    input_object_checksum = input_checksum,
    generating_script = "scripts/R/generate_qc_plots.R",
    snakemake_rule = snakemake_rule,
    git_commit = git_commit
  )
}

# Generate local figure index TSV
fig_index_df <- do.call(rbind, lapply(fig_index_rows, as.data.frame))
fig_index_tsv_path <- file.path(out_dir, "figure_index.tsv")
write.table(fig_index_df, fig_index_tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)
log_info(sprintf("Local figure index TSV saved to %s", fig_index_tsv_path), dataset = dataset_id, stage = "qc_plots")

# Generate Machine-Readable QC Table
log_info("Calculating descriptive QC statistics...", dataset = dataset_id, stage = "qc_plots")
qc_stats <- list()
for (col in required_cols) {
  vals <- meta[[col]]
  q <- quantile(vals, probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99))
  
  qc_stats[[col]] <- data.frame(
    dataset_id = dataset_id,
    metric = col,
    n_cells = length(vals),
    mean = mean(vals),
    sd = sd(vals),
    min = min(vals),
    max = max(vals),
    p1 = q[["1%"]],
    p5 = q[["5%"]],
    p10 = q[["10%"]],
    p25 = q[["25%"]],
    p50 = q[["50%"]],
    p75 = q[["75%"]],
    p90 = q[["90%"]],
    p95 = q[["95%"]],
    p99 = q[["99%"]],
    stringsAsFactors = FALSE
  )
}
qc_stats_df <- do.call(rbind, qc_stats)
qc_table_path <- file.path(out_dir, "qc_metrics_summary.tsv")
write.table(qc_stats_df, qc_table_path, sep = "\t", row.names = FALSE, quote = FALSE)
log_info(sprintf("Machine-readable QC stats saved to %s", qc_table_path), dataset = dataset_id, stage = "qc_plots")

# Save provenance record
inputs <- list(dataset_rds = input_rds, script = "scripts/R/generate_qc_plots.R")
outputs <- list(
  violins_pdf = file.path(out_dir, "qc_violins.pdf"),
  violins_png = file.path(out_dir, "qc_violins.png"),
  scatter_pdf = file.path(out_dir, "qc_scatter.pdf"),
  scatter_png = file.path(out_dir, "qc_scatter.png"),
  histograms_pdf = file.path(out_dir, "qc_histograms.pdf"),
  histograms_png = file.path(out_dir, "qc_histograms.png"),
  metrics_pdf = file.path(out_dir, "qc_metrics.pdf"),
  metrics_png = file.path(out_dir, "qc_metrics.png"),
  fig_index = fig_index_tsv_path,
  qc_table = qc_table_path
)
parameters <- list(dataset_id = dataset_id)
prov_record <- record_provenance("QCPlotGeneration", inputs, outputs, parameters, dataset = dataset_id)
save_provenance_json(prov_record, out_prov)
log_info(sprintf("Provenance record saved to %s", out_prov), dataset = dataset_id, stage = "qc_plots")

log_info("QC plot generation completed successfully.", dataset = dataset_id, stage = "qc_plots")
