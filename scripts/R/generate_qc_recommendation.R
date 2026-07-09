# scripts/R/generate_qc_recommendation.R
options(stringsAsFactors = FALSE)
set.seed(42)

library(Seurat)
library(yaml)
library(jsonlite)

source("scripts/R/logging_utils.R")
source("scripts/R/provenance_utils.R")

setup_strict_logging()
log_info("Starting QC recommendation generation...", stage = "qc_recommendation")

# Argument Parsing
args <- commandArgs(trailingOnly = TRUE)
input_rds <- ""
dataset_id <- ""
config_path <- "config/config.yaml"
out_md <- ""
out_prov <- ""

i <- 1
while (i <= length(args)) {
  if (args[i] == "--input") {
    input_rds <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--dataset-id") {
    dataset_id <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--config") {
    config_path <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-md") {
    out_md <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-prov") {
    out_prov <- args[i+1]
    i <- i + 2
  } else {
    stop(sprintf("Unknown argument: %s", args[i]))
  }
}

if (input_rds == "" || dataset_id == "" || out_md == "" || out_prov == "") {
  stop("Missing required arguments.")
}

log_info(sprintf("Loading Seurat object from %s...", input_rds), dataset = dataset_id, stage = "qc_recommendation")
subset_obj <- readRDS(input_rds)

log_info(sprintf("Loading config file from %s...", config_path), dataset = dataset_id, stage = "qc_recommendation")
config <- read_yaml(config_path)

# Extract default QC thresholds from config
qc_cfg <- config$qc
min_features <- qc_cfg$min_features
max_features <- qc_cfg$max_features
min_counts <- qc_cfg$min_counts
max_counts <- qc_cfg$max_counts
max_percent_mt <- qc_cfg$max_percent_mt
max_percent_ribo <- qc_cfg$max_percent_ribo

meta <- subset_obj@meta.data
total_cells <- nrow(meta)

# Compute metrics
# Percentiles
get_percentiles <- function(v) {
  quantile(v, probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99))
}
q_counts <- get_percentiles(meta$nCount_RNA)
q_features <- get_percentiles(meta$nFeature_RNA)
q_mt <- get_percentiles(meta$percent.mt)
q_ribo <- get_percentiles(meta$percent.ribo)

# Outliers and expected filtering
low_feat_cells <- sum(meta$nFeature_RNA < min_features)
high_feat_cells <- sum(meta$nFeature_RNA > max_features)
low_count_cells <- sum(meta$nCount_RNA < min_counts)
high_count_cells <- sum(meta$nCount_RNA > max_counts)
high_mt_cells <- sum(meta$percent.mt > max_percent_mt)
high_ribo_cells <- sum(meta$percent.ribo > max_percent_ribo)

# Total filtered cells (combining all filters)
filtered_out <- (meta$nFeature_RNA < min_features) | 
                 (meta$nFeature_RNA > max_features) | 
                 (meta$nCount_RNA < min_counts) | 
                 (meta$nCount_RNA > max_counts) | 
                 (meta$percent.mt > max_percent_mt) | 
                 (meta$percent.ribo > max_percent_ribo)
total_removed <- sum(filtered_out)
pct_removed <- (total_removed / total_cells) * 100

# Write the QC recommendation markdown
dir.create(dirname(out_md), recursive = TRUE, showWarnings = FALSE)

# Generate Markdown Content
md_content <- c(
  sprintf("# QC Recommendation Report - %s", dataset_id),
  sprintf("*Generated on: %s*", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. Observed Distributions Summary",
  "",
  "The table below summarizes the observed distributions and key percentiles for the cell quality control metrics prior to filtering.",
  "",
  "| Metric | Min | Mean | Median | Max | 1% | 5% | 10% | 90% | 95% | 99% |",
  "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
  sprintf("| **nCount_RNA** | %d | %.1f | %d | %d | %d | %d | %d | %d | %d | %d |", 
          as.integer(min(meta$nCount_RNA)), mean(meta$nCount_RNA), as.integer(median(meta$nCount_RNA)), as.integer(max(meta$nCount_RNA)),
          as.integer(q_counts["1%"]), as.integer(q_counts["5%"]), as.integer(q_counts["10%"]), as.integer(q_counts["90%"]), as.integer(q_counts["95%"]), as.integer(q_counts["99%"])),
  sprintf("| **nFeature_RNA** | %d | %.1f | %d | %d | %d | %d | %d | %d | %d | %d |", 
          as.integer(min(meta$nFeature_RNA)), mean(meta$nFeature_RNA), as.integer(median(meta$nFeature_RNA)), as.integer(max(meta$nFeature_RNA)),
          as.integer(q_features["1%"]), as.integer(q_features["5%"]), as.integer(q_features["10%"]), as.integer(q_features["90%"]), as.integer(q_features["95%"]), as.integer(q_features["99%"])),
  sprintf("| **percent.mt** | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% |", 
          min(meta$percent.mt), mean(meta$percent.mt), median(meta$percent.mt), max(meta$percent.mt),
          q_mt["1%"], q_mt["5%"], q_mt["10%"], q_mt["90%"], q_mt["95%"], q_mt["99%"]),
  sprintf("| **percent.ribo** | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% |", 
          min(meta$percent.ribo), mean(meta$percent.ribo), median(meta$percent.ribo), max(meta$percent.ribo),
          q_ribo["1%"], q_ribo["5%"], q_ribo["10%"], q_ribo["90%"], q_ribo["95%"], q_ribo["99%"]),
  "",
  "## 2. Potential Outliers",
  "",
  "Based on statistical distributions (e.g. median +/- 3 Median Absolute Deviations (MAD) or standard deviation thresholds):",
  "",
  sprintf("- **Library Depth (nCount_RNA)**: The median value is %d. Outliers are typically defined as cells with counts below %d or above %d.", 
          as.integer(median(meta$nCount_RNA)), 
          as.integer(median(meta$nCount_RNA) - 3 * mad(meta$nCount_RNA)), 
          as.integer(median(meta$nCount_RNA) + 3 * mad(meta$nCount_RNA))),
  sprintf("- **Genes Detected (nFeature_RNA)**: The median value is %d. Outliers are typically cells with fewer than %d genes (likely empty droplets or low-quality cells) or more than %d genes (potential doublets).", 
          as.integer(median(meta$nFeature_RNA)), 
          as.integer(median(meta$nFeature_RNA) - 3 * mad(meta$nFeature_RNA)), 
          as.integer(median(meta$nFeature_RNA) + 3 * mad(meta$nFeature_RNA))),
  sprintf("- **Mitochondrial Fraction (percent.mt)**: The median is %.2f%%, and the 95th percentile is %.2f%%. Cells exceeding %.2f%% show elevated mitochondrial expression indicating cellular stress or lysis.", 
          median(meta$percent.mt), q_mt["95%"], q_mt["95%"]),
  "",
  "## 3. Recommended Thresholds and Rationale",
  "",
  "Below are the recommended thresholds driven by the project configuration:",
  "",
  "| QC Metric | Recommended Threshold | Rationale |",
  "| --- | --- | --- |",
  sprintf("| **Min Features (nFeature_RNA)** | `> %d` | Exclude low-complexity droplets/dead cells that do not contain sufficient biological signal. |", min_features),
  sprintf("| **Max Features (nFeature_RNA)** | `< %d` | Exclude potential doublets or multi-cell aggregates. |", max_features),
  sprintf("| **Min Counts (nCount_RNA)** | `> %d` | Ensure sufficient library depth for robust gene expression estimation. |", min_counts),
  sprintf("| **Max Counts (nCount_RNA)** | `< %d` | Exclude cells with abnormally high UMI counts, indicating technical artifacts or doublets. |", max_counts),
  sprintf("| **Max percent.mt** | `< %.1f%%` | Standard filter to remove dying or damaged cells which release cytoplasmic RNA and retain mitochondrial transcripts. |", max_percent_mt),
  sprintf("| **Max percent.ribo** | `< %.1f%%` | Eliminate cells with extremely high ribosomal expression, which may represent technical bias or specific translation stress. |", max_percent_ribo),
  "",
  "## 4. Expected Filtering Impact",
  "",
  "Here is the estimated impact of applying each of the recommended thresholds independently and in combination:",
  "",
  "| Filter Metric | Threshold | Cells Excluded | % Excluded |",
  "| --- | --- | --- | --- |",
  sprintf("| **Min Features** | `< %d` | %d | %.2f%% |", min_features, low_feat_cells, (low_feat_cells / total_cells) * 100),
  sprintf("| **Max Features** | `> %d` | %d | %.2f%% |", max_features, high_feat_cells, (high_feat_cells / total_cells) * 100),
  sprintf("| **Min Counts** | `< %d` | %d | %.2f%% |", min_counts, low_count_cells, (low_count_cells / total_cells) * 100),
  sprintf("| **Max Counts** | `> %d` | %d | %.2f%% |", max_counts, high_count_cells, (high_count_cells / total_cells) * 100),
  sprintf("| **Max percent.mt** | `> %.1f%%` | %d | %.2f%% |", max_percent_mt, high_mt_cells, (high_mt_cells / total_cells) * 100),
  sprintf("| **Max percent.ribo** | `> %.1f%%` | %d | %.2f%% |", max_percent_ribo, high_ribo_cells, (high_ribo_cells / total_cells) * 100),
  sprintf("| **Combined Filters** | **All Above** | **%d** | **%.2f%%** |", total_removed, pct_removed),
  "",
  "### Summary of Expected Kept Cells:",
  sprintf("- **Total cells before filtering**: %d", total_cells),
  sprintf("- **Expected cells removed**: %d (%.2f%%)", total_removed, pct_removed),
  sprintf("- **Expected cells retained**: %d (%.2f%%)", total_cells - total_removed, 100 - pct_removed),
  "",
  "## 5. Potential Biological & Technical Risks",
  "",
  "1. **Mitochondrial Bias**: MPNST is a soft tissue sarcoma and tumor cells may exhibit metabolic shifts leading to higher baseline mitochondrial transcript percentages. Applying a strict 15% threshold could selectively remove tumor sub-clones or hypoxic cells. However, 15% provides a standard compromise for human cells.",
  "2. **Ribosomal Variation**: Actively dividing tumor cells or specific cell classes (such as plasma cells) may have high ribosomal fractions. A 20% ribosomal threshold should be evaluated for potential cell-class depletion.",
  "3. **Doublet Detection Confounding**: Cells with high counts are marked for exclusion by the max counts threshold, but a formal doublet detection (Milestone 3) will be needed to distinguish true high-transcript cells from doublets.",
  "",
  "## 6. Questions Requiring Researcher Approval",
  "",
  "Please review and approve the following settings prior to Milestone 3 (M3) filtering:",
  "",
  sprintf("1. **Do you approve the uniform minimum gene limit of %d features across this dataset?**", min_features),
  sprintf("2. **Do you approve the mitochondrial limit of %.1f%% for %s?** (Excludes %d cells)", max_percent_mt, dataset_id, high_mt_cells),
  sprintf("3. **Do you approve the ribosomal limit of %.1f%% for %s?** (Excludes %d cells)", max_percent_ribo, dataset_id, high_ribo_cells),
  sprintf("4. **Do you approve the overall threshold combination which will exclude %.2f%% of cells?**", pct_removed),
  "",
  "**NOTE: NO FILTERS HAVE BEEN APPLIED AT THIS STAGE. The dataset remains completely intact.**"
)

writeLines(md_content, out_md)
log_info(sprintf("QC Recommendation report saved to %s", out_md), dataset = dataset_id, stage = "qc_recommendation")

# Save provenance record
inputs <- list(dataset_rds = input_rds, script = "scripts/R/generate_qc_recommendation.R")
outputs <- list(recommendation_md = out_md)
parameters <- list(
  dataset_id = dataset_id,
  min_features = min_features,
  max_features = max_features,
  min_counts = min_counts,
  max_counts = max_counts,
  max_percent_mt = max_percent_mt,
  max_percent_ribo = max_percent_ribo
)
prov_record <- record_provenance("QCRecommendationGeneration", inputs, outputs, parameters, dataset = dataset_id)
save_provenance_json(prov_record, out_prov)
log_info(sprintf("Provenance record saved to %s", out_prov), dataset = dataset_id, stage = "qc_recommendation")

log_info("QC recommendation generation completed successfully.", dataset = dataset_id, stage = "qc_recommendation")
