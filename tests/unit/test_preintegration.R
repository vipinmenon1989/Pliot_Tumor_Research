# tests/unit/test_preintegration.R
library(yaml)
library(Seurat)

source("scripts/R/logging_utils.R")

setup_strict_logging()
log_info("Starting Pre-Integration Baseline validation tests...", stage = "test_preintegration")

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
  log_info("No datasets argument provided. Loading config/config.yaml...", stage = "test_preintegration")
  config <- read_yaml("config/config.yaml")
  datasets <- config$datasets
}

log_info(sprintf("Active datasets to validate: %s", paste(datasets, collapse = ", ")), stage = "test_preintegration")

is_test_mode <- any(grepl("sample", datasets))
config_path <- if (is_test_mode) "config/config.test.yaml" else "config/config.yaml"
config <- read_yaml(config_path)

# Load combined object
combined_rds_path <- "results/combined/pre_integration/combined_preintegration.rds"
log_info(sprintf("Loading combined Seurat object from %s...", combined_rds_path), stage = "test_preintegration")
if (!file.exists(combined_rds_path)) {
  log_error("Combined Seurat object does not exist!", stage = "test_preintegration")
  stop("Missing combined Seurat object.")
}

combined_obj <- readRDS(combined_rds_path)

# 1. Verify Cell Counts
expected_counts <- list()
if (is_test_mode) {
  for (ds in datasets) {
    expected_counts[[ds]] <- 150
  }
} else {
  expected_counts <- list(
    MPNST_1 = 7615,
    MPNST_2 = 2284,
    MPNST_3 = 2940,
    MPNST_4 = 6877
  )
}

expected_total <- sum(unlist(expected_counts))
actual_total <- ncol(combined_obj)

log_info(sprintf("Verifying total cell count (expected: %d, actual: %d)...", expected_total, actual_total), stage = "test_preintegration")
if (actual_total != expected_total) {
  log_error(sprintf("Total cell count mismatch: expected %d, got %d", expected_total, actual_total), stage = "test_preintegration")
  stop("Cell count mismatch!")
}

for (ds in datasets) {
  ds_cells <- sum(combined_obj$sample_id == ds)
  log_info(sprintf("Dataset %s: expected %d cells, found %d cells...", ds, expected_counts[[ds]], ds_cells), stage = "test_preintegration")
  if (ds_cells != expected_counts[[ds]]) {
    log_error(sprintf("Cell count mismatch for dataset %s: expected %d, got %d", ds, expected_counts[[ds]], ds_cells), stage = "test_preintegration")
    stop("Cell count mismatch per dataset!")
  }
}

# 2. Check for Duplicate Cell Names
log_info("Verifying that there are no duplicate cell names...", stage = "test_preintegration")
if (any(duplicated(colnames(combined_obj)))) {
  log_error("Duplicate cell names detected in combined object!", stage = "test_preintegration")
  stop("Duplicate cell names!")
}

# 3. Check for Namespaced M6/M7 Clustering Metadata
log_info("Verifying namespaced resolution sweep metadata columns...", stage = "test_preintegration")

for (ds in datasets) {
  # Dynamically detect resolutions for this dataset from its namespaced columns in metadata
  col_names <- colnames(combined_obj@meta.data)
  prefix <- sprintf("preint_%s_res_", ds)
  res_cols <- col_names[grepl(sprintf("^%s", prefix), col_names)]
  resolutions <- gsub(prefix, "", res_cols)
  
  for (res in resolutions) {
    target_col <- sprintf("preint_%s_res_%s", ds, res)
    if (!(target_col %in% colnames(combined_obj@meta.data))) {
      log_error(sprintf("Missing namespaced clustering column %s!", target_col), stage = "test_preintegration")
      stop("Missing namespaced clustering column.")
    }
    
    # Assert NA values for cells NOT in this dataset
    ds_cells_indices <- combined_obj$sample_id == ds
    other_cells_indices <- !ds_cells_indices
    
    non_na_other <- sum(!is.na(combined_obj@meta.data[other_cells_indices, target_col]))
    if (non_na_other > 0) {
      log_error(sprintf("Clustering column %s has non-NA values for cells outside dataset %s!", target_col, ds), stage = "test_preintegration")
      stop("Clustering columns not properly dataset-specific!")
    }
    
    na_ds <- sum(is.na(combined_obj@meta.data[ds_cells_indices, target_col]))
    if (na_ds > 0) {
      log_error(sprintf("Clustering column %s has NA values for cells within dataset %s!", target_col, ds), stage = "test_preintegration")
      stop("Clustering column has missing values for own cells!")
    }
  }
}

# 4. Check Recommended / Alternative Clustering Fields
recommended_fields <- c(
  "preint_M6_recommended_resolution",
  "preint_M7_recommended_resolution",
  "preint_M7_alternative_resolution",
  "preint_recommended_cluster",
  "preint_alternative_cluster"
)

for (f in recommended_fields) {
  if (!(f %in% colnames(combined_obj@meta.data))) {
    log_error(sprintf("Missing required recommended clustering field: %s", f), stage = "test_preintegration")
    stop("Missing recommended field.")
  }
}

# Validate prefix formatting of unique recommended clusters
rec_clusters <- combined_obj$preint_recommended_cluster
for (ds in datasets) {
  ds_rec_clusters <- unique(rec_clusters[combined_obj$sample_id == ds])
  pattern <- sprintf("^%s_C[0-9]{2}$", ds)
  non_matching <- ds_rec_clusters[!grepl(pattern, ds_rec_clusters)]
  if (length(non_matching) > 0) {
    log_error(sprintf("Recommended cluster label formatting error in dataset %s: %s", ds, paste(non_matching, collapse = ", ")), stage = "test_preintegration")
    stop("Recommended cluster formatting mismatch!")
  }
}

# 5. Verify Reductions and Absence of Integration Embeddings
log_info("Verifying dimensional reductions...", stage = "test_preintegration")
reds <- Reductions(combined_obj)
log_info(sprintf("Observed reductions in combined object: %s", paste(reds, collapse = ", ")), stage = "test_preintegration")

if (!("pca" %in% reds)) {
  log_error("Shared PCA reduction does not exist!", stage = "test_preintegration")
  stop("Missing shared PCA.")
}

if (!("umap_preintegration" %in% reds)) {
  log_error("Shared UMAP reduction 'umap_preintegration' does not exist!", stage = "test_preintegration")
  stop("Missing shared UMAP.")
}

# Ensure no batch correction or integrated reductions exist
integrated_names <- c("harmony", "cca", "rpca", "mnn", "scvi", "bbknn", "combat")
for (iname in integrated_names) {
  matching <- reds[grepl(iname, reds, ignore.case = TRUE)]
  if (length(matching) > 0) {
    log_error(sprintf("CRITICAL ERROR: Integrated reduction '%s' detected in pre-integration baseline!", matching), stage = "test_preintegration")
    stop("Accidental integration detected!")
  }
}

# 6. Verify Reports and TSV Files Exist
log_info("Verifying existence of required reports and summaries...", stage = "test_preintegration")
required_files <- c(
  "reports/PRE_INTEGRATION_ASSESSMENT.md",
  "reports/INTEGRATION_PREPARATION.md",
  "reports/milestones/M8_REPORT.md",
  "reports/combined/pre_integration/metadata_dictionary.tsv",
  "reports/combined/pre_integration/metadata_inventory.tsv",
  "reports/combined/pre_integration/neighborhood_mixing_summary.tsv",
  "reports/combined/pre_integration/composition_dataset.tsv",
  "reports/combined/pre_integration/composition_cluster.tsv",
  "reports/combined/pre_integration/confounding_summary.tsv",
  "reports/combined/pre_integration/integrity_summary.tsv"
)

for (rf in required_files) {
  if (!file.exists(rf) || file.size(rf) == 0) {
    log_error(sprintf("Required file %s is missing or empty!", rf), stage = "test_preintegration")
    stop("Missing or empty report/summary file.")
  }
}

# 7. Verify Figure Index Contains M8 figures
log_info("Verifying global figure index registrations...", stage = "test_preintegration")
fig_idx_path <- "reports/FIGURE_INDEX.tsv"
if (!file.exists(fig_idx_path)) {
  log_error("Global figure index reports/FIGURE_INDEX.tsv is missing!", stage = "test_preintegration")
  stop("Missing global figure index.")
}

fig_idx <- read.table(fig_idx_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
m8_entries <- fig_idx[fig_idx$dataset == "combined" & fig_idx$processing_stage == "pre_integration", ]

log_info(sprintf("Found %d figure index entries for pre_integration stage.", nrow(m8_entries)), stage = "test_preintegration")
if (nrow(m8_entries) < 10) {
  log_error("Under-registered pre-integration figures in FIGURE_INDEX.tsv! Expected at least 10 figures.", stage = "test_preintegration")
  stop("Incomplete figure index registration.")
}

log_info("All Pre-Integration Baseline validation tests PASSED.", stage = "test_preintegration")
