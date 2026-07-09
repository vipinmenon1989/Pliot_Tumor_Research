# scripts/R/extract_dataset.R
options(stringsAsFactors = FALSE)
set.seed(42)

library(Seurat)
library(jsonlite)

source("scripts/R/logging_utils.R")
source("scripts/R/provenance_utils.R")

setup_strict_logging()
log_info("Starting dataset extraction...", stage = "extraction")

# Argument Parsing
args <- commandArgs(trailingOnly = TRUE)
input_rds <- ""
dataset_id <- ""
id_column <- "sample_id"
out_rds <- ""
out_val_md <- ""
out_prov <- ""
out_manifest <- ""

i <- 1
while (i <= length(args)) {
  if (args[i] == "--input") {
    input_rds <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--dataset-id") {
    dataset_id <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--id-column") {
    id_column <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-rds") {
    out_rds <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-val-md") {
    out_val_md <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-prov") {
    out_prov <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--out-manifest") {
    out_manifest <- args[i+1]
    i <- i + 2
  } else {
    stop(sprintf("Unknown argument: %s", args[i]))
  }
}

if (input_rds == "" || dataset_id == "" || out_rds == "" || out_val_md == "" || out_prov == "" || out_manifest == "") {
  stop("Missing required arguments.")
}

log_info(sprintf("Loading original Seurat object from %s...", input_rds), dataset = dataset_id, stage = "extraction")
seurat_obj <- readRDS(input_rds)

log_info(sprintf("Verifying cells for dataset ID %s in column %s...", dataset_id, id_column), dataset = dataset_id, stage = "extraction")
if (!id_column %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Metadata column %s not found in input object.", id_column))
}

cells_in_dataset <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[id_column]] == dataset_id]
n_cells <- length(cells_in_dataset)
log_info(sprintf("Found %d cells belonging to dataset %s.", n_cells, dataset_id), dataset = dataset_id, stage = "extraction")

if (n_cells == 0) {
  stop(sprintf("No cells found for dataset ID %s.", dataset_id))
}

# Create subset object
log_info("Subsetting Seurat object...", dataset = dataset_id, stage = "extraction")
subset_obj <- subset(seurat_obj, cells = cells_in_dataset)

# Calculate percent.ribo and add to metadata if not present
log_info("Calculating percent.ribo...", dataset = dataset_id, stage = "extraction")
# Use ^RP[SL] pattern to match ribosomal genes
subset_obj[["percent.ribo"]] <- PercentageFeatureSet(subset_obj, pattern = "^RP[SL]")

# Preserve: assays, counts, metadata, layers, reductions (inventory only), provenance
# Save reductions inventory before dropping
reductions_present <- names(seurat_obj@reductions)
log_info(sprintf("Original reductions inventoried: %s", paste(reductions_present, collapse = ", ")), dataset = dataset_id, stage = "extraction")

# Drop reductions in the output object to prevent contamination
subset_obj@reductions <- list()
log_info("Removed dimensional reductions from working object.", dataset = dataset_id, stage = "extraction")

# Save the subsetted Seurat object
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
log_info(sprintf("Saving subsetted object to %s...", out_rds), dataset = dataset_id, stage = "extraction")
saveRDS(subset_obj, file = out_rds)
log_info("Working object saved successfully.", dataset = dataset_id, stage = "extraction")

# Calculate checksum
rds_checksum <- calculate_file_checksum(out_rds)
log_info(sprintf("MD5 Checksum: %s", rds_checksum), dataset = dataset_id, stage = "extraction")

# Collect dataset validation metrics
val_metrics <- list(
  dataset_id = dataset_id,
  cells = ncol(subset_obj),
  genes = nrow(subset_obj),
  assays = names(subset_obj@assays),
  default_assay = DefaultAssay(subset_obj),
  metadata_columns = colnames(subset_obj@meta.data),
  object_size_bytes = as.numeric(object.size(subset_obj)),
  checksum = rds_checksum,
  file_path = out_rds
)

# Assay details
assay_details <- list()
for (a in names(subset_obj@assays)) {
  assay_obj <- subset_obj[[a]]
  layers <- Layers(subset_obj, assay = a)
  
  # Check if counts exist
  has_counts <- FALSE
  if (inherits(assay_obj, "Assay5")) {
    has_counts <- any(grepl("counts", layers))
  } else {
    counts_data <- tryCatch({
      GetAssayData(subset_obj, assay = a, layer = "counts")
    }, error = function(e) {
      tryCatch({
        GetAssayData(subset_obj, assay = a, slot = "counts")
      }, error = function(e2) {
        NULL
      })
    })
    has_counts <- !is.null(counts_data) && (length(counts_data) > 0)
  }
  
  # Sum features
  features_count <- nrow(assay_obj)
  
  assay_details[[a]] <- list(
    class = as.character(class(assay_obj)),
    features = features_count,
    layers = layers,
    counts_present = has_counts
  )
}
val_metrics$assay_details <- assay_details

# Save machine-readable manifest (JSON)
dir.create(dirname(out_manifest), recursive = TRUE, showWarnings = FALSE)
writeLines(toJSON(val_metrics, auto_unbox = TRUE, pretty = TRUE), out_manifest)
log_info(sprintf("Machine-readable manifest saved to %s", out_manifest), dataset = dataset_id, stage = "extraction")

# Get git commit
git_commit <- tryCatch({
  trimws(system("git rev-parse HEAD", intern = TRUE))
}, error = function(e) "NO_GIT_COMMIT")

# Generate DATASET_VALIDATION.md report
dir.create(dirname(out_val_md), recursive = TRUE, showWarnings = FALSE)
val_md_content <- c(
  sprintf("# Dataset Validation Report - %s", dataset_id),
  sprintf("*Generated on: %s*", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. Object Overview",
  "",
  "| Metric | Value |",
  "| --- | --- |",
  sprintf("| **Dataset ID** | %s |", dataset_id),
  sprintf("| **Cells** | %s |", format(val_metrics$cells, big.mark = ",")),
  sprintf("| **Genes** | %s |", format(val_metrics$genes, big.mark = ",")),
  sprintf("| **Dimensions** | %d x %d |", val_metrics$genes, val_metrics$cells),
  sprintf("| **Default Assay** | %s |", val_metrics$default_assay),
  sprintf("| **Object Size in Memory** | %.2f MB |", val_metrics$object_size_bytes / (1024^2)),
  sprintf("| **File Path** | `%s` |", val_metrics$file_path),
  sprintf("| **MD5 Checksum** | `%s` |", val_metrics$checksum),
  "",
  "## 2. Assay Structure",
  "",
  "| Assay | Class | Features | Layers | Counts Present? |",
  "| --- | --- | --- | --- | --- |"
)

for (a in names(assay_details)) {
  ad <- assay_details[[a]]
  val_md_content <- c(
    val_md_content,
    sprintf("| %s | %s | %s | %s | %s |", 
            a, ad$class, format(ad$features, big.mark = ","), 
            paste(ad$layers, collapse = ", "), 
            if (ad$counts_present) "YES" else "NO")
  )
}

val_md_content <- c(
  val_md_content,
  "",
  "## 3. Metadata Columns",
  "",
  "The following metadata columns are preserved in the Seurat object:",
  "",
  paste0("- `", val_metrics$metadata_columns, "`", collapse = "\n"),
  "",
  "## 4. Dimensional Reductions (Inventory Only)",
  "",
  "The following legacy reductions were present in the original object and have been inventoried but removed from this working object:",
  "",
  if (length(reductions_present) > 0) paste0("- `", reductions_present, "`", collapse = "\n") else "No reductions found.",
  "",
  "## 5. Provenance",
  "",
  sprintf("- **Input File**: `%s`", input_rds),
  sprintf("- **Input Checksum**: `%s`", calculate_file_checksum(input_rds)),
  sprintf("- **Extraction Code**: `scripts/R/extract_dataset.R`"),
  sprintf("- **Git Commit**: `%s`", git_commit)
)

# Write MD
writeLines(val_md_content, out_val_md)
log_info(sprintf("Validation report saved to %s", out_val_md), dataset = dataset_id, stage = "extraction")

# Record provenance
inputs <- list(original_rds = input_rds, script = "scripts/R/extract_dataset.R")
outputs <- list(dataset_rds = out_rds, manifest = out_manifest, report = out_val_md)
parameters <- list(dataset_id = dataset_id, id_column = id_column)
prov_record <- record_provenance("DatasetExtraction", inputs, outputs, parameters, dataset = dataset_id)
save_provenance_json(prov_record, out_prov)
log_info(sprintf("Provenance record saved to %s", out_prov), dataset = dataset_id, stage = "extraction")

log_info("Dataset extraction completed successfully.", dataset = dataset_id, stage = "extraction")
