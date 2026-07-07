# Lightweight R script to collect additional metrics from the Seurat object
options(stringsAsFactors = FALSE)
library(Seurat)
library(jsonlite)

input_rds <- "/local/projects-t3/lilab/vmenon/Pilot_tumor/processed_mpnst.rds"
if (!file.exists(input_rds)) {
  input_rds <- "/autofs/projects-t3/lilab/vmenon/Pilot_tumor/processed_mpnst.rds"
}

message("Loading Seurat object...")
start_time <- Sys.time()
seurat_obj <- readRDS(input_rds)
end_time <- Sys.time()
load_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
message("Loaded in ", load_time, " seconds.")

# Memory usage
mem_bytes <- object.size(seurat_obj)
mem_gb <- mem_bytes / (1024^3)
mem_mb <- mem_bytes / (1024^2)
mem_usage_str <- if (mem_gb >= 1) sprintf("%.2f GB", mem_gb) else sprintf("%.2f MB", mem_mb)

# Duplicated features and cell names
all_genes_rna <- rownames(seurat_obj[["RNA"]])
dup_genes_rna <- sum(duplicated(all_genes_rna))
dup_cells <- sum(duplicated(colnames(seurat_obj)))

# Mitochondrial and ribosomal genes
# Check for MT genes (case-insensitive)
mt_genes <- grep("^MT-", all_genes_rna, ignore.case = TRUE, value = TRUE)
mt_count <- length(mt_genes)

# Check for ribosomal genes (case-insensitive, typically RPS/RPL)
ribo_genes <- grep("^RP[SL]", all_genes_rna, ignore.case = TRUE, value = TRUE)
ribo_count <- length(ribo_genes)

# Assay details
assays_list <- Assays(seurat_obj)
assay_details <- list()
for (a in assays_list) {
  assay_obj <- seurat_obj[[a]]
  
  # Dimensions
  dims <- dim(assay_obj)
  
  # Layers/slots
  layers_found <- c()
  v5_layers <- tryCatch({
    SeuratObject::Layers(assay_obj)
  }, error = function(e) NULL)
  
  if (!is.null(v5_layers) && length(v5_layers) > 0) {
    layers_found <- v5_layers
  } else {
    if (.hasSlot(assay_obj, "counts") && !is.null(assay_obj@counts) && length(assay_obj@counts) > 0) layers_found <- c(layers_found, "counts")
    if (.hasSlot(assay_obj, "data") && !is.null(assay_obj@data) && length(assay_obj@data) > 0) layers_found <- c(layers_found, "data")
    if (.hasSlot(assay_obj, "scale.data") && !is.null(assay_obj@scale.data) && length(assay_obj@scale.data) > 0) layers_found <- c(layers_found, "scale.data")
  }
  
  # Check counts, data, scale.data
  counts_present <- "counts" %in% layers_found || any(grepl("^counts", layers_found))
  data_present <- "data" %in% layers_found || any(grepl("^data", layers_found))
  scale_present <- "scale.data" %in% layers_found || any(grepl("^scale\\.data", layers_found))
  
  # Variable features
  var_features <- tryCatch({
    VariableFeatures(seurat_obj, assay = a)
  }, error = function(e) NULL)
  var_features_count <- length(var_features)
  
  assay_details[[a]] <- list(
    name = a,
    class = class(assay_obj)[1],
    features = dims[1],
    cells = dims[2],
    layers = layers_found,
    counts_present = counts_present,
    data_present = data_present,
    scale_present = scale_present,
    variable_features_count = var_features_count
  )
}

# Collect output
result <- list(
  estimated_memory_bytes = as.numeric(mem_bytes),
  estimated_memory_usage = mem_usage_str,
  load_time_sec = load_time,
  duplicated_genes_rna = dup_genes_rna,
  duplicated_cells = dup_cells,
  mitochondrial_genes_count = mt_count,
  mitochondrial_genes_sample = head(mt_genes, 10),
  ribosomal_genes_count = ribo_count,
  ribosomal_genes_sample = head(ribo_genes, 10),
  assays = assay_details
)

writeLines(toJSON(result, auto_unbox = TRUE, pretty = TRUE), "reports/additional_metrics.json")
message("Done saving additional metrics.")
