# Generate Deterministic Synthetic Seurat Object for Smoke Tests and CI
# Implements multiple constituent datasets, samples, patients, conditions, and cell types.

options(stringsAsFactors = FALSE)
set.seed(42)

library(Seurat)
library(Matrix)
source("scripts/R/logging_utils.R")
source("scripts/R/provenance_utils.R")

setup_strict_logging()
log_info("Starting synthetic data generation...", stage = "synthetic_generation")

# Define structures
n_cells <- 600
genes <- c(
  paste0("T-Marker-", 1:20),       # T-cell specific
  paste0("Myeloid-Marker-", 1:20),  # Myeloid specific
  paste0("Tumor-Marker-", 1:20),    # Tumor specific
  paste0("MT-", 1:10),              # Mitochondrial
  paste0("RPL-", 1:5),              # Ribosomal (large)
  paste0("RPS-", 1:5),              # Ribosomal (small)
  paste0("Gene-", 1:40)             # Housekeeping/background
)
n_genes <- length(genes)

# Cell metadata mapping
cell_types <- c("T_cell", "Myeloid", "Tumor_cell")
samples <- c("sample_1", "sample_2", "sample_3", "sample_4")
patients <- c("patient_A", "patient_B")
conditions <- c("Control", "Tumor")

# Deterministically assign cells to cell types, samples, and patients
cell_type_assign <- rep(cell_types, length.out = n_cells)
sample_assign <- rep(samples, length.out = n_cells)

# Metadata mappings
# patient_A: sample_1 (Control), sample_2 (Tumor)
# patient_B: sample_3 (Control), sample_4 (Tumor)
patient_map <- c(
  "sample_1" = "patient_A",
  "sample_2" = "patient_A",
  "sample_3" = "patient_B",
  "sample_4" = "patient_B"
)

condition_map <- c(
  "sample_1" = "Control",
  "sample_2" = "Tumor",
  "sample_3" = "Control",
  "sample_4" = "Tumor"
)

patient_assign <- patient_map[sample_assign]
condition_assign <- condition_map[sample_assign]

# Generate expression matrix (counts)
# We will use poisson distribution with rates dependent on cell type and gene type
counts_mat <- matrix(0, nrow = n_genes, ncol = n_cells)
rownames(counts_mat) <- genes
colnames(counts_mat) <- paste0("cell_", 1:n_cells)

for (j in 1:n_cells) {
  ct <- cell_type_assign[j]
  sm <- sample_assign[j]

  # Base depth variation (technical effect by sample)
  base_depth <- switch(sm,
    "sample_1" = 1500,
    "sample_2" = 2000,
    "sample_3" = 1200,
    "sample_4" = 1800
  )

  # Calculate rate factors for gene classes
  rates <- rep(0.01, n_genes) # default background rate
  names(rates) <- genes

  # Set markers high in respective cell types
  if (ct == "T_cell") {
    rates[grep("T-Marker-", genes)] <- 0.15
  } else if (ct == "Myeloid") {
    rates[grep("Myeloid-Marker-", genes)] <- 0.15
  } else if (ct == "Tumor_cell") {
    rates[grep("Tumor-Marker-", genes)] <- 0.20
  }

  # Mitochondrial rate (with slight sample-specific MT contamination effect)
  mt_rate <- 0.03
  if (sm == "sample_2") {
    mt_rate <- 0.08 # higher MT in sample 2 (technical artifact)
  }
  rates[grep("^MT-", genes)] <- mt_rate

  # Ribosomal rate
  rates[grep("^RP[L|S]-", genes)] <- 0.05

  # Draw counts
  cell_rates <- rates * base_depth
  counts_mat[, j] <- rpois(n_genes, lambda = cell_rates)
}

# Convert to sparse matrix
counts_sparse <- as(counts_mat, "CsparseMatrix")

# Create metadata dataframe
meta_df <- data.frame(
  sample_id = sample_assign,
  patient_id = patient_assign,
  condition = condition_assign,
  cell_type_annotation = cell_type_assign,
  orig.ident = sample_assign,
  row.names = colnames(counts_sparse),
  stringsAsFactors = FALSE
)

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = counts_sparse,
  project = "MPNST_Synthetic",
  meta.data = meta_df,
  min.cells = 0,
  min.features = 0
)

# Calculate MT and Ribosomal percentages
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[L|S]-")

# Add some legacy columns to mimic the real object as required by section 0.3
seurat_obj$unintegrated_clusters <- as.factor(sample(1:5, n_cells, replace = TRUE))
seurat_obj$seurat_clusters <- seurat_obj$unintegrated_clusters

# Save object
dir.create("data/synthetic", recursive = TRUE, showWarnings = FALSE)
output_path <- "data/synthetic/synthetic_mpnst.rds"
saveRDS(seurat_obj, file = output_path)

log_info(sprintf("Synthetic object saved to %s", output_path), stage = "synthetic_generation")
log_info(sprintf("Dimensions: %d genes, %d cells", nrow(seurat_obj), ncol(seurat_obj)), stage = "synthetic_generation")

# Save provenance information
inputs <- list(generator_script = "scripts/R/generate_synthetic_data.R")
outputs <- list(synthetic_rds = output_path)
parameters <- list(n_cells = n_cells, n_genes = n_genes, seed = 42)
prov_record <- record_provenance("GenerateSyntheticData", inputs, outputs, parameters)
save_provenance_json(prov_record, "reports/audits/synthetic_data_provenance.json")

log_info("Synthetic data generation finished.", stage = "synthetic_generation")
