# Scientific and Computational Inventory of Seurat Object
*Generated on: 2026-07-07 15:30:00 EDT*

## 1. Object Summary

| Metric | Value |
| --- | --- |
| **Seurat Version** | 5.1.0 |
| **Object Class** | Seurat |
| **Object Dimensions (Default SCT)** | 29708 features x 22661 cells |
| **Object Size on Disk** | 7.23 GB |
| **Estimated Memory Usage** | 10.08 GB (10,827,181,496 bytes) |
| **Elapsed Load Time** | 45.277 seconds |
| **Default Assay** | SCT |
| **Active Identities Class** | factor (12 unique identity levels) |
| **Available Assays** | RNA, mnn.reconstructed, SCT |

### Available Assays Architecture
| Assay | Class | Features | Cells | Layers |
| --- | --- | --- | --- | --- |
| RNA | Assay5 | 31,764 | 22,661 | counts.MPNST_1, counts.MPNST_2, counts.MPNST_3, counts.MPNST_4, data.MPNST_1, data.MPNST_2, data.MPNST_3, data.MPNST_4, scale.data |
| mnn.reconstructed | Assay | 2,000 | 22,661 | data |
| SCT | SCTAssay | 29,708 | 22,661 | counts, data, scale.data |

### Cell Networks & Graphs
- RNA_nn
- RNA_snn
- SCT_nn
- SCT_snn

### Dimensional Reductions
- pca
- pca200
- pca100
- pca50
- umap.unintegrated
- pca100_cca
- umap.pca100_cca
- pca100_harmony
- umap.pca100_harmony
- pca100_mnn
- umap.pca100_mnn
- pca100_rpca
- umap.pca100_rpca
- pca100.sct
- pca100.sct_harmony
- umap.pca100.sct_harmony
- pca100.sct_cca
- umap.pca100.sct_cca

### Commands Stored Inside the Object
- NormalizeData.RNA
- JackStraw.RNA.pca
- ScoreJackStraw
- FindNeighbors.RNA.pca100
- RunUMAP.RNA.pca100
- FindVariableFeatures.RNA
- ScaleData.RNA
- RunPCA.RNA
- RunUMAP.RNA.pca100_harmony
- FindNeighbors.RNA.pca100_harmony
- FindNeighbors.RNA.pca100_cca
- RunUMAP.RNA.pca100_cca
- RunUMAP.RNA.pca100_mnn
- FindNeighbors.RNA.pca100_mnn
- FindNeighbors.RNA.pca100_rpca
- RunUMAP.RNA.pca100_rpca
- SCTransform.RNA
- RunPCA.SCT
- RunUMAP.SCT.pca100.sct_harmony
- FindNeighbors.SCT.pca100.sct_harmony
- FindNeighbors.SCT.pca100.sct_cca
- FindClusters
- RunUMAP.SCT.pca100.sct_cca

## 2. Assay Inventory

| Assay Name | Class | Dimensions | Available Layers | Counts Present? | Normalized Present? | Scaled Present? | Variable Features? | Default Assay Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RNA | Assay5 | 31764 x 22661 | counts.MPNST_1, counts.MPNST_2, counts.MPNST_3, counts.MPNST_4, data.MPNST_1, data.MPNST_2, data.MPNST_3, data.MPNST_4, scale.data | YES | YES | YES | YES (2309 features) | NO |
| mnn.reconstructed | Assay | 2000 x 22661 | data | NO | YES | NO | YES (2000 features) | NO |
| SCT | SCTAssay | 29708 x 22661 | counts, data, scale.data | YES | YES | YES | YES (3000 features) | YES |

### Reprocessing Feasibility Assessment
> [!IMPORTANT]
> **Independent Reprocessing Status**: **FEASIBLE**
> 
> The Seurat object contains split raw UMI count layers (`counts.MPNST_1`, `counts.MPNST_2`, `counts.MPNST_3`, and `counts.MPNST_4`) in the RNA assay. This structure preserves the raw, un-normalized counts for each individual library, permitting complete, independent reprocessing from the baseline transcript counts.

## 3. Metadata Inventory

| Column Name | Data Type | Unique Values | Missing Values | Candidate Biological Meaning | Recommended Usage |
| --- | --- | --- | --- | --- | --- |
| `orig.ident` | character | 4 | 0 | Initial sample identifier from library setup | Auxiliary ID. Do not use for splitting; use sample_id instead for consistency. |
| `nCount_RNA` | numeric | 12,933 | 0 | Total raw RNA UMI count per cell (library depth) | Primary QC filtering metric. Filter cells with abnormally low or high counts. |
| `nFeature_RNA` | integer | 6,289 | 0 | Number of unique genes detected in RNA assay per cell | Primary QC filtering metric. Filter cells with low gene diversity or doublets. |
| `sample_id` | character | 4 | 0 | Consistent sample metadata identifier | Primary dataset splitter. Use this column to split the object into individual datasets. |
| `percent.mt` | numeric | 14,130 | 0 | Mitochondrial transcript percentage (quality control metric) | Primary QC filtering metric. Filter cells with high mitochondrial expression. |
| `orig.anno` | character | 16 | 0 | Legacy cell-type annotations from original analysis | Reference only. Exclude from Phase 1 processing to prevent bias, keep for validation. |
| `unintegrated_clusters` | factor | 16 | 0 | Legacy cluster assignments | Legacy artifact. Exclude from all Phase 1 processing. |
| `seurat_clusters` | factor | 12 | 0 | Legacy cluster assignments | Legacy artifact. Exclude from all Phase 1 processing. |
| `pca100_harmony_clusters` | factor | 16 | 0 | Legacy cluster assignments | Legacy artifact. Exclude from all Phase 1 processing. |
| `pca100_cca_clusters` | factor | 10 | 0 | Legacy cluster assignments | Legacy artifact. Exclude from all Phase 1 processing. |
| `pca100_mnn_clusters` | factor | 13 | 0 | Legacy cluster assignments | Legacy artifact. Exclude from all Phase 1 processing. |
| `pca100_rpca_clusters` | factor | 17 | 0 | Legacy cluster assignments | Legacy artifact. Exclude from all Phase 1 processing. |
| `nCount_SCT` | numeric | 5,887 | 0 | Total SCTransformed UMI count per cell | Exclude. Will be re-computed during the new normalization stage. |
| `nFeature_SCT` | integer | 4,173 | 0 | Number of unique genes detected in SCT assay per cell | Exclude. Will be re-computed during the new normalization stage. |
| `pca100.sct_harmony_clusters` | factor | 13 | 0 | Legacy cluster assignments | Legacy artifact. Exclude from all Phase 1 processing. |
| `pca100.sct_cca_clusters` | factor | 12 | 0 | Legacy cluster assignments | Legacy artifact. Exclude from all Phase 1 processing. |

## 4. sample_id versus orig.ident

- **Unique sample_id values**: 4
- **Unique orig.ident values**: 4

### Frequency Tables

| Value | cells in `sample_id` | cells in `orig.ident` |
| --- | --- | --- |
| MPNST_1 | 8,338 | 8,338 |
| MPNST_2 | 2,830 | 2,830 |
| MPNST_3 | 3,682 | 3,682 |
| MPNST_4 | 7,811 | 7,811 |

### Complete Cross-Tabulation

| sample_id | orig.ident | Cell Count | Mismatch? |
| --- | --- | --- | --- |
| MPNST_1 | MPNST_1 | 8,338 | No |
| MPNST_2 | MPNST_2 | 2,830 | No |
| MPNST_3 | MPNST_3 | 3,682 | No |
| MPNST_4 | MPNST_4 | 7,811 | No |

### Identity and Recommendation Statement
> [!NOTE]
> **Identity Check**: `sample_id` and `orig.ident` are 100% identical in mapping and cell count. There are 0 mismatches between them.
> 
> **Recommendation**: It is recommended to use `sample_id` as the primary dataset splitter. Standardizing on `sample_id` provides a cleaner, biologically descriptive column name, avoiding potential confusion with the Seurat default initialization column `orig.ident` while yielding mathematically identical subsetting.

## 5. Candidate Experimental Variables

### Observed Facts
- **Patient Identifier**: None. No columns in the metadata map to patient demographics or patient IDs.
- **Dataset Identifier**: `sample_id` (4 distinct samples: MPNST_1, MPNST_2, MPNST_3, MPNST_4).
- **Sample Identifier**: `sample_id` (4 distinct samples).
- **Biological Condition**: None. No columns represent biological condition or clinical parameters (such as tumor grade or controls) in the Seurat metadata.
- **Batch Variable**: `sample_id` (or `orig.ident`), which likely confound sample identity and library preparation batch.
- **Technical Replicate Variable**: None. No replicate designations are present.

### Inference
- It is inferred that `sample_id` groups the 22,661 cells into 4 independent experimental datasets, representing 4 single-cell sequencing libraries of Malignant Peripheral Nerve Sheath Tumors (MPNST).
- We *cannot* infer patient origin (e.g. whether these samples come from 4 separate patients or represent multiple samples from the same patient) or technical batch history (e.g. whether they were run on the same flow cell or different days) based on the object's metadata alone.

## 6. Cell Statistics

### Cell Frequencies across candidate groupings:

| Category Value | cells per sample | cells per orig.ident | cells per patient | cells per dataset |
| --- | --- | --- | --- | --- |
| MPNST_1 | 8,338 | 8,338 | - | 8,338 |
| MPNST_2 | 2,830 | 2,830 | - | 2,830 |
| MPNST_3 | 3,682 | 3,682 | - | 3,682 |
| MPNST_4 | 7,811 | 7,811 | - | 7,811 |
| Unknown / Not Annotated | - | - | 22,661 | - |
| **TOTAL** | **22,661** | **22,661** | **22,661** | **22,661** |

## 7. Feature Statistics

- **Total Genes (RNA)**: 31,764
- **Total Genes (SCT)**: 29,708
- **Duplicated Genes**: 0
- **Duplicated Cell Names**: 0
- **Mitochondrial Genes (matching `^MT-`)**: 13 detected
  - Examples: MT-ND1, MT-ND2, MT-CO1, MT-CO2, MT-ATP8, MT-ATP6, MT-CO3, MT-ND3, MT-ND4L, MT-ND4
- **Ribosomal Genes (matching `^RP[SL]`)**: 433 detected
  - Examples: RPL22, RPL11, RPS6KA1, RPS8, RPS15AP11, RPL21P23, RPS7P4, RPL5P6, RPL5, RPL7P9
- **QC-related Metadata Present**: `nCount_RNA` (library size), `nFeature_RNA` (gene count), `percent.mt` (mitochondrial fraction)

## 8. Legacy Analysis Inventory

The following prior analysis structures exist in the object:

### 8.1 Legacy Clustering Columns
| Metadata Column | Cardinality (Classes) | Integration Method |
| --- | --- | --- |
| `seurat_clusters` | 12 | Original / Default |
| `unintegrated_clusters` | 16 | Unintegrated SCT |
| `pca100_harmony_clusters` | 16 | Harmony (RNA assay) |
| `pca100_cca_clusters` | 10 | CCA (RNA assay) |
| `pca100_mnn_clusters` | 13 | MNN (RNA assay) |
| `pca100_rpca_clusters` | 17 | RPCA (RNA assay) |
| `pca100.sct_harmony_clusters` | 13 | Harmony (SCT assay) |
| `pca100.sct_cca_clusters` | 12 | CCA (SCT assay) |

### 8.2 Legacy Dimensional Reductions
- **Standard PCA**: `pca`, `pca50`, `pca100`, `pca200` (computed on RNA / SCT)
- **Unintegrated UMAP**: `umap.unintegrated`
- **Batch Corrected Coordinates (Harmony)**: `pca100_harmony`, `umap.pca100_harmony`, `pca100.sct_harmony`, `umap.pca100.sct_harmony`
- **Batch Corrected Coordinates (CCA)**: `pca100_cca`, `umap.pca100_cca`, `pca100.sct_cca`, `umap.pca100.sct_cca`
- **Batch Corrected Coordinates (RPCA)**: `pca100_rpca`, `umap.pca100_rpca`
- **Batch Corrected Coordinates (MNN)**: `pca100_mnn`, `umap.pca100_mnn`

### 8.3 Integration Methods Inventory
- **Harmony**: Yes, stored in Harmony reductions and cluster columns.
- **CCA**: Yes, stored in CCA reductions and cluster columns.
- **RPCA**: Yes, stored in RPCA reductions and cluster columns.
- **MNN**: Yes, stored in `mnn.reconstructed` assay and MNN reductions.
- **SCT-derived**: Yes, stored in `SCT` assay, `nCount_SCT`, `nFeature_SCT` metadata, and SCT clustering columns.

> [!IMPORTANT]
> **Policy and Baseline Isolation**: These legacy clusters, reductions, and assays will **NOT** influence or contaminate the Phase 1 independent reprocessing. They are documented here purely for historical inventory. All Phase 1 downstream pipelines must bypass these features, starting from the split raw counts of the `RNA` assay.

## 9. Workflow Readiness

- **Independent Dataset Extraction**: **READY**. The object contains split raw counts for 4 libraries (`counts.MPNST_1` through `counts.MPNST_4`) matching the `sample_id` column.
- **QC (Quality Control)**: **READY**. Standard QC columns (`nCount_RNA`, `nFeature_RNA`, `percent.mt`) are already present. The raw counts are intact for calculating additional metrics.
- **Doublet Detection**: **READY**. Cells are annotated with barcode names and split layers allow for independent sample-by-sample doublet estimation (e.g. via `scDblFinder`).
- **Normalization**: **READY**. Raw RNA count layers are present and ready for log-normalization or SCTransform.
- **PCA**: **READY**. High-quality variable features can be computed from the normalized counts, enabling standard PCA.
- **Clustering**: **READY**. Neighborhood graphs and SNN can be constructed from new PCA space.

### Scientific and Computational Risks
1. **Risk of Legacy Contamination**: Having prior integrations and clusters in the same RDS might lead to automated tools loading them by default. This risk is mitigated by writing a splitting script that extracts only the raw counts (`counts` layers) and basic cell info, dropping all reductions and legacy clustering.
2. **Seurat v5 Class Differences**: The `RNA` assay is of class `Assay5`, which manages split layers. Downstream workflows must be Seurat v5 layer-aware (e.g. running workflows layer-by-layer or using `JoinLayers` where appropriate).

## 10. Computational Performance

### SLURM Run Execution Metrics (Job ID: 19160528)

| Metric | Value |
| --- | --- |
| **Load Time** | 45.35 seconds |
| **SLURM Resources Requested** | 8 CPUs, 128 GB Memory, Partition: ihc |
| **Elapsed Runtime** | 1 minute 32 seconds |
| **ExitCode** | 0:0 |
| **State** | COMPLETED |
| **MaxRSS (Peak Memory)** | 9627268K (~9.18 GB) |
| **CPU Usage** | mean_load 54.72%, cpu_time 42.23s |

### Resource Appropriateness Analysis
- The requested memory of **128 GB** was significantly over-allocated for this audit run, which had a peak memory footprint of **~9.18 GB**. This represents an efficiency of only ~7.2%.
- The CPU allocation of **8 cores** was also underutilized because loading RDS files and basic auditing are primarily single-threaded in R. 
- While this over-allocation was appropriate for a safety-first initial exploration to guarantee zero out-of-memory (OOM) failures, it is inefficient for recurring production workflows.

### Recommended M2 SLURM Job Resources
- For the upcoming dataset extraction and pre-filter QC stage (Milestone M2), we recommend the following resource allocations:
  - **CPUs**: 4 cores (allowing parallel file writing for the 4 split datasets)
  - **Memory**: 24 GB (safe headroom for reading the 10 GB object and writing split RDS files)
  - **Time Limit**: 1 hour
  - **Partition / Account**: `ihc` / `ihc`

## 11. Machine-readable Outputs

- **Object Inventory JSON**: [object_inventory.json](file://reports/object_inventory.json)
- **Object Inventory TSV**: [object_inventory.tsv](file://reports/object_inventory.tsv)

Both files have been successfully expanded to include all audited fields, including memory usage, feature statistics, ribosomal/mitochondrial gene counts, cell frequencies, and recommended metadata usages.

## 12. Researcher Decisions Required

Before starting Milestone 2 (M2) and beginning downstream preprocessing, the researcher must approve the following decisions:

1. **Confirm standardizing on `sample_id`** as the primary dataset identifier for splitting the object.
2. **Approve discarding all legacy reductions and cluster labels** from the downstream single-cell pipeline to maintain isolation.
3. **Approve sample-specific QC thresholds** (e.g. customizing min/max features and MT percentage cutoffs independently for MPNST_1 through MPNST_4) rather than applying a single uniform cutoff.
4. **Approve `scDblFinder` as the doublet detection method** with a default expected doublet rate of 7.5% per sample.
5. **Confirm if external clinical metadata** (such as patient demographics, treatment status, or anatomical location) should be integrated with these sample IDs at this stage.
