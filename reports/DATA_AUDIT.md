# Scientific and Computational Audit of Seurat Object
*Generated on: 2026-07-08 19:27:53 EDT*

## 1. Executive Summary
- **Object Class**: Seurat
- **Seurat Version**: 5.1.0
- **File Size on Disk**: 7.23 GB
- **Dimensions**: 29708 features (genes) x 22661 cells
- **Default Assay**: SCT
- **Raw RNA Counts Detected**: YES
- **Independent Reprocessing Possible**: YES
- **Recommended Dataset Identifier**: `sample_id`

---

## 2. Observed Facts

### 2.1 Object Architecture & Assays
The object contains the following assays and dimensional reductions:
| Assay | Class | Contents / Slots / Layers | Default? |
| --- | --- | --- | --- |
| RNA | Assay5 | counts.MPNST_1, counts.MPNST_2, counts.MPNST_3, counts.MPNST_4, data.MPNST_1, data.MPNST_2, data.MPNST_3, data.MPNST_4, scale.data | NO |
| mnn.reconstructed | Assay | data | NO |
| SCT | SCTAssay | counts, data, scale.data | YES |

### 2.2 Dimensional Reductions
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

### 2.3 Cell Networks & Graphs
- RNA_nn
- RNA_snn
- SCT_nn
- SCT_snn

### 2.4 Active Identities
Active identity class: `factor`
Total unique identities: 12
Identity levels (top 10):
- 1
- 3
- 4
- 0
- 7
- 2
- 5
- 6
- 9
- 8

### 2.5 Metadata Schema Audit
Below is the complete inventory of metadata columns present in the object:
| Column Name | Data Type | Cardinality | Missing Values | Biological / Technical Meaning |
| --- | --- | --- | --- | --- |
| orig.ident | character | 4 | 0 | Initial sample identifier from library setup |
| nCount_RNA | numeric | 12933 | 0 | Total RNA UMI count per cell (library depth) |
| nFeature_RNA | integer | 6289 | 0 | Number of unique genes detected in RNA assay per cell |
| sample_id | character | 4 | 0 | Consistent sample metadata identifier |
| percent.mt | numeric | 14130 | 0 | Mitochondrial transcript percentage (quality control metric) |
| orig.anno | character | 16 | 0 | Legacy cell-type annotations from original analysis |
| unintegrated_clusters | factor | 16 | 0 | Legacy cluster assignments |
| seurat_clusters | factor | 12 | 0 | Legacy cluster assignments |
| pca100_harmony_clusters | factor | 16 | 0 | Legacy cluster assignments |
| pca100_cca_clusters | factor | 10 | 0 | Legacy cluster assignments |
| pca100_mnn_clusters | factor | 13 | 0 | Legacy cluster assignments |
| pca100_rpca_clusters | factor | 17 | 0 | Legacy cluster assignments |
| nCount_SCT | numeric | 5887 | 0 | Total SCTransformed UMI count per cell |
| nFeature_SCT | integer | 4173 | 0 | Number of unique genes detected in SCT assay per cell |
| pca100.sct_harmony_clusters | factor | 13 | 0 | Legacy cluster assignments |
| pca100.sct_cca_clusters | factor | 12 | 0 | Legacy cluster assignments |

### 2.6 Dataset and Sample Cross-Tabulation
Evaluation of `sample_id` and `orig.ident` distributions across cells:
| sample_id | orig.ident | Cell Count |
| --- | --- | --- |
| MPNST_1 | MPNST_1 | 8338 |
| MPNST_2 | MPNST_2 | 2830 |
| MPNST_3 | MPNST_3 | 3682 |
| MPNST_4 | MPNST_4 | 7811 |

### 2.7 Cell Statistics Summary
- **Total Cells**: 22661
- **Total Genes**: 29708

#### Cells per sample_id:
| sample_id | Cells |
| --- | --- |
| MPNST_1 | 8338 |
| MPNST_2 | 2830 |
| MPNST_3 | 3682 |
| MPNST_4 | 7811 |

#### Cells per orig.ident:
| orig.ident | Cells |
| --- | --- |
| MPNST_1 | 8338 |
| MPNST_2 | 2830 |
| MPNST_3 | 3682 |
| MPNST_4 | 7811 |

### 2.8 Legacy Analysis Inventory
The following metadata columns are identified as prior analysis artifacts and **must not** be used for Phase 1 decisions:
| Legacy Column | Exists in Object? | Number of Classes |
| --- | --- | --- |
| seurat_clusters | YES | 12 |
| unintegrated_clusters | YES | 16 |
| pca100_harmony_clusters | YES | 16 |
| pca100_cca_clusters | YES | 10 |
| pca100_mnn_clusters | YES | 13 |
| pca100_rpca_clusters | YES | 17 |
| pca100.sct_harmony_clusters | YES | 13 |
| pca100.sct_cca_clusters | YES | 12 |

---

## 3. Inferences & Interpretations

### 3.1 Scientific Reprocessing Feasibility
Raw RNA UMI counts exist in the RNA assay counts slot/layer, and cell/feature dimensions are intact. Therefore, complete independent reprocessing is scientifically and computationally possible.

### 3.2 Metadata Column Meaning
- `sample_id`: Represents the distinct sequencing libraries or tumor samples. Cardinality shows how many independent samples are integrated.
- `orig.ident`: Represents the library ID or run ID. This appears to map 1-to-1 with sample_id in some cells, but must be cross-referenced.
- `percent.mt` & `percent.ribo`: Represent cell quality. High mitochondrial fraction usually indicates dying or damaged cells.
- `orig.anno`: Legacy cell-type assignments, possibly done manually or via automated classifiers in prior iterations.

---

## 4. Unknowns & Technical Gaps
- **Sample details**: Demographics, anatomical site of tumor collection, and batch details are missing from the Seurat object itself.
- **Library preparation batch**: Whether multiple libraries were processed on the same sequencing run or across different chemistry versions is unknown.
- **Patient demographics**: Patient clinical metadata is not fully annotated in the object, besides what is inferable from candidate patient variables.

---

## 5. Scientific and Computational Risks

### 5.1 Scientific Risks
- **Pseudoreplication**: Biological replicates must be handled at the patient level rather than cell level. Cell counts per patient vary significantly, which can bias findings if not handled correctly.
- **Legacy Contamination**: The presence of prior integration coordinates and clusters might tempt downstream analysts to shortcut the pipeline, violating the independent-processing directive.

### 5.2 Computational Risks
- **Memory Footprint**: The real RDS is very large (approx. 7.2 GB on disk). Deserializing this in R requires substantial RAM. All processing must be constrained to the SLURM HPC queue with proper allocations.
- **Seurat Versioning**: The object version must be compatible with the Seurat package version in `R_env`. If version mismatches occur, certain assays or functions may throw deprecation errors.

---

## 6. Recommendations & Action Plan
1. **Primary Dataset Identifier**: Standardize on `sample_id` as the primary grouping variable for splitting and preprocessing the independent datasets.
2. **Independent Preprocessing**: Split the Seurat object into sample-specific Seurat objects based on `sample_id` and run independent QC filtering on each.
3. **Discard Legacy Coordinates**: Explicitly ignore and delete all legacy cluster annotations and dimensional reductions before running new Phase 1 workflows to prevent bias.

---

## 7. Researcher Decisions Required
1. **Confirm primary dataset splitter**: Verify that splitting by `sample_id` is correct.
2. **Patient metadata mapping**: Confirm if additional patient-level clinical metadata should be integrated at this stage.
3. **QC Threshold alignment**: Approve the use of configurable, sample-specific QC thresholds rather than a single uniform cutoff.
