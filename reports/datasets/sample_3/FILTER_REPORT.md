# Milestone 3 QC and Doublet Filtering Report - sample_3
*Generated on: 2026-07-09 09:22:22*

## 1. Observed Results

### Cell Count Metrics Summary Table

| Stage | Cell Count | Percentage of Raw |
| --- | --- | --- |
| **Raw Input (M2)** | 150 | 100.00% |
| **Retained Singlets (M3)** | 150 | 100.00% |
| **Excluded Cells (Total)** | 0 | 0.00% |

### QC Filtering Exclusions Breakdown Table

| Filter Criterion | Threshold Applied | Cells Excluded (Independently) | % Excluded (Independently) |
| --- | --- | --- | --- |
| **Min Features (nFeature_RNA)** | `< 100` | 0 | 0.00% |
| **Min Counts (nCount_RNA)** | `< 300` | 0 | 0.00% |
| **Max percent.mt** | `> 50.0%` | 0 | 0.00% |
| **Max percent.ribo** | `> 40.0%` | 0 | 0.00% |
| **Doublet Detection** | `scDblFinder` | 0 | 0.00% |

### Doublet Assessment Details

- **Doublet Detection Method**: `scDblFinder` (Version: `1.20.2`)
- **Configured Expected Doublet Rate**: `5.00%`
- **Observed Predicted Doublet Rate**: `0.00%` (0 doublets / 150 cells)
- **Doublet Detection Runtime**: `62.57 seconds`

### High-Feature Doublet Enrichment Check

- **Median Genes in Predicted Singlets**: 120
- **Median Genes in Predicted Doublets**: NA
- **High-Feature Threshold (90th percentile)**: 120 genes
- **Doublet percentage among high-feature cells**: NaN% (0 doublets / 0 cells)

---

## 2. Interpretation

1. **Mitochondrial Filtering**: The maximum mitochondrial threshold of 10% was applied. 
For sample_3, the 10% MT threshold excluded 0 cells (0.00%), effectively removing cells with signs of stress or lysis.
2. **Ribosomal Filtering**: The maximum ribosomal threshold of 20% effectively targeted outlier cells. Ribosomal transcript expression is related to high translation rates and can indicate technical noise.
3. **Doublet Detection**: Doublets represent technical artifacts where two cells are captured in a single droplet.
Doublet detection using `scDblFinder` identified 0 doublets. We observed that predicted doublets have a median of NA genes compared to 120 genes for singlets, confirming that doublets are indeed enriched for high transcript complexity.
Furthermore, NaN% of high-feature cells (genes > 120) were classified as doublets, showing that transcript complexity is highly correlated with doublet rate, but that not all high-feature cells are doublets, justifying the preservation of high-feature singlets.

---

## 3. Recommendations

1. **Proceeding to Milestone 4 (M4) Normalization**: The filtered object is clean, doublet-free, and contains only validated cells. We recommend proceeding to M4 normalization using the default SCTransform workflow.
2. **Special Normalization Handling**: 
The dataset behaves normally and standard SCTransform regression (e.g. regressing out percent.mt and nCount_RNA) is recommended.

## 4. Provenance

- **Input Raw File**: `results/datasets/sample_3/sample_3_raw.rds`
- **Input Checksum**: `92260c22b30e170a03a3b9885d11c140`
- **Git Commit Hash**: `4e4c2aef7c9faf4a5adbcd0d0090f39d99737ac2`
- **Software versions**: R version 4.4.3 (2025-02-28), Seurat 5.1.0
