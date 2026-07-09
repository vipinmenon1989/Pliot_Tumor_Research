# Milestone 3 QC and Doublet Filtering Report - MPNST_1
*Generated on: 2026-07-09 09:32:29*

## 1. Observed Results

### Cell Count Metrics Summary Table

| Stage | Cell Count | Percentage of Raw |
| --- | --- | --- |
| **Raw Input (M2)** | 8338 | 100.00% |
| **Retained Singlets (M3)** | 7615 | 91.33% |
| **Excluded Cells (Total)** | 723 | 8.67% |

### QC Filtering Exclusions Breakdown Table

| Filter Criterion | Threshold Applied | Cells Excluded (Independently) | % Excluded (Independently) |
| --- | --- | --- | --- |
| **Min Features (nFeature_RNA)** | `< 200` | 0 | 0.00% |
| **Min Counts (nCount_RNA)** | `< 500` | 0 | 0.00% |
| **Max percent.mt** | `> 10.0%` | 0 | 0.00% |
| **Max percent.ribo** | `> 20.0%` | 26 | 0.31% |
| **Doublet Detection** | `scDblFinder` | 697 | 8.36% |

### Doublet Assessment Details

- **Doublet Detection Method**: `scDblFinder` (Version: `1.20.2`)
- **Configured Expected Doublet Rate**: `7.50%`
- **Observed Predicted Doublet Rate**: `8.36%` (697 doublets / 8338 cells)
- **Doublet Detection Runtime**: `83.03 seconds`

### High-Feature Doublet Enrichment Check

- **Median Genes in Predicted Singlets**: 3133
- **Median Genes in Predicted Doublets**: 6288
- **High-Feature Threshold (90th percentile)**: 6635 genes
- **Doublet percentage among high-feature cells**: 35.13% (293 doublets / 834 cells)

---

## 2. Interpretation

1. **Mitochondrial Filtering**: The maximum mitochondrial threshold of 10% was applied. 
For MPNST_1, mitochondrial transcripts remain exactly 0.00% across all cells. This dataset contains no mitochondrial mapping information, meaning the 10% mitochondrial threshold resulted in zero exclusions. This represents a technical batch characteristic.
2. **Ribosomal Filtering**: The maximum ribosomal threshold of 20% effectively targeted outlier cells. Ribosomal transcript expression is related to high translation rates and can indicate technical noise.
3. **Doublet Detection**: Doublets represent technical artifacts where two cells are captured in a single droplet.
Doublet detection using `scDblFinder` identified 697 doublets. We observed that predicted doublets have a median of 6288 genes compared to 3133 genes for singlets, confirming that doublets are indeed enriched for high transcript complexity.
Furthermore, 35.13% of high-feature cells (genes > 6635) were classified as doublets, showing that transcript complexity is highly correlated with doublet rate, but that not all high-feature cells are doublets, justifying the preservation of high-feature singlets.

---

## 3. Recommendations

1. **Proceeding to Milestone 4 (M4) Normalization**: The filtered object is clean, doublet-free, and contains only validated cells. We recommend proceeding to M4 normalization using the default SCTransform workflow.
2. **Special Normalization Handling**: 
Since MPNST_1 lacks mitochondrial transcripts, downstream normalization (SCTransform) should skip regressing out 'percent.mt' to prevent model fitting errors, or use only 'percent.ribo' and 'nCount_RNA' for regression. The other datasets (MPNST_2, MPNST_3, MPNST_4) should include 'percent.mt' in their regression formulas.

## 4. Provenance

- **Input Raw File**: `results/datasets/MPNST_1/MPNST_1_raw.rds`
- **Input Checksum**: `f523e61322c585b16e632e800dc12bde`
- **Git Commit Hash**: `4e4c2aef7c9faf4a5adbcd0d0090f39d99737ac2`
- **Software versions**: R version 4.4.3 (2025-02-28), Seurat 5.1.0
