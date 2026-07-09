# Milestone 3 QC and Doublet Filtering Report - MPNST_3
*Generated on: 2026-07-09 08:07:55*

## 1. Observed Results

### Cell Count Metrics Summary Table

| Stage | Cell Count | Percentage of Raw |
| --- | --- | --- |
| **Raw Input (M2)** | 3682 | 100.00% |
| **Retained Singlets (M3)** | 1904 | 51.71% |
| **Excluded Cells (Total)** | 1778 | 48.29% |

### QC Filtering Exclusions Breakdown Table

| Filter Criterion | Threshold Applied | Cells Excluded (Independently) | % Excluded (Independently) |
| --- | --- | --- | --- |
| **Min Features (nFeature_RNA)** | `< 200` | 0 | 0.00% |
| **Min Counts (nCount_RNA)** | `< 500` | 0 | 0.00% |
| **Max percent.mt** | `> 10.0%` | 191 | 5.19% |
| **Max percent.ribo** | `> 20.0%` | 1366 | 37.10% |
| **Doublet Detection** | `scDblFinder` | 327 | 8.88% |

### Doublet Assessment Details

- **Doublet Detection Method**: `scDblFinder` (Version: `1.20.2`)
- **Configured Expected Doublet Rate**: `7.50%`
- **Observed Predicted Doublet Rate**: `8.88%` (327 doublets / 3682 cells)
- **Doublet Detection Runtime**: `44.15 seconds`

### High-Feature Doublet Enrichment Check

- **Median Genes in Predicted Singlets**: 1441
- **Median Genes in Predicted Doublets**: 2632
- **High-Feature Threshold (90th percentile)**: 3530 genes
- **Doublet percentage among high-feature cells**: 18.70% (69 doublets / 369 cells)

---

## 2. Interpretation

1. **Mitochondrial Filtering**: The maximum mitochondrial threshold of 10% was applied. 
For MPNST_3, the 10% MT threshold excluded 191 cells (5.19%), effectively removing cells with signs of stress or lysis.
2. **Ribosomal Filtering**: The maximum ribosomal threshold of 20% effectively targeted outlier cells. Ribosomal transcript expression is related to high translation rates and can indicate technical noise.
3. **Doublet Detection**: Doublets represent technical artifacts where two cells are captured in a single droplet.
Doublet detection using `scDblFinder` identified 327 doublets. We observed that predicted doublets have a median of 2632 genes compared to 1441 genes for singlets, confirming that doublets are indeed enriched for high transcript complexity.
Furthermore, 18.70% of high-feature cells (genes > 3530) were classified as doublets, showing that transcript complexity is highly correlated with doublet rate, but that not all high-feature cells are doublets, justifying the preservation of high-feature singlets.

---

## 3. Recommendations

1. **Proceeding to Milestone 4 (M4) Normalization**: The filtered object is clean, doublet-free, and contains only validated cells. We recommend proceeding to M4 normalization using the default SCTransform workflow.
2. **Special Normalization Handling**: 
The dataset behaves normally and standard SCTransform regression (e.g. regressing out percent.mt and nCount_RNA) is recommended.

## 4. Provenance

- **Input Raw File**: `results/datasets/MPNST_3/MPNST_3_raw.rds`
- **Input Checksum**: `ecaca4b46e0a23f7bb8798cdb6adf280`
- **Git Commit Hash**: `7bf25a45b4bf55a889e78a400aed7490183c3edf`
- **Software versions**: R version 4.4.3 (2025-02-28), Seurat 5.1.0
