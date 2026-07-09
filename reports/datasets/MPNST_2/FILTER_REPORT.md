# Milestone 3 QC and Doublet Filtering Report - MPNST_2
*Generated on: 2026-07-09 08:07:50*

## 1. Observed Results

### Cell Count Metrics Summary Table

| Stage | Cell Count | Percentage of Raw |
| --- | --- | --- |
| **Raw Input (M2)** | 2830 | 100.00% |
| **Retained Singlets (M3)** | 1502 | 53.07% |
| **Excluded Cells (Total)** | 1328 | 46.93% |

### QC Filtering Exclusions Breakdown Table

| Filter Criterion | Threshold Applied | Cells Excluded (Independently) | % Excluded (Independently) |
| --- | --- | --- | --- |
| **Min Features (nFeature_RNA)** | `< 200` | 0 | 0.00% |
| **Min Counts (nCount_RNA)** | `< 500` | 0 | 0.00% |
| **Max percent.mt** | `> 10.0%` | 487 | 17.21% |
| **Max percent.ribo** | `> 20.0%` | 796 | 28.13% |
| **Doublet Detection** | `scDblFinder` | 226 | 7.99% |

### Doublet Assessment Details

- **Doublet Detection Method**: `scDblFinder` (Version: `1.20.2`)
- **Configured Expected Doublet Rate**: `7.50%`
- **Observed Predicted Doublet Rate**: `7.99%` (226 doublets / 2830 cells)
- **Doublet Detection Runtime**: `39.08 seconds`

### High-Feature Doublet Enrichment Check

- **Median Genes in Predicted Singlets**: 2290
- **Median Genes in Predicted Doublets**: 3996
- **High-Feature Threshold (90th percentile)**: 4541 genes
- **Doublet percentage among high-feature cells**: 21.55% (61 doublets / 283 cells)

---

## 2. Interpretation

1. **Mitochondrial Filtering**: The maximum mitochondrial threshold of 10% was applied. 
For MPNST_2, the 10% MT threshold excluded 487 cells (17.21%), effectively removing cells with signs of stress or lysis.
2. **Ribosomal Filtering**: The maximum ribosomal threshold of 20% effectively targeted outlier cells. Ribosomal transcript expression is related to high translation rates and can indicate technical noise.
3. **Doublet Detection**: Doublets represent technical artifacts where two cells are captured in a single droplet.
Doublet detection using `scDblFinder` identified 226 doublets. We observed that predicted doublets have a median of 3996 genes compared to 2290 genes for singlets, confirming that doublets are indeed enriched for high transcript complexity.
Furthermore, 21.55% of high-feature cells (genes > 4541) were classified as doublets, showing that transcript complexity is highly correlated with doublet rate, but that not all high-feature cells are doublets, justifying the preservation of high-feature singlets.

---

## 3. Recommendations

1. **Proceeding to Milestone 4 (M4) Normalization**: The filtered object is clean, doublet-free, and contains only validated cells. We recommend proceeding to M4 normalization using the default SCTransform workflow.
2. **Special Normalization Handling**: 
The dataset behaves normally and standard SCTransform regression (e.g. regressing out percent.mt and nCount_RNA) is recommended.

## 4. Provenance

- **Input Raw File**: `results/datasets/MPNST_2/MPNST_2_raw.rds`
- **Input Checksum**: `cd15a88b7f70a6492a1f05c94725cd8c`
- **Git Commit Hash**: `7bf25a45b4bf55a889e78a400aed7490183c3edf`
- **Software versions**: R version 4.4.3 (2025-02-28), Seurat 5.1.0
