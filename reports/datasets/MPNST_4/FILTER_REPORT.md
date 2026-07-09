# Milestone 3 QC and Doublet Filtering Report - MPNST_4
*Generated on: 2026-07-09 09:30:06*

## 1. Observed Results

### Cell Count Metrics Summary Table

| Stage | Cell Count | Percentage of Raw |
| --- | --- | --- |
| **Raw Input (M2)** | 7811 | 100.00% |
| **Retained Singlets (M3)** | 3378 | 43.25% |
| **Excluded Cells (Total)** | 4433 | 56.75% |

### QC Filtering Exclusions Breakdown Table

| Filter Criterion | Threshold Applied | Cells Excluded (Independently) | % Excluded (Independently) |
| --- | --- | --- | --- |
| **Min Features (nFeature_RNA)** | `< 200` | 0 | 0.00% |
| **Min Counts (nCount_RNA)** | `< 500` | 0 | 0.00% |
| **Max percent.mt** | `> 10.0%` | 2343 | 30.00% |
| **Max percent.ribo** | `> 20.0%` | 1793 | 22.95% |
| **Doublet Detection** | `scDblFinder` | 717 | 9.18% |

### Doublet Assessment Details

- **Doublet Detection Method**: `scDblFinder` (Version: `1.20.2`)
- **Configured Expected Doublet Rate**: `7.50%`
- **Observed Predicted Doublet Rate**: `9.18%` (717 doublets / 7811 cells)
- **Doublet Detection Runtime**: `72.85 seconds`

### High-Feature Doublet Enrichment Check

- **Median Genes in Predicted Singlets**: 2093
- **Median Genes in Predicted Doublets**: 3323
- **High-Feature Threshold (90th percentile)**: 3443 genes
- **Doublet percentage among high-feature cells**: 40.13% (313 doublets / 780 cells)

---

## 2. Interpretation

1. **Mitochondrial Filtering**: The maximum mitochondrial threshold of 10% was applied. 
For MPNST_4, the 10% MT threshold excluded 2343 cells (30.00%), effectively removing cells with signs of stress or lysis.
2. **Ribosomal Filtering**: The maximum ribosomal threshold of 20% effectively targeted outlier cells. Ribosomal transcript expression is related to high translation rates and can indicate technical noise.
3. **Doublet Detection**: Doublets represent technical artifacts where two cells are captured in a single droplet.
Doublet detection using `scDblFinder` identified 717 doublets. We observed that predicted doublets have a median of 3323 genes compared to 2093 genes for singlets, confirming that doublets are indeed enriched for high transcript complexity.
Furthermore, 40.13% of high-feature cells (genes > 3443) were classified as doublets, showing that transcript complexity is highly correlated with doublet rate, but that not all high-feature cells are doublets, justifying the preservation of high-feature singlets.

---

## 3. Recommendations

1. **Proceeding to Milestone 4 (M4) Normalization**: The filtered object is clean, doublet-free, and contains only validated cells. We recommend proceeding to M4 normalization using the default SCTransform workflow.
2. **Special Normalization Handling**: 
The dataset behaves normally and standard SCTransform regression (e.g. regressing out percent.mt and nCount_RNA) is recommended.

## 4. Provenance

- **Input Raw File**: `results/datasets/MPNST_4/MPNST_4_raw.rds`
- **Input Checksum**: `f32811cf69595dcab9c5f5067efccbaa`
- **Git Commit Hash**: `4e4c2aef7c9faf4a5adbcd0d0090f39d99737ac2`
- **Software versions**: R version 4.4.3 (2025-02-28), Seurat 5.1.0
