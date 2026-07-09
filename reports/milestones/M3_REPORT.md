# Milestone 3 (M3) Execution Report — QC Filtering and Doublet Assessment

*Generated on: 2026-07-09 08:15:00 EDT*

## 1. Execution Summary
- **Authorized Milestone**: Milestone 3 (M3) — QC Filtering and Doublet Assessment
- **Input Objects**: 4 dataset-specific raw Seurat objects (M2 outputs)
- **Output Objects Created**: 4 high-quality filtered Seurat objects (retaining only singlets that pass QC)
- **QC Thresholds Applied**:
  - `min_features` = 200
  - `min_counts` = 500
  - `max_percent_mt` = 10.0%
  - `max_percent_ribo` = 20.0%
  - `max_features` = None (unlimited)
  - `max_counts` = None (unlimited)
- **HPC Execution Details**:
  - **SLURM Job ID**: 19161422
  - **Node**: `ihc-grid-1-1-1`
  - **State**: COMPLETED (ExitCode 0:0)
  - **Elapsed Time**: 3 minutes 24 seconds
  - **Requested Resources**: 8 CPUs, 64 GB RAM, 12 hours
  - **Peak Memory (MaxRSS)**: 38,676,744K (~36.88 GB) — Highly efficient execution.
- **Git Commit Hash**: `7bf25a45b4bf55a889e78a400aed7490183c3edf`

---

## 2. Cross-Dataset Comparison Summary

The following table summarizes the QC filtering and doublet detection statistics across all four constituent datasets:

| Dataset ID | Cells Raw (M2) | Cells Filtered (M3) | Cells Removed | Percent Removed | Predicted Doublets | Doublet % | Median Genes (Singlets) | Median Genes (Doublets) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **MPNST_1** | 8,338 | 7,615 | 723 | 8.67% | 697 | 8.36% | 3,254 | 5,420 |
| **MPNST_2** | 2,830 | 1,502 | 1,328 | 46.93% | 226 | 7.99% | 2,347 | 4,213 |
| **MPNST_3** | 3,682 | 1,904 | 1,778 | 48.29% | 327 | 8.88% | 1,489 | 3,115 |
| **MPNST_4** | 7,811 | 3,378 | 4,433 | 56.75% | 717 | 9.18% | 2,058 | 3,923 |

### QC Exclusion Drivers (Percent cells removed independently):
- **`MPNST_1`**: Min Features: 0.00% | Min Counts: 0.00% | Max MT: 0.00% | Max Ribo: 0.31% | Doublet: 8.36%
- **`MPNST_2`**: Min Features: 0.00% | Min Counts: 0.00% | Max MT: 17.21% | Max Ribo: 28.13% | Doublet: 7.99%
- **`MPNST_3`**: Min Features: 0.00% | Min Counts: 0.00% | Max MT: 5.19% | Max Ribo: 37.10% | Doublet: 8.88%
- **`MPNST_4`**: Min Features: 0.00% | Min Counts: 0.00% | Max MT: 30.00% | Max Ribo: 22.95% | Doublet: 9.18%

---

## 3. Doublet Assessment & High-Feature Enrichment Analysis

### Selected Method and Software Version:
- **Software**: `scDblFinder` (Bioconductor package, version `1.20.2`)
- **Justification**: Highly standard, robust, compatible with Seurat v5 count matrices via temporary conversion to `SingleCellExperiment` (after joining raw layers of the RNA assay), and preserves deterministic reproducibility.
- **Parameters**: Expected doublet rate `0.075` (7.5%), using default parameters (Poisson doublet simulation, PCA and kNN classification).

### High-Feature Enrichment Evidence:
We analyzed the enrichment of predicted doublets in the high-complexity cells (top 10% highest gene counts):
- **MPNST_1**: 35.13% of the highest-feature cells are predicted doublets (293 doublets / 834 cells).
- **MPNST_2**: 21.55% of the highest-feature cells are predicted doublets (61 doublets / 283 cells).
- **MPNST_3**: 18.70% of the highest-feature cells are predicted doublets (69 doublets / 369 cells).
- **MPNST_4**: 40.13% of the highest-feature cells are predicted doublets (313 doublets / 780 cells).

### Scientific Interpretation:
Predicted doublets are significantly enriched in the high-feature population (median gene counts for doublets are roughly 1.7x to 2.1x higher than singlets). However, the majority of high-feature cells (60% to 81%) are classified as **singlets**. This demonstrates that:
1. High transcript complexity alone is not a sufficient indicator of doublet status; actively transcribing tumor cells can naturally express many features.
2. Applying a hard maximum feature threshold would have inappropriately discarded a large number of high-complexity singlet cells, resulting in a loss of valuable biological signal.
3. The use of a formal doublet classification model (`scDblFinder`) is scientifically justified and necessary to separate technical artifacts from biological complexity.

---

## 4. Key Scientific and Technical Concerns

1. **Massive Cell Depletion in MPNST_2, MPNST_3, and MPNST_4**:
   - Applying uniform QC thresholds (`max_percent_mt = 10%` and `max_percent_ribo = 20%`) resulted in extremely high exclusion rates: **46.93%** in MPNST_2, **48.29%** in MPNST_3, and **56.75%** in MPNST_4!
   - In MPNST_3, the **37.10% ribosomal exclusion** is particularly severe, as is the **30.00% mitochondrial exclusion** in MPNST_4.
   - Such high exclusion rates risk losing entire cell populations, depleting biologically relevant tumor subclones, and introducing technical bias.
2. **Technical Discrepancy of MPNST_1**:
   - `MPNST_1` has exactly **0.00% mitochondrial transcripts** in all cells, meaning it suffered 0.00% MT-based exclusion. This represents a major batch effect or mapping difference that must be handled carefully in downstream integration.

---

## 5. Recommendations for Milestone 4 (M4) Normalization

1. **Approved Datasets for Normalization**:
   - All four filtered datasets (`MPNST_1` to `MPNST_4`) are technically ready for normalization.
2. **SCTransform Formula Adjustments**:
   - **`MPNST_1`**: Since mitochondrial transcripts are absent, SCTransform **must NOT** attempt to regress out `percent.mt` (this would cause numerical singularity errors or model failure). Regression should only include `nCount_RNA` and `percent.ribo`.
   - **`MPNST_2`, `MPNST_3`, `MPNST_4`**: Standard SCTransform is recommended, regressing out both `percent.mt` and `nCount_RNA`.
3. **Threshold Re-evaluation**:
   - Before freezing M3 and proceeding, the researcher should consider whether the mitochondrial threshold for `MPNST_4` (currently 10%) should be relaxed (e.g. to 15% or 20%) to avoid losing 30% of the dataset, and whether the ribosomal threshold for `MPNST_3` (currently 20%) should be relaxed (e.g. to 25% or 30%) to retain more cells.

---

## 6. Generated Artifacts
- **Filtered RDS Objects**:
  - `results/datasets/MPNST_1/MPNST_1_filtered.rds` (MD5: `39542f88f51f16c7858368736d322090`)
  - `results/datasets/MPNST_2/MPNST_2_filtered.rds` (MD5: `23fbba71da436673dabbdf285f6aaac6`)
  - `results/datasets/MPNST_3/MPNST_3_filtered.rds` (MD5: `3d9024896629a60d5245f8f1e2839721`)
  - `results/datasets/MPNST_4/MPNST_4_filtered.rds` (MD5: `c02c0996f78167196c671ce05c1f120e`)
- **Dataset-Specific Filter Reports**:
  - [MPNST_1 FILTER_REPORT.md](file://reports/datasets/MPNST_1/FILTER_REPORT.md)
  - [MPNST_2 FILTER_REPORT.md](file://reports/datasets/MPNST_2/FILTER_REPORT.md)
  - [MPNST_3 FILTER_REPORT.md](file://reports/datasets/MPNST_3/FILTER_REPORT.md)
  - [MPNST_4 FILTER_REPORT.md](file://reports/datasets/MPNST_4/FILTER_REPORT.md)
- **Publication-Quality Post-filter Figures** (violins, scatter, density, histograms, doublet status, filtering summaries in PDF/PNG).
- **Global Figure Index**: [reports/FIGURE_INDEX.tsv](file://reports/FIGURE_INDEX.tsv) (fully updated with M3 figures).
- **Milestone Progress Log**: [PROGRESS.md](file://PROGRESS.md) (M3 marked complete).
- **Changelog**: [CHANGELOG.md](file://CHANGELOG.md) (M2 frozen, M3 recorded).
