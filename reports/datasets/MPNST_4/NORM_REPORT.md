# Milestone 4 Normalization and Variable Feature Report - MPNST_4
*Generated on: 2026-07-09 19:59:28*

## 1. Parameters & Configuration

- **Input Filtered Object**: `results/datasets/MPNST_4/MPNST_4_filtered_specific.rds`
- **Output Normalized Object**: `results/datasets/MPNST_4/MPNST_4_normalized.rds`
- **Normalization Method**: `SCTransform`
- **Requested Variable Features**: `3000`
- **Mitochondrial Regression**: `YES (percent.mt)`
- **Seed Used**: `42`
- **Execution Time**: `50.51 seconds`

## 2. Dataset Characteristics

- **Total Cells**: `6877`
- **Total Raw Genes**: `29708`

## 3. High Variance Features

Top 20 highly variable features selected under `SCTransform` (sorted by residual_variance):



## 4. Diagnostics & Visualizations
The following diagnostic plots were generated to assess normalization quality:
- [Scatter Plot (HVF)](file:///reports/datasets/MPNST_4/var_features_scatter.png)
- [Distribution Plot (HVF Metric)](file:///reports/datasets/MPNST_4/var_features_distribution.png)
- [Expression Violins (Top 6)](file:///reports/datasets/MPNST_4/top_features_violins.png)

---
## 5. Software & Environment Provenance
- **Seurat Version**: `5.4.0`
- **sctransform Version**: `0.4.3`
- **Git Commit Hash**: `31cc849c7c8aaf27e8f2b594fb40d0c54989a579`
