# Milestone 4 Normalization and Variable Feature Report - sample_3
*Generated on: 2026-07-09 18:35:55*

## 1. Parameters & Configuration

- **Input Filtered Object**: `results/datasets/sample_3/sample_3_filtered_specific.rds`
- **Output Normalized Object**: `results/datasets/sample_3/sample_3_normalized.rds`
- **Normalization Method**: `LogNormalize`
- **Requested Variable Features**: `1000`
- **Mitochondrial Regression**: `YES (percent.mt)`
- **Seed Used**: `42`
- **Execution Time**: `0.42 seconds`

## 2. Dataset Characteristics

- **Total Cells**: `150`
- **Total Raw Genes**: `120`

## 3. High Variance Features

Top 20 highly variable features selected under `LogNormalize` (sorted by variance.standardized):

1. **RPL-4** (variance.standardized: 1.2644, mean: 60.0733, variance: 82.0550)
2. **MT-7** (variance.standardized: 1.2445, mean: 36.4000, variance: 42.7651)
3. **Gene-14** (variance.standardized: 1.2269, mean: 12.1733, variance: 14.3590)
4. **MT-3** (variance.standardized: 1.1979, mean: 35.8667, variance: 42.9888)
5. **Gene-25** (variance.standardized: 1.1962, mean: 11.6600, variance: 13.5950)
6. **Gene-26** (variance.standardized: 1.1787, mean: 12.3267, variance: 14.7181)
7. **Gene-31** (variance.standardized: 1.1701, mean: 12.3200, variance: 14.5681)
8. **Gene-16** (variance.standardized: 1.1673, mean: 12.3867, variance: 14.9367)
9. **Gene-30** (variance.standardized: 1.1662, mean: 11.7800, variance: 13.0721)
10. **Gene-34** (variance.standardized: 1.1313, mean: 11.9667, variance: 12.5022)
11. **Gene-36** (variance.standardized: 1.1299, mean: 11.3667, variance: 12.5559)
12. **Gene-12** (variance.standardized: 1.1150, mean: 11.9333, variance: 12.3043)
13. **Gene-40** (variance.standardized: 1.1116, mean: 12.0733, variance: 12.5382)
14. **T-Marker-12** (variance.standardized: 1.1047, mean: 66.3733, variance: 6086.2489)
15. **RPS-1** (variance.standardized: 1.0982, mean: 59.3067, variance: 75.2073)
16. **Gene-38** (variance.standardized: 1.0924, mean: 11.9667, variance: 12.0727)
17. **Gene-8** (variance.standardized: 1.0918, mean: 11.5733, variance: 12.4208)
18. **RPS-3** (variance.standardized: 1.0896, mean: 60.4867, variance: 71.0434)
19. **Gene-23** (variance.standardized: 1.0835, mean: 12.7067, variance: 15.6315)
20. **RPL-1** (variance.standardized: 1.0754, mean: 60.0333, variance: 69.9251)

## 4. Diagnostics & Visualizations
The following diagnostic plots were generated to assess normalization quality:
- [Scatter Plot (HVF)](file:///reports/datasets/sample_3/var_features_scatter.png)
- [Distribution Plot (HVF Metric)](file:///reports/datasets/sample_3/var_features_distribution.png)
- [Expression Violins (Top 6)](file:///reports/datasets/sample_3/top_features_violins.png)

---
## 5. Software & Environment Provenance
- **Seurat Version**: `5.4.0`
- **sctransform Version**: `0.4.3`
- **Git Commit Hash**: `31cc849c7c8aaf27e8f2b594fb40d0c54989a579`
