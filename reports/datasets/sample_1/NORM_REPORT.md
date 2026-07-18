# Milestone 4 Normalization and Variable Feature Report - sample_1
*Generated on: 2026-07-09 18:36:23*

## 1. Parameters & Configuration

- **Input Filtered Object**: `results/datasets/sample_1/sample_1_filtered_specific.rds`
- **Output Normalized Object**: `results/datasets/sample_1/sample_1_normalized.rds`
- **Normalization Method**: `LogNormalize`
- **Requested Variable Features**: `1000`
- **Mitochondrial Regression**: `YES (percent.mt)`
- **Seed Used**: `42`
- **Execution Time**: `0.55 seconds`

## 2. Dataset Characteristics

- **Total Cells**: `150`
- **Total Raw Genes**: `120`

## 3. High Variance Features

Top 20 highly variable features selected under `LogNormalize` (sorted by variance.standardized):

1. **RPL-5** (variance.standardized: 1.2685, mean: 74.5533, variance: 95.7253)
2. **Gene-20** (variance.standardized: 1.2160, mean: 14.9267, variance: 18.1758)
3. **RPL-4** (variance.standardized: 1.2031, mean: 75.2133, variance: 90.4374)
4. **Gene-7** (variance.standardized: 1.1807, mean: 14.6400, variance: 17.7219)
5. **Gene-23** (variance.standardized: 1.1712, mean: 14.9000, variance: 17.5134)
6. **Gene-18** (variance.standardized: 1.1621, mean: 14.6667, variance: 17.4855)
7. **Gene-2** (variance.standardized: 1.1560, mean: 14.7400, variance: 17.4420)
8. **MT-2** (variance.standardized: 1.1549, mean: 44.9400, variance: 54.1105)
9. **Gene-6** (variance.standardized: 1.1523, mean: 14.9067, variance: 17.2261)
10. **RPS-1** (variance.standardized: 1.1426, mean: 75.2800, variance: 86.1358)
11. **Gene-34** (variance.standardized: 1.1305, mean: 14.1600, variance: 14.9407)
12. **MT-3** (variance.standardized: 1.1291, mean: 44.6000, variance: 53.2886)
13. **Gene-8** (variance.standardized: 1.1241, mean: 14.9667, variance: 16.9183)
14. **Gene-40** (variance.standardized: 1.1225, mean: 15.3933, variance: 17.6094)
15. **MT-8** (variance.standardized: 1.1179, mean: 45.3200, variance: 51.9506)
16. **Gene-24** (variance.standardized: 1.0883, mean: 15.3533, variance: 17.0220)
17. **Gene-1** (variance.standardized: 1.0759, mean: 15.1533, variance: 16.6139)
18. **Gene-38** (variance.standardized: 1.0754, mean: 15.2267, variance: 16.6865)
19. **RPS-3** (variance.standardized: 1.0633, mean: 75.0133, variance: 79.5300)
20. **Myeloid-Marker-8** (variance.standardized: 1.0627, mean: 83.3667, variance: 9578.7975)

## 4. Diagnostics & Visualizations
The following diagnostic plots were generated to assess normalization quality:
- [Scatter Plot (HVF)](file:///reports/datasets/sample_1/var_features_scatter.png)
- [Distribution Plot (HVF Metric)](file:///reports/datasets/sample_1/var_features_distribution.png)
- [Expression Violins (Top 6)](file:///reports/datasets/sample_1/top_features_violins.png)

---
## 5. Software & Environment Provenance
- **Seurat Version**: `5.4.0`
- **sctransform Version**: `0.4.3`
- **Git Commit Hash**: `31cc849c7c8aaf27e8f2b594fb40d0c54989a579`
