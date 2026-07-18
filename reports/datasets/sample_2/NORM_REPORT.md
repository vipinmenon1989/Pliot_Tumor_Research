# Milestone 4 Normalization and Variable Feature Report - sample_2
*Generated on: 2026-07-09 18:35:01*

## 1. Parameters & Configuration

- **Input Filtered Object**: `results/datasets/sample_2/sample_2_filtered_specific.rds`
- **Output Normalized Object**: `results/datasets/sample_2/sample_2_normalized.rds`
- **Normalization Method**: `LogNormalize`
- **Requested Variable Features**: `1000`
- **Mitochondrial Regression**: `YES (percent.mt)`
- **Seed Used**: `42`
- **Execution Time**: `0.43 seconds`

## 2. Dataset Characteristics

- **Total Cells**: `150`
- **Total Raw Genes**: `120`

## 3. High Variance Features

Top 20 highly variable features selected under `LogNormalize` (sorted by variance.standardized):

1. **MT-6** (variance.standardized: 1.5784, mean: 161.2933, variance: 173.9402)
2. **Tumor-Marker-7** (variance.standardized: 1.4395, mean: 147.6267, variance: 33166.9335)
3. **MT-3** (variance.standardized: 1.4364, mean: 161.0867, variance: 173.0998)
4. **Tumor-Marker-19** (variance.standardized: 1.4167, mean: 147.6267, variance: 32641.7791)
5. **Tumor-Marker-13** (variance.standardized: 1.4036, mean: 147.5800, variance: 32820.0439)
6. **MT-5** (variance.standardized: 1.3604, mean: 161.6467, variance: 128.6327)
7. **Tumor-Marker-20** (variance.standardized: 1.3535, mean: 147.4467, variance: 33010.9602)
8. **Gene-8** (variance.standardized: 1.3333, mean: 20.2067, variance: 26.9436)
9. **Tumor-Marker-9** (variance.standardized: 1.2698, mean: 147.2533, variance: 32913.2911)
10. **Gene-38** (variance.standardized: 1.2601, mean: 20.2600, variance: 25.7776)
11. **Gene-37** (variance.standardized: 1.2498, mean: 19.8867, variance: 24.0877)
12. **Tumor-Marker-5** (variance.standardized: 1.2297, mean: 147.1600, variance: 32803.6521)
13. **Tumor-Marker-1** (variance.standardized: 1.2215, mean: 147.1800, variance: 32386.2828)
14. **Gene-34** (variance.standardized: 1.2133, mean: 20.1200, variance: 24.0795)
15. **Gene-12** (variance.standardized: 1.2061, mean: 19.6333, variance: 24.5157)
16. **MT-10** (variance.standardized: 1.1933, mean: 160.6733, variance: 171.8724)
17. **Gene-25** (variance.standardized: 1.1722, mean: 19.8267, variance: 22.7483)
18. **Gene-31** (variance.standardized: 1.1490, mean: 20.6133, variance: 23.9703)
19. **Gene-15** (variance.standardized: 1.1396, mean: 20.0600, variance: 22.3655)
20. **Gene-35** (variance.standardized: 1.1343, mean: 20.2400, variance: 23.0964)

## 4. Diagnostics & Visualizations
The following diagnostic plots were generated to assess normalization quality:
- [Scatter Plot (HVF)](file:///reports/datasets/sample_2/var_features_scatter.png)
- [Distribution Plot (HVF Metric)](file:///reports/datasets/sample_2/var_features_distribution.png)
- [Expression Violins (Top 6)](file:///reports/datasets/sample_2/top_features_violins.png)

---
## 5. Software & Environment Provenance
- **Seurat Version**: `5.4.0`
- **sctransform Version**: `0.4.3`
- **Git Commit Hash**: `31cc849c7c8aaf27e8f2b594fb40d0c54989a579`
