# Milestone 4 Normalization and Variable Feature Report - sample_4
*Generated on: 2026-07-09 18:35:01*

## 1. Parameters & Configuration

- **Input Filtered Object**: `results/datasets/sample_4/sample_4_filtered_specific.rds`
- **Output Normalized Object**: `results/datasets/sample_4/sample_4_normalized.rds`
- **Normalization Method**: `LogNormalize`
- **Requested Variable Features**: `1000`
- **Mitochondrial Regression**: `YES (percent.mt)`
- **Seed Used**: `42`
- **Execution Time**: `0.52 seconds`

## 2. Dataset Characteristics

- **Total Cells**: `150`
- **Total Raw Genes**: `120`

## 3. High Variance Features

Top 20 highly variable features selected under `LogNormalize` (sorted by variance.standardized):

1. **RPL-2** (variance.standardized: 1.3197, mean: 91.1800, variance: 121.6654)
2. **Gene-14** (variance.standardized: 1.1988, mean: 18.1200, variance: 22.3613)
3. **Gene-2** (variance.standardized: 1.1781, mean: 17.5467, variance: 21.2428)
4. **Gene-39** (variance.standardized: 1.1752, mean: 18.0133, variance: 22.1206)
5. **RPL-1** (variance.standardized: 1.1721, mean: 90.1733, variance: 101.0167)
6. **Gene-33** (variance.standardized: 1.1713, mean: 17.8267, variance: 21.7281)
7. **Gene-13** (variance.standardized: 1.1583, mean: 18.4133, variance: 21.0495)
8. **Gene-5** (variance.standardized: 1.1433, mean: 17.9000, variance: 21.3658)
9. **MT-10** (variance.standardized: 1.1397, mean: 55.3600, variance: 70.2319)
10. **Gene-19** (variance.standardized: 1.1174, mean: 17.9667, variance: 21.0391)
11. **Gene-36** (variance.standardized: 1.1067, mean: 19.0000, variance: 20.5369)
12. **MT-3** (variance.standardized: 1.0988, mean: 54.1933, variance: 62.9892)
13. **Gene-18** (variance.standardized: 1.0961, mean: 17.3000, variance: 19.5403)
14. **T-Marker-9** (variance.standardized: 1.0940, mean: 100.2733, variance: 13808.4416)
15. **Gene-15** (variance.standardized: 1.0792, mean: 18.3467, variance: 19.5434)
16. **Gene-38** (variance.standardized: 1.0720, mean: 17.7200, variance: 19.5989)
17. **RPS-2** (variance.standardized: 1.0710, mean: 88.7800, variance: 90.0788)
18. **MT-2** (variance.standardized: 1.0649, mean: 53.3267, variance: 57.7114)
19. **Gene-12** (variance.standardized: 1.0599, mean: 18.4000, variance: 19.2483)
20. **Gene-37** (variance.standardized: 1.0526, mean: 17.6867, variance: 19.1696)

## 4. Diagnostics & Visualizations
The following diagnostic plots were generated to assess normalization quality:
- [Scatter Plot (HVF)](file:///reports/datasets/sample_4/var_features_scatter.png)
- [Distribution Plot (HVF Metric)](file:///reports/datasets/sample_4/var_features_distribution.png)
- [Expression Violins (Top 6)](file:///reports/datasets/sample_4/top_features_violins.png)

---
## 5. Software & Environment Provenance
- **Seurat Version**: `5.4.0`
- **sctransform Version**: `0.4.3`
- **Git Commit Hash**: `31cc849c7c8aaf27e8f2b594fb40d0c54989a579`
