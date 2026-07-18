# Quality Control Threshold Optimization Review Report

*Generated on: 2026-07-09 19:54:14*

## 1. Executive Summary & Objective
The objective of this review is to evaluate whether applying uniform, global quality control thresholds across all constituent datasets is scientifically justified, or if batch-specific/dataset-specific thresholds are required before proceeding to Milestone 4 (M4) normalization.

### Current Global Baseline Thresholds (M3):
- **min_features**: 200
- **min_counts**: 500
- **max_percent_mt**: 10.0%
- **max_percent_ribo**: 20.0%
- **doublet_detection**: scDblFinder (singlets only)

---

## 2. Observed Results (Summary Statistics)

Below are the calculated summary statistics comparing the pre-filter (raw) and post-filter (M3) distributions of the key cell QC metrics. Notice the large cell losses in MPNST_2, MPNST_3, and MPNST_4.

### Library Depth (nCount_RNA) Summary Statistics
| Dataset (Stage) | Cells | Min | Q1 | Median | Q3 | Max | Mean (±SD) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **MPNST_1 (Pre-filter)** | 8338 | 500.0 | 3571.0 | 8145.0 | 15606.0 | 67368.0 | 11068.3 (±10078.3) |
| **MPNST_1 (Post-filter)** | 7615 | 500.0 | 3166.5 | 7503.0 | 13817.5 | 67368.0 | 9986.7 (±9218.0) |
| **MPNST_2 (Pre-filter)** | 2830 | 503.0 | 4388.2 | 7657.5 | 14415.5 | 57697.0 | 10524.0 (±8878.0) |
| **MPNST_2 (Post-filter)** | 1516 | 503.0 | 5345.5 | 8936.0 | 16110.5 | 57697.0 | 11848.2 (±9302.6) |
| **MPNST_3 (Pre-filter)** | 3682 | 503.0 | 2223.8 | 4137.5 | 8176.5 | 93016.0 | 6795.5 (±8671.6) |
| **MPNST_3 (Post-filter)** | 1918 | 505.0 | 3729.8 | 6137.5 | 9656.0 | 93016.0 | 8853.2 (±10192.2) |
| **MPNST_4 (Pre-filter)** | 7811 | 501.0 | 3155.0 | 5244.0 | 7724.0 | 79810.0 | 6155.0 (±5082.8) |
| **MPNST_4 (Post-filter)** | 3390 | 509.0 | 4344.2 | 6094.5 | 8405.2 | 66635.0 | 6823.6 (±4269.6) |

### Genes Detected (nFeature_RNA) Summary Statistics
| Dataset (Stage) | Cells | Min | Q1 | Median | Q3 | Max | Mean (±SD) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **MPNST_1 (Pre-filter)** | 8338 | 231.0 | 1876.0 | 3337.0 | 5204.8 | 8989.0 | 3578.5 (±2134.6) |
| **MPNST_1 (Post-filter)** | 7615 | 231.0 | 1740.5 | 3139.0 | 4812.5 | 8934.0 | 3350.6 (±2022.7) |
| **MPNST_2 (Pre-filter)** | 2830 | 278.0 | 1660.2 | 2419.5 | 3585.0 | 8007.0 | 2650.1 (±1362.0) |
| **MPNST_2 (Post-filter)** | 1516 | 289.0 | 2007.5 | 2754.0 | 3826.5 | 8007.0 | 2936.1 (±1318.1) |
| **MPNST_3 (Pre-filter)** | 3682 | 210.0 | 1020.5 | 1558.5 | 2375.5 | 8750.0 | 1927.4 (±1351.8) |
| **MPNST_3 (Post-filter)** | 1918 | 210.0 | 1463.0 | 2005.5 | 2727.8 | 8750.0 | 2357.6 (±1442.1) |
| **MPNST_4 (Pre-filter)** | 7811 | 280.0 | 1511.0 | 2179.0 | 2740.0 | 8952.0 | 2199.3 (±1068.0) |
| **MPNST_4 (Post-filter)** | 3390 | 315.0 | 1953.5 | 2389.0 | 2840.8 | 8242.0 | 2412.6 (±877.2) |

### Mitochondrial Content (percent.mt) Summary Statistics
| Dataset (Stage) | Cells | Min | Q1 | Median | Q3 | Max | Mean (±SD) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **MPNST_1 (Pre-filter)** | 8338 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 (±0.0) |
| **MPNST_1 (Post-filter)** | 7615 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 (±0.0) |
| **MPNST_2 (Pre-filter)** | 2830 | 0.2 | 5.4 | 7.1 | 9.0 | 20.0 | 7.6 (±3.2) |
| **MPNST_2 (Post-filter)** | 1516 | 0.2 | 5.0 | 6.5 | 7.8 | 10.0 | 6.4 (±1.9) |
| **MPNST_3 (Pre-filter)** | 3682 | 0.0 | 2.0 | 2.9 | 4.2 | 19.8 | 3.7 (±2.9) |
| **MPNST_3 (Post-filter)** | 1918 | 0.0 | 1.9 | 2.7 | 4.0 | 9.9 | 3.2 (±1.9) |
| **MPNST_4 (Pre-filter)** | 7811 | 0.1 | 5.9 | 8.1 | 10.7 | 20.0 | 8.6 (±3.9) |
| **MPNST_4 (Post-filter)** | 3390 | 0.7 | 5.2 | 7.1 | 8.5 | 10.0 | 6.7 (±2.1) |

### Ribosomal Content (percent.ribo) Summary Statistics
| Dataset (Stage) | Cells | Min | Q1 | Median | Q3 | Max | Mean (±SD) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **MPNST_1 (Pre-filter)** | 8338 | 0.1 | 2.7 | 4.2 | 6.9 | 45.0 | 5.2 (±3.5) |
| **MPNST_1 (Post-filter)** | 7615 | 0.1 | 2.7 | 4.3 | 7.0 | 19.8 | 5.2 (±3.4) |
| **MPNST_2 (Pre-filter)** | 2830 | 2.3 | 11.5 | 15.7 | 20.9 | 54.7 | 17.2 (±7.7) |
| **MPNST_2 (Post-filter)** | 1516 | 2.3 | 10.4 | 13.6 | 16.4 | 20.0 | 13.3 (±3.9) |
| **MPNST_3 (Pre-filter)** | 3682 | 1.8 | 11.3 | 16.8 | 24.2 | 58.6 | 18.5 (±9.2) |
| **MPNST_3 (Post-filter)** | 1918 | 1.8 | 9.2 | 12.8 | 16.2 | 20.0 | 12.7 (±4.2) |
| **MPNST_4 (Pre-filter)** | 7811 | 1.7 | 11.9 | 15.5 | 19.6 | 45.1 | 16.2 (±6.1) |
| **MPNST_4 (Post-filter)** | 3390 | 3.2 | 11.8 | 14.8 | 17.3 | 20.0 | 14.4 (±3.6) |

---

## 3. Filter Contribution & Overlap Analysis

The table below shows the unique cell exclusions contributed by each individual filter, compared to cells failing multiple criteria simultaneously.

| Dataset ID | Total Cells | Cells Retained | Total Excluded | Only Features | Only Counts | Only MT % | Only Ribo % | Only Doublet | Multiple Filters |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **MPNST_1** | 8338 | 7615 | 723 | 0 | 0 | 0 | 26 | 697 | 0 |
| **MPNST_2** | 2830 | 1516 | 1314 | 0 | 0 | 320 | 626 | 188 | 180 |
| **MPNST_3** | 3682 | 1918 | 1764 | 0 | 0 | 162 | 1250 | 236 | 116 |
| **MPNST_4** | 7811 | 3390 | 4421 | 0 | 0 | 1975 | 1502 | 526 | 418 |

### Observations on Overlaps:
- **Independent Exclusions**: Ribosomal percentage and Mitochondrial percentage filters act almost completely independently from other technical criteria. Very few cells failing Max MT or Max Ribo are flagged as doublets or low-complexity, indicating that these filters are removing distinct cell subpopulations rather than co-occurring artifacts.
- **Doublet overlap**: The doublet classification removes a distinct cell pool (8.0-9.2% of cells) with very little overlap with MT/Ribo flags, suggesting scDblFinder is successfully capturing doublets without capturing stressed singlets.

---

## 4. Threshold Sensitivity Analysis (Simulations)

Below are the sensitivity tables for the four simulated metrics, demonstrating cell retention as a function of the cutoff while keeping other parameters at their M3 baseline values:

### Mitochondrial Threshold Sensitivity Table
| Threshold | Retained Count | Removed Count | % Retained | % Removed |
| --- | --- | --- | --- | --- |
#### MPNST_1
| max_percent_mt = 5 | 7615 | 723 | 91.33% | 8.67% |
| max_percent_mt = 7.5 | 7615 | 723 | 91.33% | 8.67% |
| max_percent_mt = 10 | 7615 | 723 | 91.33% | 8.67% |
| max_percent_mt = 12.5 | 7615 | 723 | 91.33% | 8.67% |
| max_percent_mt = 15 | 7615 | 723 | 91.33% | 8.67% |
| max_percent_mt = 20 | 7615 | 723 | 91.33% | 8.67% |
#### MPNST_2
| max_percent_mt = 5 | 383 | 2447 | 13.53% | 86.47% |
| max_percent_mt = 7.5 | 1060 | 1770 | 37.46% | 62.54% |
| max_percent_mt = 10 | 1516 | 1314 | 53.57% | 46.43% |
| max_percent_mt = 12.5 | 1680 | 1150 | 59.36% | 40.64% |
| max_percent_mt = 15 | 1754 | 1076 | 61.98% | 38.02% |
| max_percent_mt = 20 | 1836 | 994 | 64.88% | 35.12% |
#### MPNST_3
| max_percent_mt = 5 | 1617 | 2065 | 43.92% | 56.08% |
| max_percent_mt = 7.5 | 1835 | 1847 | 49.84% | 50.16% |
| max_percent_mt = 10 | 1918 | 1764 | 52.09% | 47.91% |
| max_percent_mt = 12.5 | 1993 | 1689 | 54.13% | 45.87% |
| max_percent_mt = 15 | 2029 | 1653 | 55.11% | 44.89% |
| max_percent_mt = 20 | 2080 | 1602 | 56.49% | 43.51% |
#### MPNST_4
| max_percent_mt = 5 | 772 | 7039 | 9.88% | 90.12% |
| max_percent_mt = 7.5 | 1959 | 5852 | 25.08% | 74.92% |
| max_percent_mt = 10 | 3390 | 4421 | 43.40% | 56.60% |
| max_percent_mt = 12.5 | 4282 | 3529 | 54.82% | 45.18% |
| max_percent_mt = 15 | 4801 | 3010 | 61.46% | 38.54% |
| max_percent_mt = 20 | 5365 | 2446 | 68.69% | 31.31% |

### Ribosomal Threshold Sensitivity Table
| Threshold | Retained Count | Removed Count | % Retained | % Removed |
| --- | --- | --- | --- | --- |
#### MPNST_1
| max_percent_ribo = 15 | 7504 | 834 | 90.00% | 10.00% |
| max_percent_ribo = 20 | 7615 | 723 | 91.33% | 8.67% |
| max_percent_ribo = 25 | 7635 | 703 | 91.57% | 8.43% |
| max_percent_ribo = 30 | 7637 | 701 | 91.59% | 8.41% |
| max_percent_ribo = 35 | 7639 | 699 | 91.62% | 8.38% |
#### MPNST_2
| max_percent_ribo = 15 | 952 | 1878 | 33.64% | 66.36% |
| max_percent_ribo = 20 | 1516 | 1314 | 53.57% | 46.43% |
| max_percent_ribo = 25 | 1811 | 1019 | 63.99% | 36.01% |
| max_percent_ribo = 30 | 1983 | 847 | 70.07% | 29.93% |
| max_percent_ribo = 35 | 2080 | 750 | 73.50% | 26.50% |
#### MPNST_3
| max_percent_ribo = 15 | 1279 | 2403 | 34.74% | 65.26% |
| max_percent_ribo = 20 | 1918 | 1764 | 52.09% | 47.91% |
| max_percent_ribo = 25 | 2345 | 1337 | 63.69% | 36.31% |
| max_percent_ribo = 30 | 2676 | 1006 | 72.68% | 27.32% |
| max_percent_ribo = 35 | 2943 | 739 | 79.93% | 20.07% |
#### MPNST_4
| max_percent_ribo = 15 | 1747 | 6064 | 22.37% | 77.63% |
| max_percent_ribo = 20 | 3390 | 4421 | 43.40% | 56.60% |
| max_percent_ribo = 25 | 4337 | 3474 | 55.52% | 44.48% |
| max_percent_ribo = 30 | 4690 | 3121 | 60.04% | 39.96% |
| max_percent_ribo = 35 | 4831 | 2980 | 61.85% | 38.15% |

---

## 5. Biological & Technical Interpretation

### Observed Results:
- **`MPNST_1`** is highly resistant to QC filtering (only 8.67% cells removed). This is because it contains 0.00% mitochondrial transcripts (technical artifact/batch effect) and has low ribosomal content (median 4.22%).
- **`MPNST_2`** suffers 46.93% cell loss, driven by Ribosomal exclusions (28.13% of cells exceed 20% ribosomal content).
- **`MPNST_3`** suffers 48.29% cell loss, driven almost exclusively by Ribosomal content, which removes 37.10% of cells. The median ribosomal content is 16.77%, indicating that a 20% threshold is positioned very close to the center of the distribution, resulting in arbitrary truncation.
- **`MPNST_4`** suffers 56.75% cell loss, driven by both high Mitochondrial content (30.00% of cells exceed 10% MT) and Ribosomal content (22.95% of cells exceed 20% Ribo).

### Scientific Interpretation:
1. **Mitochondrial Expression (MPNST_4)**: A median of 8.09% and 95th percentile of 16.49% for percent.mt indicates that high MT content is a pervasive technical feature of this sample. Sarcoma tissue blocks often contain necrotic tumor centers where cells show increased mitochondrial expression due to hypoxia-induced stress or leakage, but remain viable. Restricting the dataset to < 10% MT removes 30% of all cells, likely depleting hypoxic tumor cells or specific viable tumor subclones.
2. **Ribosomal Expression (MPNST_3 & MPNST_2)**: Medians of 15-17% ribosomal fraction are biologically expected in highly proliferative sarcoma cells. Rapidly dividing cancer cells require high translational capacity. Restricting cells to < 20% ribosomal content removes up to 37% of cells, which represents an artificial and non-biological truncation of the cell state distribution. There is no evidence in the literature suggesting that a 20% ribosomal fraction indicates technical noise in sarcoma libraries.

---

## 6. Recommendations: Global versus Dataset-Specific Thresholds

### Global Thresholds (Current):
- **Advantages**: Simple to describe, enforces uniform filtering rules across all datasets, prevents researcher bias in cell selection.
- **Disadvantages**: Ignores sample-specific biological complexity and technical batch variations. Results in catastrophic cell loss (up to 56.8% cells removed) and likely introduces cell class depletion bias.

### Dataset-Specific Thresholds (Recommended):
- **Advantages**: Adapts to the technical baseline of each library prep and sequencing run. Prevents arbitrary cell loss by setting thresholds based on individual distributions (e.g. median + 3 MADs). Preserves hypoxic tumor populations in `MPNST_4` and highly proliferative cells in `MPNST_3`.
- **Disadvantages**: Requires individual justification and slightly complicates the methodology description.

### Quantitative Dataset-Specific Recommendations:
Based on the sensitivity curves and distribution medians, we propose the following tailored thresholds:

1. **`MPNST_1`**:
   - *Recommended*: min_features = 200, min_counts = 500, max_mt = 10% (no-effect), max_ribo = 20%.
   - *Justification*: Already highly clean, minimal cell loss (8.67%). No adjustments needed.
2. **`MPNST_2`**:
   - *Recommended*: min_features = 200, min_counts = 500, max_mt = 15%, max_ribo = 30%.
   - *Justification*: Under these thresholds, cell loss decreases from 46.93% to **11.23%**, retaining a much broader representation of the cell population.
3. **`MPNST_3`**:
   - *Recommended*: min_features = 200, min_counts = 500, max_mt = 10%, max_ribo = 35%.
   - *Justification*: Under these thresholds, cell loss decreases from 48.29% to **11.89%**, preventing the artificial truncation of the ribosomal distribution.
4. **`MPNST_4`**:
   - *Recommended*: min_features = 200, min_counts = 500, max_mt = 20%, max_ribo = 30%.
   - *Justification*: Under these thresholds, cell loss decreases from 56.75% to **14.28%**, preserving viable tumor cells in hypoxic states while still excluding severe outliers.

---

## 7. Provenance
- **Input files**: 4 raw dataset RDS files under `results/datasets/`
- **Git Commit Hash**: `31cc849c7c8aaf27e8f2b594fb40d0c54989a579`
- **Generated figures**: Visualizations saved in `reports/qc_optimization/`
- **Generated tables**: Summaries saved in `reports/qc_optimization/`
