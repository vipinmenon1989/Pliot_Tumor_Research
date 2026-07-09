# QC Recommendation Report - MPNST_4
*Generated on: 2026-07-08 19:31:27*

## 1. Observed Distributions Summary

The table below summarizes the observed distributions and key percentiles for the cell quality control metrics prior to filtering.

| Metric | Min | Mean | Median | Max | 1% | 5% | 10% | 90% | 95% | 99% |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **nCount_RNA** | 501 | 6155.0 | 5244 | 79810 | 555 | 862 | 1462 | 11216 | 14130 | 24161 |
| **nFeature_RNA** | 280 | 2199.3 | 2179 | 8952 | 383 | 530 | 803 | 3443 | 4008 | 5480 |
| **percent.mt** | 0.07% | 8.59% | 8.09% | 19.98% | 1.49% | 3.00% | 3.88% | 14.23% | 16.49% | 19.11% |
| **percent.ribo** | 1.69% | 16.17% | 15.52% | 45.10% | 5.36% | 7.55% | 9.09% | 24.03% | 27.18% | 34.10% |

## 2. Potential Outliers

Based on statistical distributions (e.g. median +/- 3 Median Absolute Deviations (MAD) or standard deviation thresholds):

- **Library Depth (nCount_RNA)**: The median value is 5244. Outliers are typically defined as cells with counts below -4745 or above 15233.
- **Genes Detected (nFeature_RNA)**: The median value is 2179. Outliers are typically cells with fewer than -551 genes (likely empty droplets or low-quality cells) or more than 4909 genes (potential doublets).
- **Mitochondrial Fraction (percent.mt)**: The median is 8.09%, and the 95th percentile is 16.49%. Cells exceeding 16.49% show elevated mitochondrial expression indicating cellular stress or lysis.

## 3. Recommended Thresholds and Rationale

Below are the recommended thresholds driven by the project configuration:

| QC Metric | Recommended Threshold | Rationale |
| --- | --- | --- |
| **Min Features (nFeature_RNA)** | `> 200` | Exclude low-complexity droplets/dead cells that do not contain sufficient biological signal. |
| **Max Features (nFeature_RNA)** | `< 6000` | Exclude potential doublets or multi-cell aggregates. |
| **Min Counts (nCount_RNA)** | `> 500` | Ensure sufficient library depth for robust gene expression estimation. |
| **Max Counts (nCount_RNA)** | `< 50000` | Exclude cells with abnormally high UMI counts, indicating technical artifacts or doublets. |
| **Max percent.mt** | `< 15.0%` | Standard filter to remove dying or damaged cells which release cytoplasmic RNA and retain mitochondrial transcripts. |
| **Max percent.ribo** | `< 20.0%` | Eliminate cells with extremely high ribosomal expression, which may represent technical bias or specific translation stress. |

## 4. Expected Filtering Impact

Here is the estimated impact of applying each of the recommended thresholds independently and in combination:

| Filter Metric | Threshold | Cells Excluded | % Excluded |
| --- | --- | --- | --- |
| **Min Features** | `< 200` | 0 | 0.00% |
| **Max Features** | `> 6000` | 49 | 0.63% |
| **Min Counts** | `< 500` | 0 | 0.00% |
| **Max Counts** | `> 50000` | 11 | 0.14% |
| **Max percent.mt** | `> 15.0%` | 621 | 7.95% |
| **Max percent.ribo** | `> 20.0%` | 1793 | 22.95% |
| **Combined Filters** | **All Above** | **2403** | **30.76%** |

### Summary of Expected Kept Cells:
- **Total cells before filtering**: 7811
- **Expected cells removed**: 2403 (30.76%)
- **Expected cells retained**: 5408 (69.24%)

## 5. Potential Biological & Technical Risks

1. **Mitochondrial Bias**: MPNST is a soft tissue sarcoma and tumor cells may exhibit metabolic shifts leading to higher baseline mitochondrial transcript percentages. Applying a strict 15% threshold could selectively remove tumor sub-clones or hypoxic cells. However, 15% provides a standard compromise for human cells.
2. **Ribosomal Variation**: Actively dividing tumor cells or specific cell classes (such as plasma cells) may have high ribosomal fractions. A 20% ribosomal threshold should be evaluated for potential cell-class depletion.
3. **Doublet Detection Confounding**: Cells with high counts are marked for exclusion by the max counts threshold, but a formal doublet detection (Milestone 3) will be needed to distinguish true high-transcript cells from doublets.

## 6. Questions Requiring Researcher Approval

Please review and approve the following settings prior to Milestone 3 (M3) filtering:

1. **Do you approve the uniform minimum gene limit of 200 features across this dataset?**
2. **Do you approve the mitochondrial limit of 15.0% for MPNST_4?** (Excludes 621 cells)
3. **Do you approve the ribosomal limit of 20.0% for MPNST_4?** (Excludes 1793 cells)
4. **Do you approve the overall threshold combination which will exclude 30.76% of cells?**

**NOTE: NO FILTERS HAVE BEEN APPLIED AT THIS STAGE. The dataset remains completely intact.**
