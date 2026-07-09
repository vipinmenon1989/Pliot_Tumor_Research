# QC Recommendation Report - MPNST_3
*Generated on: 2026-07-08 19:30:05*

## 1. Observed Distributions Summary

The table below summarizes the observed distributions and key percentiles for the cell quality control metrics prior to filtering.

| Metric | Min | Mean | Median | Max | 1% | 5% | 10% | 90% | 95% | 99% |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **nCount_RNA** | 503 | 6795.5 | 4137 | 93016 | 580 | 981 | 1376 | 13155 | 21140 | 48469 |
| **nFeature_RNA** | 210 | 1927.4 | 1558 | 8750 | 383 | 564 | 734 | 3530 | 5053 | 6990 |
| **percent.mt** | 0.00% | 3.73% | 2.89% | 19.84% | 0.47% | 1.26% | 1.53% | 6.87% | 10.16% | 15.90% |
| **percent.ribo** | 1.75% | 18.48% | 16.77% | 58.60% | 4.24% | 6.61% | 7.91% | 32.17% | 36.06% | 41.69% |

## 2. Potential Outliers

Based on statistical distributions (e.g. median +/- 3 Median Absolute Deviations (MAD) or standard deviation thresholds):

- **Library Depth (nCount_RNA)**: The median value is 4137. Outliers are typically defined as cells with counts below -6236 or above 14511.
- **Genes Detected (nFeature_RNA)**: The median value is 1558. Outliers are typically cells with fewer than -1250 genes (likely empty droplets or low-quality cells) or more than 4367 genes (potential doublets).
- **Mitochondrial Fraction (percent.mt)**: The median is 2.89%, and the 95th percentile is 10.16%. Cells exceeding 10.16% show elevated mitochondrial expression indicating cellular stress or lysis.

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
| **Max Features** | `> 6000` | 94 | 2.55% |
| **Min Counts** | `< 500` | 0 | 0.00% |
| **Max Counts** | `> 50000` | 32 | 0.87% |
| **Max percent.mt** | `> 15.0%` | 55 | 1.49% |
| **Max percent.ribo** | `> 20.0%` | 1366 | 37.10% |
| **Combined Filters** | **All Above** | **1507** | **40.93%** |

### Summary of Expected Kept Cells:
- **Total cells before filtering**: 3682
- **Expected cells removed**: 1507 (40.93%)
- **Expected cells retained**: 2175 (59.07%)

## 5. Potential Biological & Technical Risks

1. **Mitochondrial Bias**: MPNST is a soft tissue sarcoma and tumor cells may exhibit metabolic shifts leading to higher baseline mitochondrial transcript percentages. Applying a strict 15% threshold could selectively remove tumor sub-clones or hypoxic cells. However, 15% provides a standard compromise for human cells.
2. **Ribosomal Variation**: Actively dividing tumor cells or specific cell classes (such as plasma cells) may have high ribosomal fractions. A 20% ribosomal threshold should be evaluated for potential cell-class depletion.
3. **Doublet Detection Confounding**: Cells with high counts are marked for exclusion by the max counts threshold, but a formal doublet detection (Milestone 3) will be needed to distinguish true high-transcript cells from doublets.

## 6. Questions Requiring Researcher Approval

Please review and approve the following settings prior to Milestone 3 (M3) filtering:

1. **Do you approve the uniform minimum gene limit of 200 features across this dataset?**
2. **Do you approve the mitochondrial limit of 15.0% for MPNST_3?** (Excludes 55 cells)
3. **Do you approve the ribosomal limit of 20.0% for MPNST_3?** (Excludes 1366 cells)
4. **Do you approve the overall threshold combination which will exclude 40.93% of cells?**

**NOTE: NO FILTERS HAVE BEEN APPLIED AT THIS STAGE. The dataset remains completely intact.**
