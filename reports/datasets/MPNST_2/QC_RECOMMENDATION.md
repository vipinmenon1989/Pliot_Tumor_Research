# QC Recommendation Report - MPNST_2
*Generated on: 2026-07-09 09:28:28*

## 1. Observed Distributions Summary

The table below summarizes the observed distributions and key percentiles for the cell quality control metrics prior to filtering.

| Metric | Min | Mean | Median | Max | 1% | 5% | 10% | 90% | 95% | 99% |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **nCount_RNA** | 503 | 10524.0 | 7657 | 57697 | 662 | 1398 | 2209 | 22697 | 28581 | 40392 |
| **nFeature_RNA** | 278 | 2650.1 | 2419 | 8007 | 435 | 736 | 963 | 4541 | 5126 | 6212 |
| **percent.mt** | 0.19% | 7.58% | 7.10% | 19.99% | 1.68% | 3.27% | 4.10% | 11.83% | 14.08% | 18.11% |
| **percent.ribo** | 2.35% | 17.17% | 15.67% | 54.71% | 4.97% | 7.27% | 8.80% | 28.42% | 32.50% | 40.47% |

## 2. Potential Outliers

Based on statistical distributions (e.g. median +/- 3 Median Absolute Deviations (MAD) or standard deviation thresholds):

- **Library Depth (nCount_RNA)**: The median value is 7657. Outliers are typically defined as cells with counts below -11276 or above 26591.
- **Genes Detected (nFeature_RNA)**: The median value is 2419. Outliers are typically cells with fewer than -1823 genes (likely empty droplets or low-quality cells) or more than 6662 genes (potential doublets).
- **Mitochondrial Fraction (percent.mt)**: The median is 7.10%, and the 95th percentile is 14.08%. Cells exceeding 14.08% show elevated mitochondrial expression indicating cellular stress or lysis.

## 3. Recommended Thresholds and Rationale

Below are the recommended thresholds driven by the project configuration:

| QC Metric | Recommended Threshold | Rationale |
| --- | --- | --- |
| **Min Features (nFeature_RNA)** | `> 200` | Exclude low-complexity droplets/dead cells that do not contain sufficient biological signal. |
| **Max Features (nFeature_RNA)** | `< 999999999` | Exclude potential doublets or multi-cell aggregates. |
| **Min Counts (nCount_RNA)** | `> 500` | Ensure sufficient library depth for robust gene expression estimation. |
| **Max Counts (nCount_RNA)** | `< 999999999` | Exclude cells with abnormally high UMI counts, indicating technical artifacts or doublets. |
| **Max percent.mt** | `< 10.0%` | Standard filter to remove dying or damaged cells which release cytoplasmic RNA and retain mitochondrial transcripts. |
| **Max percent.ribo** | `< 20.0%` | Eliminate cells with extremely high ribosomal expression, which may represent technical bias or specific translation stress. |

## 4. Expected Filtering Impact

Here is the estimated impact of applying each of the recommended thresholds independently and in combination:

| Filter Metric | Threshold | Cells Excluded | % Excluded |
| --- | --- | --- | --- |
| **Min Features** | `< 200` | 0 | 0.00% |
| **Max Features** | `> 999999999` | 0 | 0.00% |
| **Min Counts** | `< 500` | 0 | 0.00% |
| **Max Counts** | `> 999999999` | 0 | 0.00% |
| **Max percent.mt** | `> 10.0%` | 487 | 17.21% |
| **Max percent.ribo** | `> 20.0%` | 796 | 28.13% |
| **Combined Filters** | **All Above** | **1126** | **39.79%** |

### Summary of Expected Kept Cells:
- **Total cells before filtering**: 2830
- **Expected cells removed**: 1126 (39.79%)
- **Expected cells retained**: 1704 (60.21%)

## 5. Potential Biological & Technical Risks

1. **Mitochondrial Bias**: MPNST is a soft tissue sarcoma and tumor cells may exhibit metabolic shifts leading to higher baseline mitochondrial transcript percentages. Applying a strict 15% threshold could selectively remove tumor sub-clones or hypoxic cells. However, 15% provides a standard compromise for human cells.
2. **Ribosomal Variation**: Actively dividing tumor cells or specific cell classes (such as plasma cells) may have high ribosomal fractions. A 20% ribosomal threshold should be evaluated for potential cell-class depletion.
3. **Doublet Detection Confounding**: Cells with high counts are marked for exclusion by the max counts threshold, but a formal doublet detection (Milestone 3) will be needed to distinguish true high-transcript cells from doublets.

## 6. Questions Requiring Researcher Approval

Please review and approve the following settings prior to Milestone 3 (M3) filtering:

1. **Do you approve the uniform minimum gene limit of 200 features across this dataset?**
2. **Do you approve the mitochondrial limit of 10.0% for MPNST_2?** (Excludes 487 cells)
3. **Do you approve the ribosomal limit of 20.0% for MPNST_2?** (Excludes 796 cells)
4. **Do you approve the overall threshold combination which will exclude 39.79% of cells?**

**NOTE: NO FILTERS HAVE BEEN APPLIED AT THIS STAGE. The dataset remains completely intact.**
