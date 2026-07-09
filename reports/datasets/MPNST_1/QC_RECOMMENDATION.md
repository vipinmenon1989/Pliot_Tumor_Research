# QC Recommendation Report - MPNST_1
*Generated on: 2026-07-09 09:26:25*

## 1. Observed Distributions Summary

The table below summarizes the observed distributions and key percentiles for the cell quality control metrics prior to filtering.

| Metric | Min | Mean | Median | Max | 1% | 5% | 10% | 90% | 95% | 99% |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **nCount_RNA** | 500 | 11068.3 | 8145 | 67368 | 532 | 721 | 1085 | 24962 | 31762 | 44875 |
| **nFeature_RNA** | 231 | 3578.5 | 3337 | 8989 | 377 | 502 | 701 | 6635 | 7304 | 8388 |
| **percent.mt** | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| **percent.ribo** | 0.08% | 5.17% | 4.22% | 45.00% | 0.46% | 1.13% | 1.69% | 10.01% | 12.04% | 16.10% |

## 2. Potential Outliers

Based on statistical distributions (e.g. median +/- 3 Median Absolute Deviations (MAD) or standard deviation thresholds):

- **Library Depth (nCount_RNA)**: The median value is 8145. Outliers are typically defined as cells with counts below -16240 or above 32530.
- **Genes Detected (nFeature_RNA)**: The median value is 3337. Outliers are typically cells with fewer than -3992 genes (likely empty droplets or low-quality cells) or more than 10666 genes (potential doublets).
- **Mitochondrial Fraction (percent.mt)**: The median is 0.00%, and the 95th percentile is 0.00%. Cells exceeding 0.00% show elevated mitochondrial expression indicating cellular stress or lysis.

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
| **Max percent.mt** | `> 10.0%` | 0 | 0.00% |
| **Max percent.ribo** | `> 20.0%` | 26 | 0.31% |
| **Combined Filters** | **All Above** | **26** | **0.31%** |

### Summary of Expected Kept Cells:
- **Total cells before filtering**: 8338
- **Expected cells removed**: 26 (0.31%)
- **Expected cells retained**: 8312 (99.69%)

## 5. Potential Biological & Technical Risks

1. **Mitochondrial Bias**: MPNST is a soft tissue sarcoma and tumor cells may exhibit metabolic shifts leading to higher baseline mitochondrial transcript percentages. Applying a strict 15% threshold could selectively remove tumor sub-clones or hypoxic cells. However, 15% provides a standard compromise for human cells.
2. **Ribosomal Variation**: Actively dividing tumor cells or specific cell classes (such as plasma cells) may have high ribosomal fractions. A 20% ribosomal threshold should be evaluated for potential cell-class depletion.
3. **Doublet Detection Confounding**: Cells with high counts are marked for exclusion by the max counts threshold, but a formal doublet detection (Milestone 3) will be needed to distinguish true high-transcript cells from doublets.

## 6. Questions Requiring Researcher Approval

Please review and approve the following settings prior to Milestone 3 (M3) filtering:

1. **Do you approve the uniform minimum gene limit of 200 features across this dataset?**
2. **Do you approve the mitochondrial limit of 10.0% for MPNST_1?** (Excludes 0 cells)
3. **Do you approve the ribosomal limit of 20.0% for MPNST_1?** (Excludes 26 cells)
4. **Do you approve the overall threshold combination which will exclude 0.31% of cells?**

**NOTE: NO FILTERS HAVE BEEN APPLIED AT THIS STAGE. The dataset remains completely intact.**
