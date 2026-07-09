# QC Recommendation Report - sample_2
*Generated on: 2026-07-09 09:22:59*

## 1. Observed Distributions Summary

The table below summarizes the observed distributions and key percentiles for the cell quality control metrics prior to filtering.

| Metric | Min | Mean | Median | Max | 1% | 5% | 10% | 90% | 95% | 99% |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **nCount_RNA** | 9943 | 10865.8 | 10259 | 12424 | 9992 | 10077 | 10102 | 12262 | 12293 | 12368 |
| **nFeature_RNA** | 120 | 120.0 | 120 | 120 | 120 | 120 | 120 | 120 | 120 | 120 |
| **percent.mt** | 12.36% | 14.85% | 15.40% | 16.62% | 12.59% | 12.82% | 12.97% | 16.08% | 16.28% | 16.53% |
| **percent.ribo** | 7.76% | 9.25% | 9.59% | 10.48% | 7.80% | 8.01% | 8.09% | 10.10% | 10.17% | 10.41% |

## 2. Potential Outliers

Based on statistical distributions (e.g. median +/- 3 Median Absolute Deviations (MAD) or standard deviation thresholds):

- **Library Depth (nCount_RNA)**: The median value is 10259. Outliers are typically defined as cells with counts below 9656 or above 10862.
- **Genes Detected (nFeature_RNA)**: The median value is 120. Outliers are typically cells with fewer than 120 genes (likely empty droplets or low-quality cells) or more than 120 genes (potential doublets).
- **Mitochondrial Fraction (percent.mt)**: The median is 15.40%, and the 95th percentile is 16.28%. Cells exceeding 16.28% show elevated mitochondrial expression indicating cellular stress or lysis.

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
| **Min Features** | `< 200` | 150 | 100.00% |
| **Max Features** | `> 999999999` | 0 | 0.00% |
| **Min Counts** | `< 500` | 0 | 0.00% |
| **Max Counts** | `> 999999999` | 0 | 0.00% |
| **Max percent.mt** | `> 10.0%` | 150 | 100.00% |
| **Max percent.ribo** | `> 20.0%` | 0 | 0.00% |
| **Combined Filters** | **All Above** | **150** | **100.00%** |

### Summary of Expected Kept Cells:
- **Total cells before filtering**: 150
- **Expected cells removed**: 150 (100.00%)
- **Expected cells retained**: 0 (0.00%)

## 5. Potential Biological & Technical Risks

1. **Mitochondrial Bias**: MPNST is a soft tissue sarcoma and tumor cells may exhibit metabolic shifts leading to higher baseline mitochondrial transcript percentages. Applying a strict 15% threshold could selectively remove tumor sub-clones or hypoxic cells. However, 15% provides a standard compromise for human cells.
2. **Ribosomal Variation**: Actively dividing tumor cells or specific cell classes (such as plasma cells) may have high ribosomal fractions. A 20% ribosomal threshold should be evaluated for potential cell-class depletion.
3. **Doublet Detection Confounding**: Cells with high counts are marked for exclusion by the max counts threshold, but a formal doublet detection (Milestone 3) will be needed to distinguish true high-transcript cells from doublets.

## 6. Questions Requiring Researcher Approval

Please review and approve the following settings prior to Milestone 3 (M3) filtering:

1. **Do you approve the uniform minimum gene limit of 200 features across this dataset?**
2. **Do you approve the mitochondrial limit of 10.0% for sample_2?** (Excludes 150 cells)
3. **Do you approve the ribosomal limit of 20.0% for sample_2?** (Excludes 0 cells)
4. **Do you approve the overall threshold combination which will exclude 100.00% of cells?**

**NOTE: NO FILTERS HAVE BEEN APPLIED AT THIS STAGE. The dataset remains completely intact.**
