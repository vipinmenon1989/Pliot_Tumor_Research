# QC Recommendation Report - sample_1
*Generated on: 2026-07-09 09:23:46*

## 1. Observed Distributions Summary

The table below summarizes the observed distributions and key percentiles for the cell quality control metrics prior to filtering.

| Metric | Min | Mean | Median | Max | 1% | 5% | 10% | 90% | 95% | 99% |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **nCount_RNA** | 6632 | 7388.8 | 6937 | 8560 | 6705 | 6789 | 6827 | 8414 | 8472 | 8536 |
| **nFeature_RNA** | 120 | 120.0 | 120 | 120 | 120 | 120 | 120 | 120 | 120 | 120 |
| **percent.mt** | 4.79% | 6.14% | 6.37% | 7.14% | 4.93% | 5.15% | 5.30% | 6.81% | 6.88% | 7.03% |
| **percent.ribo** | 8.24% | 10.25% | 10.56% | 12.12% | 8.31% | 8.55% | 8.79% | 11.27% | 11.41% | 11.78% |

## 2. Potential Outliers

Based on statistical distributions (e.g. median +/- 3 Median Absolute Deviations (MAD) or standard deviation thresholds):

- **Library Depth (nCount_RNA)**: The median value is 6937. Outliers are typically defined as cells with counts below 6501 or above 7372.
- **Genes Detected (nFeature_RNA)**: The median value is 120. Outliers are typically cells with fewer than 120 genes (likely empty droplets or low-quality cells) or more than 120 genes (potential doublets).
- **Mitochondrial Fraction (percent.mt)**: The median is 6.37%, and the 95th percentile is 6.88%. Cells exceeding 6.88% show elevated mitochondrial expression indicating cellular stress or lysis.

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
| **Max percent.mt** | `> 10.0%` | 0 | 0.00% |
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
2. **Do you approve the mitochondrial limit of 10.0% for sample_1?** (Excludes 0 cells)
3. **Do you approve the ribosomal limit of 20.0% for sample_1?** (Excludes 0 cells)
4. **Do you approve the overall threshold combination which will exclude 100.00% of cells?**

**NOTE: NO FILTERS HAVE BEEN APPLIED AT THIS STAGE. The dataset remains completely intact.**
