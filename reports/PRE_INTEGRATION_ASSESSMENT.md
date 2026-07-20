# Milestone 8 Combined Pre-Integration Assessment
*Generated on: 2026-07-19 23:25:04*

---

## 1. Executive Summary
This report provides a structural, non-integrated baseline assessment of the combined four MPNST Phase 1 datasets before any batch correction or integration is applied.

---

## 2. Dataset Composition
The cell count distribution across the combined datasets:

| Dataset ID | Cells | Percentage (%) |
| --- | :---: | :---: |
| **MPNST_1** | 7,615 | 38.62% |
| **MPNST_2** | 2,284 | 11.58% |
| **MPNST_3** | 2,940 | 14.91% |
| **MPNST_4** | 6,877 | 34.88% |
| **Total** | **19,716** | **100.00%** |

---

## 3. Observed Results

### 3.1 Combined Object Construction
The combined Seurat object was constructed by merging the four constituent dataset Seurat objects (`MPNST_1`, `MPNST_2`, `MPNST_3`, and `MPNST_4`). This resulted in a single object containing exactly **19,716 cells** and **31,764 genes** (default assay: `SCT`).

### 3.2 Shared Non-Integrated Dimensional Reductions
A single unified expression representation was computed by running `SCTransform` globally on the merged raw counts (regressing out `percent.mt` and selecting 3,000 variable features). We then computed:
1. **One shared non-integrated PCA** (retaining 30 dimensions).
2. **One shared neighbor graph** (using 30 dimensions).
3. **One shared pre-integration UMAP** (reduction name: `umap_preintegration`, seed: `42`).

*No batch correction, covariate alignment, or embedding concatenation was performed.*

### 3.3 Spatial Structure and Dataset Segregation
Visual inspection of the shared UMAP reveals complete spatial segregation. The cells from the four datasets form distinct, non-overlapping neighborhoods in the common embedding space, with minimal to zero spatial overlap.

### 3.4 Neighborhood Mixing Diagnostics
We calculated the fraction of 15 nearest neighbors in the shared PCA space that belong to the same dataset for each cell:

| Dataset | Mean Same-Dataset Neighbor Fraction | Mean Neighborhood Dataset Entropy |
| --- | :---: | :---: |
| **MPNST_1** | 0.9883 | 0.0328 |
| **MPNST_2** | 0.9102 | 0.2158 |
| **MPNST_3** | 0.9690 | 0.0751 |
| **MPNST_4** | 0.9779 | 0.0626 |
| **Global Mean** | **0.9613** | **0.0966** |

The global mean same-dataset neighbor fraction is **96.13%**, confirming that cells occupy neighborhoods composed almost exclusively of cells from their own constituent dataset.

### 3.5 QC Technical Covariates
Continuous UMAP plots show that the spatial separation is not driven by technical QC variables (`nCount_RNA`, `nFeature_RNA`). The distribution of library sizes and gene detection rates does not correlate with the major axes of UMAP separation, except for `percent.mt` in `MPNST_4`, where the high-stress cells form a distinct sub-neighborhood.

### 3.6 Confounding Analysis
Auditing of the metadata confirmed that clinical variables (such as anatomical site, tumor subtype, biological condition, primary vs. metastatic status, and patient demographics) are absent from the objects. Thus, the dataset/sample identities serve as the primary batch covariates. The technical batch effect and biological patient variance are 100% confounded.

---

## 4. Interpretation
The global segregation observed in the shared UMAP embedding suggests the presence of substantial dataset/sample-specific differences. Since each dataset corresponds to a separate clinical sample, these differences are likely a combined product of:
1. **Biological Heterogeneity**: Real tumor-to-tumor biological differences, donor-specific microenvironments, and cell-type composition variations.
2. **Potential Technical Contributors**: Unmeasured batch effects arising from separate library preparation chemistry, capture batches, sequencing runs, or operator handling.

---

## 5. Limitations
- **Confounding**: In the absence of clinical covariates or multi-sample patient designs, it is mathematically impossible to separate technical batch effects from true patient-specific biological differences.
- **Unmeasured Technical Variables**: Technical variables (sequencing run, library chemistry, operator, capture channel) are unmeasured, meaning we cannot definitively attribute segregation to technical batches.

---

## 6. Questions for Phase 2
1. **Integration Method Performance**: How do Harmony, Seurat CCA, and Seurat RPCA perform in aligning these four datasets?
2. **Biological Signal Preservation**: Does dataset integration collapse real biological boundaries, such as the `MPNST_4` stress-response cell state?
3. **Downstream Agreement**: Do integrated and non-integrated clustering assignments lead to consistent biological annotations and marker genes?
