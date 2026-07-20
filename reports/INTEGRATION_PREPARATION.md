# Milestone 8 Integration Preparation Report
*Generated on: 2026-07-19 23:25:04*

## 1. Integration Design Decisions

This document outlines the scientific preparation and design parameters for future batch correction and dataset integration, answering the key architectural questions based on the Milestone 8 pre-integration baseline.

---

### Question 1: Is strong dataset/sample-associated structure present?
**Answer**: Yes. Both visual inspection of the shared UMAP and quantitative nearest-neighbor diagnostics (mean same-dataset neighbor fraction of **96.13%**) confirm complete spatial segregation of the four constituent datasets (`MPNST_1`, `MPNST_2`, `MPNST_3`, and `MPNST_4`).

### Question 2: Is separation global or population-specific?
**Answer**: The separation is global. The spatial segregation is observed across the entire UMAP embedding, affecting all major transcriptomic zones and cell classes (rather than being restricted to a single cell lineage or specific cluster).

### Question 3: Which technical variables are actually observed?
**Answer**: None. There are no technical covariates (such as library preparation date, sequencing run, operator, flow cell ID, or capture batch) explicitly annotated in the metadata of the constituent Seurat objects. The variables `sample_id` and `orig.ident` act as the sole proxies for the combined batch/biological variables.

### Question 4: Which technical factors are only plausible/unmeasured?
**Answer**: Library preparation chemistry, cell capture batch (e.g., Chromium chip channel), sequencing chemistry/batch, operator effects, and instrument variations are unmeasured but highly plausible technical contributors to the observed dataset segregation.

### Question 5: Which biological variables are observed?
**Answer**: The raw transcriptomic counts (RNA assay) and the legacy cell type annotations (`orig.anno`) are the only biological variables observed. The legacy annotation contains 16 cell classes (e.g., Malignant NC-MES, Immature SC-like, Macrophages).

### Question 6: What confounding is demonstrable?
**Answer**: The dataset identity (`sample_id`) is 100% confounded with the individual patient/sample identity. Because each dataset represents exactly one clinical sample, technical batch variance and patient-specific biological variance are completely inseparable.

### Question 7: What confounding remains unknown?
**Answer**: Any confounding associated with patient demographics, tumor site (primary vs. metastatic), anatomical location (e.g., sciatic nerve vs. retroperitoneum), clinical treatment history (pre- vs. post-chemotherapy), or tumor subtype remains entirely unknown due to the lack of clinical metadata.

### Question 8: What biological signal could integration remove?
**Answer**: Aggressive integration could over-smooth and erase:
1. **Patient-specific tumor subclones**: Unique malignant programs or copy-number variant (CNV) driven states.
2. **Cellular composition differences**: Real differences in immune infiltration (e.g., macrophage polarization) or stromal density.
3. **Malignant functional states**: The tumor-specific stress-response and hypoxia-like programs observed in `MPNST_4`.

### Question 9: Should the non-integrated M8 baseline remain a permanent reference?
**Answer**: Yes. The non-integrated shared PCA and UMAP coordinates must remain a permanent reference. Any post-integration coordinates generated in Phase 2 must be evaluated against this baseline to ensure that biological boundaries are not collapsed and that artificial cell states are not introduced.

### Question 10: Should integration methods be evaluated in Phase 2?
**Answer**: Yes. Strong dataset-associated structure is present before integration. To perform unified multi-sample downstream tasks (such as joint cell-type annotation, comparative differential expression, and pathway analysis), integration should be evaluated in Phase 2 against the frozen non-integrated baseline.

### Question 11: Which methods may reasonably be benchmarked later?
**Answer**: We recommend benchmarking:
1. **Harmony**: A soft k-means clustering-based integration method that is fast and suitable for large datasets.
2. **Seurat CCA (Canonical Correlation Analysis)**: Suitable for identifying shared correlation structures when cell types are shared.
3. **Seurat RPCA (Reciprocal PCA)**: A faster, less aggressive alternative to CCA that projects datasets into each other's PCA space, useful when cell type composition varies significantly.

### Question 12: What metrics should assess integration quality?
**Answer**:
1. **Dataset Mixing**:
   - **iLISI (integration Local Inverse Simpson Index)**: Measures the local mixing of datasets in the neighborhood of each cell (target: high value).
   - **kBET (k-nearest neighbor Batch Effect Test)**: Statistically checks if the local neighborhood distribution is representative of the global dataset composition.
2. **Biological Conservation**:
   - **cLISI (cell-type Local Inverse Simpson Index)**: Evaluates whether biological cell-type boundaries are preserved (target: low value, indicating cell types do not mix artificially).
   - **Marker Retention**: Verifies that the top marker genes identified in M7 remain highly specific to their respective cell types post-integration.
