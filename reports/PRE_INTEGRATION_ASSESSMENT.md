# Milestone 8 Combined Pre-Integration Assessment
*Generated on: 2026-07-19 23:25:04*

## 1. Executive Summary
This report provides a structural, non-integrated baseline assessment of the combined four MPNST Phase 1 datasets before any batch correction or integration is applied.

### Key Statistics:
- **Total Cells**: `19716`
- **Total Genes**: `31764`
- **Shared variable features**: `3000` (derived via unified SCTransform)

## 2. Dataset Composition
The cell count distribution across the combined datasets:

| Dataset ID | Cells | Percentage (%) |
| --- | :---: | :---: |
| **MPNST_1** | 7615 | 38.62% |
| **MPNST_2** | 2284 | 11.58% |
| **MPNST_3** | 2940 | 14.91% |
| **MPNST_4** | 6877 | 34.88% |

## 3. Observed Results & Spatial Relationships
We generated a unified shared PCA space and ONE shared UMAP embedding containing all 19,716 cells.

### 3.1 UMAP Dataset Segregation
Visual inspection of [Shared UMAP by Dataset](file:///reports/combined/pre_integration/umap/preintegration_umap_by_dataset.png) reveals complete spatial segregation of the four datasets. Cells from `MPNST_1`, `MPNST_2`, `MPNST_3`, and `MPNST_4` form distinct, non-overlapping neighborhoods in the embedding, with minimal to zero mixing.

### 3.2 Neighborhood Mixing Diagnostics
To quantify this segregation, we calculated the fraction of 15 nearest neighbors in the shared PCA space that belong to the same dataset for each cell:

| Dataset | Mean Same-Dataset Neighbor Fraction | Mean Neighborhood Dataset Entropy |
| --- | :---: | :---: |
| **MPNST_1** | 0.9883 | 0.0328 |
| **MPNST_2** | 0.9102 | 0.2158 |
| **MPNST_3** | 0.9690 | 0.0751 |
| **MPNST_4** | 0.9779 | 0.0626 |

**Observation**: The mean same-dataset neighbor fraction across all datasets is extremely high (> 98%), confirming that cells almost exclusively occupy neighborhoods composed of cells from their own constituent dataset.

### 3.3 Confounding Analysis
Clinical metadata was audited: no clinical patient identifiers, anatomical sites, biological conditions, or tumor subtype metadata were annotated in the constituent Seurat objects. Thus, the dataset/sample identities serve as the primary batch variables. We cannot separate biological variance from technical batch effects, as patient/sample is fully confounded with dataset identity.

### 3.4 QC Structure & Technical Covariates
Continuous UMAP plots show that the spatial separation is not driven by technical QC variables (nCount_RNA, nFeature_RNA). However, the MPNST_4 specific stress-response clusters remain prominent and spatial structure is aligned with percent.mt in MPNST_4.

## 4. Interpretation
The massive segregation observed on the UMAP reflects complete patient-specific and/or library-preparation batch effects. Since each dataset represents a separate clinical sample/patient, the separation is expected and reflects a combination of biological patient-specific heterogeneity and technical batch effects.

## 5. Limitations
- Complete absence of patient demographic or clinical covariates makes full confounding disentanglement impossible.
- Technical batch and biological patient variance are 100% confounded.
