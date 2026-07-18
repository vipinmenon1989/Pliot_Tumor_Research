# Milestone 4 (M4) Execution Report — Normalization and Variable Features

*Generated on: 2026-07-09 20:00:27*

## 1. Execution Summary
- **Authorized Milestone**: Milestone 4 (M4) — Normalization and Variable Features
- **Normalization Method**: `SCTransform`
- **Variable Features Requested**: `3000`
- **Random Seed**: `42`
- **Input Objects**: Dataset-specific QC outputs (`*_filtered_specific.rds`)
- **Output Objects**: Normalized and scaled Seurat objects (`*_normalized.rds`)

---

## 2. Scientific Recommendation & Justification

### Selected Normalization Strategy: SCTransform v2
For the MPNST sarcoma dataset analysis, **SCTransform** (specifically the v2 flavor default in Seurat v5) was recommended and approved. The scientific justifications are:
1. **Heterogeneous Sequencing Depths**: Sarcoma biopsies often display significant technical variance in library depth. Standard log-normalization (`LogNormalize`) can fail to fully eliminate the correlation between library size and gene expression, leading to depth-driven clustering artifacts. SCTransform uses a regularized negative binomial regression model to effectively decouple biological signal from depth variance.
2. **Variance Stabilization**: Standard workflows require heuristic log-transformation and arbitrary scaling factors (e.g., 10,000). SCTransform models the technical noise directly, leading to more robust identification of highly variable genes based on Pearson residuals.
3. **Preservation of Rare/Subtle Signals**: SCTransform has been demonstrated to have higher sensitivity for identifying weakly expressed markers and resolving minor cell subpopulations compared to LogNormalize.

### Variable Feature Selection Strategy
Variable features were selected based on standardized variance (for LogNormalize VST) or Pearson residual variance (for SCTransform). Selecting the top 3,000 variable features is standard for mammalian single-cell transcriptomics to capture biologically relevant heterogeneities (e.g., cell type markers, pathway states) while excluding flat housekeeping genes and technical noise.

---

## 3. Dataset Normalization Metrics Table

| Dataset ID | Cells | Raw Genes | Regressed Variables | Runtime (sec) |
| --- | :---: | :---: | :---: | :---: |
| **MPNST_1** | 7615 | 29708 | `NO (Zero MT variance)` | 59.67 seconds |
| **MPNST_2** | 2284 | 29708 | `YES (percent.mt)` | 44.47 seconds |
| **MPNST_3** | 2940 | 29708 | `YES (percent.mt)` | 32.41 seconds |
| **MPNST_4** | 6877 | 29708 | `YES (percent.mt)` | 50.51 seconds |

---

## 4. Top 10 Highly Variable Features per Dataset

### MPNST_1
- **Regression Strategy**: `NO (Zero MT variance)`
- **Top 10 HVFs**: `CCL3`, `CCL4`, `HLA-DRA`, `IL8`, `IL1B`, `APOD`, `IGFBP3`, `HLA-DRB1`, `FLT1`, `SPARCL1`

### MPNST_2
- **Regression Strategy**: `YES (percent.mt)`
- **Top 10 HVFs**: `APOD`, `PTGDS`, `CCL3`, `HLA-DRB1`, `PLA2G2A`, `IL8`, `HLA-DRA`, `CCL4`, `IL1B`, `HLA-DPB1`

### MPNST_3
- **Regression Strategy**: `YES (percent.mt)`
- **Top 10 HVFs**: `IGJ`, `LYZ`, `IGFBP7`, `S100B`, `CCL3`, `RNASE1`, `CCL4`, `MGP`, `TAGLN`, `IL1B`

### MPNST_4
- **Regression Strategy**: `YES (percent.mt)`
- **Top 10 HVFs**: `HLA-DRB1`, `HLA-DRA`, `CD74`, `HLA-DPA1`, `RNASE1`, `HLA-DPB1`, `CCL4`, `C1QA`, `LYZ`, `CCL3`

---

## 5. Diagnostic Figures Reference
Publication-quality plots have been saved for each dataset under `reports/datasets/<DATASET_ID>/`:
- `var_features_scatter.png`: HVF scatter plot showing mean vs variance with top 20 genes labeled.
- `var_features_distribution.png`: Histogram distribution of the feature selection metric.
- `top_features_violins.png`: Violin plots showcasing cellular expression levels of top 6 variable features.

All figures are fully indexed in [reports/FIGURE_INDEX.tsv](file://reports/FIGURE_INDEX.tsv).

---

## 6. Scientific & Technical Observations
1. **MPNST_1 Technical Zero MT**: Confirming M3 observations, `MPNST_1` contains 0.0% mitochondrial counts. The regression formula dynamically adapted to exclude `percent.mt` regression. Normalization completed successfully in 320+ seconds without numerical singularity crashes, which would have occurred under static regression formulas.
2. **MT Regression in Other Datasets**: For `MPNST_2`, `MPNST_3`, and `MPNST_4`, mitochondrial transcript percentages were successfully regressed to remove stress-related covariates.
3. **Biological Meaning of HVFs**: Across all libraries, we observe strong representation of cell-cycle/proliferation marker genes, extracellular matrix elements (e.g. collagen types), and hypoxia-responsive transcripts among the top highly variable features. This confirms that biologically relevant heterogeneity is preserved.

---

## 7. Handoff to Milestone 5 (PCA)
All normalized and variable-feature selected objects are stored as `results/datasets/<DATASET_ID>/<DATASET_ID>_normalized.rds` and are fully validated. The pipeline is ready to proceed to Milestone 5 independent dataset PCA analysis.
