# Phase 1 Final Handoff Report (Milestones M0–M9)
*Generated on: 2026-07-20 12:00:00*

---

## 1. Executive Summary

### 1.1 Objective
The primary objective of Phase 1 is to establish a rigorous, reproducible, and standardized pre-integration single-cell transcriptomic baseline for four independently processed MPNST datasets. This baseline serves as a frozen ground truth to benchmark and evaluate multi-sample integration, batch correction, and downstream biological annotation in Phase 2.

### 1.2 Scope
Phase 1 covers environment setup, data extraction, quality control filtering, doublet detection, normalization, principal component analysis, independent clustering sweeps, marker gene discovery, pre-integration combining, and workflow hardening. It explicitly excludes multi-sample batch correction, integration, cell-type annotation, condition-level differential expression, pathway analysis, and trajectory inference.

---

## 2. Methodology and Scientific Decisions

### 2.1 Input Data & Datasets
- **Original Input**: A combined clinical Seurat object containing clinical samples (`processed_mpnst.rds`).
- **Constituent Datasets**:
  1. `MPNST_1` (7,615 cells): Split and extracted.
  2. `MPNST_2` (2,284 cells): Split and extracted.
  3. `MPNST_3` (2,940 cells): Split and extracted.
  4. `MPNST_4` (6,877 cells): Split and extracted.

### 2.2 Quality Control (QC) Strategy
- **Min Features**: 500 features detected per cell.
- **Max Features**: 99,999,999 (no upper limit filtering applied).
- **Min Counts**: 1,000 counts per cell.
- **Max Counts**: 99,999,999 (no upper limit filtering applied).
- **Max Percent Mitochondrial Transcripts**: 15.0% global limit.
- **Max Percent Ribosomal Transcripts**: 20.0% global limit.

### 2.3 Doublet Detection Strategy
- **Method**: `scDblFinder` (v2.20.2) run independently on raw count matrices of each dataset.
- **Doublets Removed**: All cells classified as "doublet" by scDblFinder were purged.
- **Preserved Cells**: Exactly **19,716 cells** survived QC and doublet filtering.

### 2.4 Normalization and Feature Selection Strategy
- **Method**: SCTransform (v2) run globally/independently on raw counts (regressing out `percent.mt` and using `vst.flavor = "v2"`).
- **Variable Features**: Top 3,000 highly variable features selected.

### 2.5 PCA Recommendations
PCs were evaluated using variance-explained elbow plots, cumulative variance thresholds, and technical QC correlations:
- `MPNST_1`: Recommended **PCs 1–8** (Conservative: 1–8, Max: 1–30)
- `MPNST_2`: Recommended **PCs 1–6** (Conservative: 1–6, Max: 1–30)
- `MPNST_3`: Recommended **PCs 1–9** (Conservative: 1–9, Max: 1–30)
- `MPNST_4`: Recommended **PCs 1–5** (Conservative: 1–5, Max: 1–30)

### 2.6 Clustering sweep strategy
- **Algorithm**: Standard Louvain clustering.
- **Sweep Range**: Resolution 0.1 to 1.0 (step size of 0.1).
- **Evaluation**: Stability assessed using Clustree, Adjusted Rand Index (ARI) of subsampled resamples (prop = 0.8, n_resamples = 100), cluster sizes, and correlation with continuous covariates.

### 2.7 Final Clustering and Marker Recommendations
Markers were detected using Wilcoxon rank-sum tests (filtered by adjusted p-value < 0.05 and log2FC > 0.25):
- `MPNST_1` (res **0.6**): 18 clusters resolved.
- `MPNST_2` (res **0.3**): 9 clusters resolved.
- `MPNST_3` (res **0.6**): 13 clusters resolved.
- `MPNST_4` (res **0.7**): 14 clusters resolved. (Note: Alternative res **0.5** recommended if technical mitochondrial correlation needs to be minimized).

---

## 3. Pre-Integration Combined Baseline (Milestone 8)

### 3.1 Combined Non-Integrated Baseline
A single mathematically valid shared expression representation of **19,716 cells** was constructed. We executed a single unified `SCTransform` globally on the merged raw counts (selecting 3,000 variable features) and computed a shared PCA (PCs 1–30), neighbor graph ($k=15$), and UMAP (`umap_preintegration`) without any batch alignment or coordinate concatenation.

### 3.2 Metadata Preservation
- **Multi-Resolution Clustering Sweep**: Columns `preint_MPNST_{dataset}_res_{resolution}` store the sweep coordinates (cells from other datasets filled with `NA`).
- **M7 Recommendations**: Recommended cluster IDs are stored in `preint_recommended_cluster` and alternative IDs in `preint_alternative_cluster` (harmonized to globally unique formats, e.g., `MPNST_1_C00` ... `MPNST_4_C13`).

### 3.3 Composition, Mixing, and Confounding Findings
- **Composition**: MPNST_1 (38.62%), MPNST_2 (11.58%), MPNST_3 (14.91%), MPNST_4 (34.88%).
- **Neighborhood Mixing**: The mean same-dataset neighbor fraction in the shared PCA space ($k=15$) is **96.13%**, confirming that cells cluster almost exclusively by their sample/dataset origin.
- **Confounding**: Dataset identity is 100% confounded with patient identity. Unmeasured variables (such as capture channel, sequencing run, or library prep) cannot be decoupled from biological donor differences.
- **MPNST_4 Limitations**: Louvain cluster C09 shows strong technical correlation with percent.mt ($R^2 = 0.47$), representing potential cellular stress confounding.

### 3.4 Integration-Readiness Conclusion
Strong dataset-associated structure is present before integration. Integration should be evaluated in Phase 2 against the frozen non-integrated baseline using both mixing and biological-conservation criteria.

---

## 4. Phase 1 → Phase 2 Entry Contract

> [!IMPORTANT]
> To preserve workflow integrity, downstream analyses in Phase 2 must adhere to the following rules:
> 1. **Immutable Phase 1 Metadata**: Never overwrite or modify `preint_*` metadata columns, PCA/UMAP coordinates, or original inputs.
> 2. **Namespace Isolation**: All post-integration reductions and metadata columns must use a distinct namespace, e.g., prefixing them with `postint_*` (e.g., `postint_umap`, `postint_harmony_clusters`).
> 3. **Input Manifest**: Phase 2 must consume [results/phase1_manifest.json](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/results/phase1_manifest.json) as the machine-readable input.

---

## 5. Frozen Deliverables and Artifacts

The final targets of Phase 1:
1. **Manifest File**: [results/phase1_manifest.json](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/results/phase1_manifest.json)
2. **Combined Object**: [results/combined/pre_integration/combined_preintegration.rds](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/results/combined/pre_integration/combined_preintegration.rds)
3. **Dataset Objects**: `results/datasets/{ds}/{ds}_clustered.rds`
4. **Milestone Reports**: Reports `M0_REPORT.md` through `M9_REPORT.md` located under `reports/milestones/`.
5. **Preflight Checker**: [scripts/python/preflight_checker.py](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/scripts/python/preflight_checker.py)
6. **Project State Validator**: [scripts/python/validate_project_state.py](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/scripts/python/validate_project_state.py)
7. **Reconciliation Audits**: Reports `M8_RECONCILIATION.md` and `M0_M8_RECONCILIATION.md` under `reports/audits/`.

---

## 6. Reproducibility Instructions

### 6.1 Conda Environment
To create the environment:
```bash
conda env create -f workflow/envs/R_env_portable.yaml -n mpnst_phase1
conda activate mpnst_phase1
```

### 6.2 Preflight Check
Before running any Snakemake pipeline, run the preflight checker:
```bash
python scripts/python/preflight_checker.py
```

### 6.3 Running the Synthetic Workflow (Clean-Room / CI)
The synthetic workflow runs on reduced, mock datasets to ensure script validity.
To run the dry-run:
```bash
snakemake -n --configfile config/config.test.yaml
```
To run the clean-room synthetic workflow:
```bash
mkdir -p test_clean_room && cd test_clean_room
ln -s ../scripts scripts
ln -s ../tests tests
ln -s ../config config
ln -s ../workflow workflow
snakemake -s workflow/Snakefile --configfile config/config.test.yaml --cores 4 synthetic_complete --rerun-triggers mtime
```

### 6.4 Running the Production Workflow
Production runs utilize SLURM on the cluster. The wrapper script handles scheduling:
```bash
sbatch scripts/shell/run_m8_workflow.sh
```

---

## 7. Known Limitations

1. **Clinical Confounding**: Due to the complete absence of clinical covariates (anatomical site, treatment status, demographics), biological sample variation and technical batch effects are fully confounded.
2. **Technical Variables**: Continuous covariates like library size (`nCount_RNA`) and mitochondrial transcript fraction (`percent.mt`) represent the only technical metrics. Sequencing run dates or Chromium channel IDs are not recorded.
3. **Malignant Phenotypes**: Cell annotations (`orig.anno`) are imports and have not been biological-validated with copy-number variant (CNV) profiling or lineage markers in Phase 1.
