# Milestone 8 (M8) Consolidated Execution Report
*Generated on: 2026-07-20 12:00:00*

---

## 1. Authorization & Context

### 1.1 Objective
The scientific objective of Milestone 8 (M8) is to construct a mathematically valid shared non-integrated expression space from the four Phase 1 constituent datasets (`MPNST_1`, `MPNST_2`, `MPNST_3`, and `MPNST_4`), followed by a single shared PCA and UMAP, establishing a ground-truth baseline reference before any integration or batch correction is applied in Phase 2.

### 1.2 Input Objects
The input objects are the four independently clustered Seurat objects from Milestone 6:
1. `results/datasets/MPNST_1/MPNST_1_clustered.rds`
2. `results/datasets/MPNST_2/MPNST_2_clustered.rds`
3. `results/datasets/MPNST_3/MPNST_3_clustered.rds`
4. `results/datasets/MPNST_4/MPNST_4_clustered.rds`

### 1.3 Input Checksums
The SHA256 checksums of the input objects:
- **MPNST_1**: `ba7894ac4ba7f576b816ce3854027e18e3d9a447bcbde0247fb2253c51df578f`
- **MPNST_2**: `595e893d813b340cc7ae56429233ee8ec7c27cd85885355fdb0a8a1dc174a94f`
- **MPNST_3**: `bf937b48be7c7699587fcaa511ba2d33f943616863b57aee58e18dd757480c96`
- **MPNST_4**: `3305ffa19fd6d89408bfb003fec4a3c21807cf19eaadad62a11a096593530896`

---

## 2. Combined Object Construction

### 2.1 Construction Strategy
We loaded the constituent clustered objects, stripped any legacy assay/reduction columns (except the raw count assays), checked cell names and prepended dataset ID (to avoid name collisions), merged the Seurat objects, namespaced metadata resolution sweep columns, and generated globally unique recommended/alternative cluster labels.

### 2.2 Expected vs. Actual Cell Counts
- **Expected Cell Count**: `19,716` cells
- **Actual Combined Cell Count**: `19,716` cells
- **Cell Loss**: 0 cells lost.
- **Duplicated Cell Names**: 0 duplicates found.

### 2.3 Feature-Space Strategy
We derived 3,000 highly variable features globally from the merged raw counts.

### 2.4 Assay/Layer Configuration
- **Assay**: `SCT` (default) containing `counts`, `data`, and `scale.data` layers.
- **Assay**: `RNA` containing raw `counts` and normalized `data` layers.

### 2.5 Normalization Compatibility Strategy
We executed a unified global `SCTransform` (v2) on the merged raw counts (regressing out `percent.mt` and using `vst.flavor = "v2"`) to calculate comparable Pearson residuals across all datasets.

---

## 3. Metadata Preservation

### 3.1 M6 Multi-Resolution Clustering Sweep
All ten frozen M6 clustering assignments (resolutions 0.1 to 1.0) were preserved in the combined object for each dataset using a namespaced format:
- Format: `preint_MPNST_{dataset}_res_{resolution}` (e.g., `preint_MPNST_1_res_0.1` to `preint_MPNST_4_res_1.0`).
- Cells from other datasets are filled with `NA` in these columns to prevent mixing.

### 3.2 M7 Recommendations
We preserved all M7-specific recommendations:
- `preint_M6_recommended_resolution`: Original computationally recommended resolution.
- `preint_M7_recommended_resolution`: Finalized recommended resolution.
- `preint_M7_alternative_resolution`: Proposed alternative resolution.
- `preint_recommended_cluster`: Harmonized recommended cluster IDs (globally unique, e.g., `MPNST_1_C00` ... `MPNST_4_C13`).
- `preint_alternative_cluster`: Harmonized alternative cluster IDs (globally unique).

---

## 4. Dimensional Reduction Methodology

### 4.1 Shared PCA
- **Methodology**: Standard PCA run on the 3,000 unified variable features using the `scale.data` layer of the combined `SCT` assay.
- **Dimensions**: PCs 1-30.
- **Reduction Name**: `pca`.

### 4.2 Shared Neighbors
- **Methodology**: SNN graph construction on PCs 1-30.
- **Reduction Name**: `combined_preintegration_snn`.

### 4.3 Shared UMAP
- **Methodology**: Native R `uwot` UMAP on PCs 1-30 using the cosine metric.
- **Parameters**: Seed `42`, default parameters.
- **Reduction Name**: `umap_preintegration`.
- **Confirmation**: A single coordinate system contains coordinates for all 19,716 cells. Dataset-specific UMAP coordinates were **NOT** concatenated.

---

## 5. Metadata Audit

### 5.1 Available Metadata
- **Provenance**: `orig.ident`, `sample_id` (representing the four dataset/sample batches).
- **Technical**: `nCount_RNA`, `nFeature_RNA`, `percent.mt`, `percent.ribo`, `doublet_class`, `doublet_score`.
- **Biological**: `orig.anno` (legacy annotations).

### 5.2 Unavailable Metadata
- **Patient Identity**: Not recorded / unavailable.
- **Tumor Metadata**: Tumor type, tumor subtype, primary/metastatic status, anatomical site, biological condition, and library/capture batch are not recorded in the input objects.

---

## 6. Scientific Assessment

### 6.1 Composition Findings
- **MPNST_1**: 7,615 cells (38.62%)
- **MPNST_2**: 2,284 cells (11.58%)
- **MPNST_3**: 2,940 cells (14.91%)
- **MPNST_4**: 6,877 cells (34.88%)

### 6.2 Neighborhood Mixing Diagnostics
We evaluated dataset mixing in the shared PCA space ($k=15$ neighbors):
- **Mean Same-Dataset Fraction**: **96.13%** (MPNST_1: 98.83%, MPNST_2: 91.02%, MPNST_3: 96.90%, MPNST_4: 97.79%).
- **Mean Neighborhood Dataset Entropy**: **0.0966**.

### 6.3 Confounding Analysis
Dataset identity (`sample_id`) is 100% confounded with patient/sample identity. We cannot separate biological variance from technical batch effects.

### 6.4 Technical Concerns
In `MPNST_4`, a sub-neighborhood is strongly aligned with high `percent.mt` (mitochondrial transcripts), suggesting potential cellular stress bias.

### 6.5 Biological Interpretation Limitations
No clinical metadata exists, meaning donor biological differences cannot be decoupled from potential technical library batches.

### 6.6 Integration-Readiness Conclusion
Strong dataset-associated structure is present before integration. Integration should be evaluated in Phase 2 against the frozen non-integrated baseline using both mixing and biological-conservation criteria.

---

## 7. Execution and Computational Resources

### 7.1 SLURM Job Details
- **JobID**: `19403199`
- **Requested Resources**: 8 CPUs, 64 GB RAM, partition `ihc`, node `ihc-grid-1-1-1`.
- **ExitCode**: `0:0` (Success)
- **Elapsed Time**: `17m 17s`
- **MaxRSS**: `~32.00 GB`

### 7.2 Validation Status
- **Snakemake validation status**: Passed (`logs/workflow/test_preintegration.done` target completed).
- **Tests run**: `tests/unit/test_preintegration.R` verified cell counts, cell name uniqueness, namespaced columns, dimensional reductions, and absence of integrated embeddings.

---

## 8. Deliverables & Outputs

### 8.1 Files Created
- Combined pre-integration Seurat object: `results/combined/pre_integration/combined_preintegration.rds` (SHA256: `c3fdce8b61602725977f989d8bbf10030ffe48ef6140372a13d2c59903a8dc66`)
- Metadata Dictionary: `reports/combined/pre_integration/metadata_dictionary.tsv`
- Metadata Inventory: `reports/combined/pre_integration/metadata_inventory.tsv`
- Composition TSVs: `reports/combined/pre_integration/composition_dataset.tsv`, `composition_cluster.tsv`
- Confounding Summary: `reports/combined/pre_integration/confounding_summary.tsv`
- Neighbor mixing: `reports/combined/pre_integration/neighborhood_mixing_summary.tsv`
- Figures: 10 PCA and UMAP figures saved under `reports/combined/pre_integration/pca/` and `umap/` (PDF and PNG)
- Narrative reports: `reports/PRE_INTEGRATION_ASSESSMENT.md`, `reports/INTEGRATION_PREPARATION.md`

### 8.2 Files Modified
- `workflow/Snakefile`
- `reports/FIGURE_INDEX.tsv`
- `PROGRESS.md`
- `CHANGELOG.md`
- `README.md`

---

## 9. Warnings and Limitations
- **Buffering Warning**: Standard output in R is buffered when redirected; log file updates may be delayed until chunk thresholds are reached or the script terminates.
- **NFS Caching**: Direct recursive file searches on NFS mounts may return stale or missing directory contents until an explicit metadata lookup (such as `ls`) invalidates the attribute cache.
- **No integration or batch correction has yet been performed.**
