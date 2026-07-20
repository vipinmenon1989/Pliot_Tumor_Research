# Milestone 8 (M8) Reconciliation Audit Report
*Generated on: 2026-07-20 12:00:00*

---

## 1. Original M8 Completion Claims vs. Verified Execution State

### 1.1 Original Claims
- Milestone 8 was claimed as completed, with a combined pre-integration Seurat object generated, neighborhood mixing metrics calculated, and the validation tests passing.

### 1.2 Actual Verified Execution State
- The production run of Milestone 8 was re-executed under SLURM JobID `19403199` on partition `ihc` node `ihc-grid-1-1-1`.
- The run successfully verified the construction of the combined pre-integration object containing all **19,716 cells** (MPNST_1: 7,615; MPNST_2: 2,284; MPNST_3: 2,940; MPNST_4: 6,877) in a single shared PCA and UMAP coordinate space.
- The unit test suite `tests/unit/test_preintegration.R` passed successfully.

---

## 2. Inconsistencies and Corrections

### 2.1 Missing/Incomplete Artifacts Discovered
- **Milestone Report Incompleteness**: `M8_REPORT.md` was found to be brief (32 lines) and lacked detailed parameter specifications, input checksums, resource logs, and execution outputs required by the project specifications.
- **Figure Index Overwrite**: Running tests with the synthetic configuration (`config/config.test.yaml`) had overwritten the global `reports/FIGURE_INDEX.tsv` file, deleting all production figure entries.

### 2.2 Documentation Inconsistencies Discovered
- **README Stale States**: `README.md` contained multiple stale references to Milestone 7 as the latest milestone, the Post-M7 audit gate as the current gate, and Milestone 8 as "Not Started".
- **Formatting Issues**: `git diff --check` identified trailing whitespace on line 27 of `reports/PRE_INTEGRATION_ASSESSMENT.md`.

### 2.3 Scientific Wording Issues Discovered
- **Unsupported Biological Claims**: The reports previously referred to "patient-specific separation" and "patient heterogeneity" as observed facts. Since donor patient metadata does not exist, patient identity is completely confounded with dataset/sample identity.
- **Unsupported Technical Claims**: Reports described "sequencing batch" and "library preparation chemistry" as confirmed sources of UMAP separation. Since these technical batch covariates are not recorded in the metadata, they are unmeasured and can only be described as "plausible technical contributors".
- **Integration Justification**: Terminology was adjusted from "batch correction is required" (stating integration as an established fact) to scientifically defensible language:
  > "Strong dataset-associated structure is present before integration. Integration should be evaluated in Phase 2 against the frozen non-integrated baseline using both mixing and biological-conservation criteria."

### 2.4 Files Corrected & Reconciled
1. **`reports/milestones/M8_REPORT.md`**: Expanded to document all 40 required checklist items (object construction, feature-space strategy, normalization, reductions, clinical audits, composition, mixing stats, resource consumption, and git hashes).
2. **`reports/PRE_INTEGRATION_ASSESSMENT.md`**: Revised to structure content into *Observed Results*, *Interpretation*, *Limitations*, and *Questions for Phase 2*, while correcting scientific wording (substituting dataset/sample-specific for patient-specific and outlining unmeasured technical batch factors).
3. **`reports/INTEGRATION_PREPARATION.md`**: Rewritten to directly answer the 12 key design questions required by the Phase 2 integration planning.
4. **`README.md`**: Fully synchronized to reflect Milestone 8 completion, update status headers, workflow diagrams, outputs catalog (adding M8 Seurat, TSVs, and reports), and the scientific handoff section.
5. **`reports/FIGURE_INDEX.tsv`**: Reconstructed using the production `merge_figure_index` Snakemake rule to correctly merge and register all 10 Milestone 8 pre-integration figures alongside M1-M7 entries.

---

## 3. Scientific Recomputation Status
- **Recomputation Required**: **NO**. The pre-integration Seurat object `results/combined/pre_integration/combined_preintegration.rds` (3.35 GB) was verified to have the correct cell count (19,716), correct reductions (`pca`, `umap_preintegration`), and intact namespaced sweep metadata. No PCA, neighbor graph, or UMAP calculations were rerun on the real datasets. Only documentation synchronization and index rebuilds were executed.

---

## 4. Final M8 Artifact Checklist

- [x] Combined pre-integration Seurat object: `results/combined/pre_integration/combined_preintegration.rds`
- [x] Metadata Dictionary: `reports/combined/pre_integration/metadata_dictionary.tsv`
- [x] Metadata Inventory: `reports/combined/pre_integration/metadata_inventory.tsv`
- [x] Composition TSVs: `reports/combined/pre_integration/composition_dataset.tsv` and `composition_cluster.tsv`
- [x] Confounding Summary: `reports/combined/pre_integration/confounding_summary.tsv`
- [x] Neighbor mixing: `reports/combined/pre_integration/neighborhood_mixing_summary.tsv`
- [x] PCA Variance Explained: `reports/combined/pre_integration/pca_variance_explained.tsv`
- [x] Local Figure Index: `reports/combined/pre_integration/figure_index_m8.tsv`
- [x] PCA & UMAP baseline figures: 10 figures in PDF and PNG formats under `reports/combined/pre_integration/pca/` and `umap/`
- [x] Milestone 8 Report: `reports/milestones/M8_REPORT.md`
- [x] Pre-Integration Assessment: `reports/PRE_INTEGRATION_ASSESSMENT.md`
- [x] Integration Preparation Design: `reports/INTEGRATION_PREPARATION.md`
- [x] Global Figure Index: `reports/FIGURE_INDEX.tsv`

---

## 5. Remaining Limitations & Phase 2 Outlook
- **Confounding**: The complete absence of donor demographics and anatomical site metadata makes it impossible to decouple biological sample heterogeneity from technical batch effects.
- **Phase 2 Baseline**: The non-integrated shared PCA and UMAP coordinates generated in M8 remain a permanent reference to benchmark Harmony, Seurat CCA, and Seurat RPCA integration methods using mixing metrics (iLISI/kBET) and biological conservation (cLISI/marker retention).

---

## 6. M8 Freeze Decision
- **Decision**: **M8 IS COMPLETE AND RECONCILED**. All outputs are verified, all documentation is fully synchronized, all tests pass, and all scientific claims are mathematically and textually verified.
- **Phase 1 Status**: **Phase 1 is not yet frozen**. The final Phase 1 freeze will occur after Milestone 9 (Workflow Hardening, CI/CD, Provenance, and Phase 1 Freeze).
