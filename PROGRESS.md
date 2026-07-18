# MPNST Phase 1 Project Progress

This document tracks the progress of the Phase 1 independent dataset analysis, reproducible workflow engineering, marker discovery, and pre-integration handoff.

| Milestone Status Summary

| Milestone | Description | Status | Target Date | Completion Date |
| --- | --- | --- | --- | --- |
| **M0** | Infrastructure, Environment, and SLURM Safety | **Completed** | 2026-07-06 | 2026-07-06 |
| **M1** | Real-Data Audit | **Completed** | 2026-07-07 | 2026-07-07 |
| **M2** | Dataset Extraction and Pre-Filter QC | **Completed** | 2026-07-08 | 2026-07-08 |
| **M3** | QC Filtering and Doublet Assessment | **Completed** | 2026-07-09 | 2026-07-09 |
| **M4** | Normalization and Variable Features | **Completed** | 2026-07-09 | 2026-07-09 |
| **M5** | PCA and PC Evaluation | **Completed** | 2026-07-10 | 2026-07-10 |
| **M6** | Clustering Resolution Sweep | Not Started | - | - |
| **M7** | Marker Discovery and Dataset Recommendations | Not Started | - | - |
| **M8** | Combined Pre-Integration Baseline | Not Started | - | - |
| **M9** | Workflow Hardening, CI/CD, Provenance, and Phase 1 Freeze | Not Started | - | - |

---

## Detailed Milestone Logs

### M0 — Infrastructure, Environment, and SLURM Safety
- **Status**: Completed (2026-07-06)
- **Conda Environment**: Successfully installed Snakemake without affecting critical packages. Generated a portable env specification `workflow/envs/R_env_portable.yaml` and a comparison diff.
- **Repository Setup**: Initialized directory skeleton. Moved verification checks to `scripts/shell/verify_milestone.sh`.
- **Configuration**: Created base config files and schema validation.
- **Synthetic Data**: Generated and verified synthetic data via Snakemake workflow execution from clean state.
- **SLURM Integration**: Submitted a trivial smoke job (JobID `19160128`) on `ihc` partition and verified successful execution log.
- **Whitespace / Quality Check**: Cleaned all trailing whitespaces. Verified that `git diff --cached --check` passes.
- **Git Commits**:
  - Initial M0 Implementation Commit: f6da0c537dc618925b6fed673a249e2601599d89
  - Documentation Commit: e1ab6d11c4234c498728db47c5b0e7962c6d5a94
  - M0 Cleanup and README Documentation Commit: 89ca4ea3aacae9d5a5da60f73c1a8e358fec7ae3

### M1 — Real-Data Audit
- **Status**: Completed & Extended to Complete Scientific Inventory (2026-07-07)
- **Codebase Update**: Developed `scripts/R/audit_object.R` to run read-only scientific audits, supporting both synthetic and split-layer Seurat v5 structures. Developed `scripts/R/inspect_additional.R` to retrieve advanced object statistics (memory size, duplicate genes/cells, and mitochondrial/ribosomal gene lists). Developed `tests/unit/test_audit.R` to validate audit outputs.
- **Synthetic Validation**: Ran synthetic data validation through local Snakemake execution. All unit tests passed successfully.
- **SLURM Production Run**: Submitted the real dataset audit via sbatch (JobID `19160528`) on partition `ihc` node `ihc-grid-1-1-1`. Run additional inspection via JobID `19160550`.
- **HPC Execution Metrics**:
  - State: COMPLETED (ExitCode 0:0)
  - Elapsed: 00:01:32 (SLURM job), 76.94 seconds (Snakemake rule)
  - MaxRSS: 9627268K (~9.18 GB)
  - CPU Usage: mean_load 54.72%, cpu_time 42.23s
  - Estimated object memory footprint: 10.08 GB
- **Key Scientific Findings**:
  - Object dimensions: 29,708 features x 22,661 cells in default SCT assay (RNA assay contains 31,764 features x 22,661 cells).
  - Raw counts: Verified split raw RNA counts exist in the Seurat object (split layers `counts.MPNST_1` through `counts.MPNST_4` under `RNA` assay). Reprocessing from raw UMI counts is scientifically feasible.
  - Dataset Identifier: Recommended using `sample_id` as the primary dataset splitter. Both `sample_id` and `orig.ident` have identical cardinality (4 samples: MPNST_1 (8338 cells), MPNST_2 (2830 cells), MPNST_3 (3682 cells), MPNST_4 (7811 cells)) and map 100% identically (0 mismatches).
  - Feature Statistics: 0 duplicated genes, 0 duplicated cell names. 13 mitochondrial genes matching `^MT-` and 433 ribosomal genes matching `^RP[SL]` identified.
  - Legacy Inventory: Successfully logged all legacy clusters (SCT, Harmony, CCA, RPCA, MNN) and dimensional reductions. Legacy columns have been marked for exclusion from Phase 1 processing to prevent data leakage.
- **Artifacts Generated / Expanded**:
  - [DATA_AUDIT.md](file://reports/DATA_AUDIT.md) (Expanded with full evidence tables and all 12 requested sections)
  - [DATASET_STRUCTURE.tsv](file://reports/DATASET_STRUCTURE.tsv)
  - [object_inventory.json](file://reports/object_inventory.json) (Expanded with all machine-readable metrics)
  - [object_inventory.tsv](file://reports/object_inventory.tsv) (Expanded with recommended usages)
  - [M1_REPORT.md](file://reports/milestones/M1_REPORT.md) (Revised milestone report)
  - [m1_audit_provenance.json](file://reports/audits/m1_audit_provenance.json)

### M2 — Dataset Extraction and Pre-Filter QC
- **Status**: Completed (2026-07-08)
- **Codebase Update**: Developed `scripts/R/extract_dataset.R` to split the Seurat object by `sample_id` and generate validation reports/manifests, and updated it with a robust counts checker. Developed `scripts/R/generate_qc_plots.R` to produce 4 types of publication-quality diagnostic plots (PDF and PNG) per sample and compile local figure indices. Developed `scripts/R/generate_qc_recommendation.R` to assess QC metrics against thresholds and write markdown reports. Added unit test `tests/unit/test_extraction.R` to validate extracted object integrity.
- **Workflow Integration**: Updated `workflow/Snakefile` with rules for dataset extraction, plotting, recommendations, figure index merging, and extraction validation tests. Enabled sample list configuration in `config.yaml` and `config.test.yaml`.
- **Synthetic Validation**: Successfully executed synthetic tests locally, verifying that all rules and tests run to completion and pass.
- **SLURM Production Run**: Submitted the production run to SLURM (JobID `19161116`) on partition `ihc` node `ihc-grid-1-1-1`.
- **HPC Execution Metrics**:
  - State: COMPLETED (ExitCode 0:0)
  - Elapsed: 00:06:35
  - MaxRSS: 66935856K (~63.83 GB) — Efficiency of 99.7% of the 64 GB limit.
- **Key Scientific Findings**:
  - Successfully extracted four datasets (`MPNST_1` to `MPNST_4`).
  - Discovered that `MPNST_1` has 0.00% mitochondrial transcripts across all 8,338 cells.
  - Calculated percent.ribo across all cells using `^RP[SL]`.
  - Assessed expected filtering impact of default thresholds (combined filter excludes 16.21% for MPNST_1, 8.02% for MPNST_2, 5.43% for MPNST_3, and 4.07% for MPNST_4).
- **Artifacts Generated**:
  - `results/datasets/MPNST_*/` (Split Seurat objects and extraction provenance logs)
  - `reports/datasets/MPNST_*/DATASET_VALIDATION.md`
  - `reports/datasets/MPNST_*/QC_RECOMMENDATION.md`
  - `reports/datasets/MPNST_*/manifest.json`
  - `reports/datasets/MPNST_*/` QC plots (violins, scatter, histograms, metrics in PDF and PNG) and metrics summary TSVs
  - `reports/FIGURE_INDEX.tsv`
  - `reports/milestones/M2_REPORT.md`

### M3 — QC Filtering and Doublet Assessment
- **Status**: Completed (2026-07-09)
- **Codebase Update**:
  - Developed `scripts/R/filter_and_detect_doublets.R` to run scDblFinder doublet detection and apply QC filters, supporting a custom output suffix to enable dual-strategy evaluation.
  - Developed `scripts/R/qc_optimization_review.R` to run threshold sensitivity simulations and grid intersection overlap plots.
  - Developed `scripts/R/compare_qc_strategies.R` to quantitatively and visually compare global vs dataset-specific strategies.
  - Created unit tests `tests/unit/test_filtering.R` and `tests/unit/test_filtering_specific.R` to check compliance for both strategies.
- **Workflow Integration**: Updated `workflow/Snakefile` with rules for both Global (Strategy A) and Dataset-Specific (Strategy B) filtering, optimization review, comparison analysis, unit tests, and global figure index merging. Added wildcard constraints to prevent filename ambiguities.
- **Synthetic Validation**: Successfully executed synthetic tests locally, verifying that all rules and tests run to completion and pass.
- **SLURM Production Runs**:
  - Global Run: Submitted to SLURM (JobID `19161422` on partition `ihc` node `ihc-grid-1-1-1`, COMPLETED, 3m 24s, MaxRSS 36.88 GB).
  - Specific & Comparison Run: Submitted to SLURM (JobID `19161426` on partition `ihc` node `ihc-grid-1-1-1`, COMPLETED, 9m 40s, MaxRSS 63.67 GB).
- **Key Scientific Findings**:
  - Flat global thresholds lead to a catastrophic cell loss of **46.9%** in `MPNST_2`, **48.3%** in `MPNST_3`, and **56.8%** in `MPNST_4`.
  - Under proposed dataset-specific thresholds, cell retention increases to **80.7%** (`MPNST_2`), **79.8%** (`MPNST_3`), and **88.0%** (`MPNST_4`), avoiding artificial truncation of biological expression profiles (such as ribosomal and mitochondrial fractions naturally elevated in sarcomas) while successfully removing doublets (~8.0-9.2% of cells) and low-complexity cells.
- **Artifacts Generated**:
  - `results/datasets/MPNST_*/MPNST_*_filtered.rds` & `MPNST_*_filtered_specific.rds` (Filtered Seurat objects for both strategies)
  - `reports/datasets/MPNST_*/FILTER_REPORT.md` & `FILTER_REPORT_SPECIFIC.md` (Dataset filter reports)
  - `reports/datasets/MPNST_*/manifest_filtered.json` & `manifest_filtered_specific.json`
  - `reports/datasets/MPNST_*/` post-filter QC plots (violins, scatter, density, histograms, doublet summaries, filtering summaries for both strategies)
  - `reports/qc_optimization/` plots and statistics (distribution comparisons, sensitivity curves, and grid intersection overlaps)
  - `reports/qc_comparison/` plots and tables (retention barplot, mt/ribo violin comparison, composite distributions, and strategy comparison summaries)
  - `reports/QC_OPTIMIZATION_REPORT.md` and `reports/QC_COMPARISON_REPORT.md`
  - `reports/FIGURE_INDEX.tsv` (Authoritative merged index of all 32 generated figure paths)
  - `reports/milestones/M3_REPORT.md`

### M4 — Normalization and Variable Features
- **Status**: Completed (2026-07-09)
- **Codebase Update**:
  - Developed `scripts/R/normalize_and_find_features.R` to run Seurat v5 SCTransform v2 or LogNormalize dynamically, handling zero-mitochondrial libraries (like `MPNST_1`) and setting parallelization global limits options for the `future` package.
  - Developed `scripts/python/generate_m4_report.py` to compile normalization runtimes and HVFs into a consolidated summary report.
  - Created unit test `tests/unit/test_normalization.R` to validate normalized objects and diagnostic outputs.
- **Workflow Integration**: Extended `workflow/Snakefile` with rules `normalize_and_find_features`, `test_normalization`, and `generate_m4_report` using dataset-specific QC inputs (`*_filtered_specific.rds`). Updated global figure index rules.
- **Synthetic Validation**: Successfully validated the normalization workflow on synthetic datasets locally.
- **SLURM Production Run**: Submitted the production run to SLURM (JobID `19174235`) on partition `ihc` node `ihc-grid-1-1-1`.
- **HPC Execution Metrics**:
  - State: COMPLETED (ExitCode 0:0)
  - Elapsed: 00:04:17
  - MaxRSS: 25287072K (~24.12 GB) — Highly efficient memory utilization under SLURM.
- **Key Scientific Findings**:
  - Standardized SCTransform v2 normalization decoupled sequencing depth covariance across MPNST constituent libraries.
  - `MPNST_1` had no mitochondrial transcripts, and the workflow dynamically bypassed mitochondrial regression, completing successfully.
  - Highly variable features (HVFs) across libraries successfully captured sarcoma-related biology, including extracellular matrix elements (collagens, APOD), cell-cycle markers, and macrophage-associated chemokine markers (CCL3, CCL4).
- **Artifacts Generated**:
  - `results/datasets/MPNST_*/MPNST_*_normalized.rds` (Normalized Seurat objects)
  - `reports/datasets/MPNST_*/NORM_REPORT.md` (Dataset-specific normalization reports)
  - `reports/datasets/MPNST_*/variable_features.tsv` (List of 3,000 highly variable features with metrics)
  - `reports/datasets/MPNST_*/` diagnostic plots (scatter, distribution, and top 6 expression violins in PDF and PNG)
  - `results/datasets/MPNST_*/normalization_provenance.json` (Digests and session specifications)
  - `reports/milestones/M4_REPORT.md` (Consolidated Milestone 4 report)

### M5 — PCA and PC Evaluation
- **Status**: Completed (2026-07-10)
- **Codebase Update**:
  - Developed `scripts/R/run_pca_and_evaluation.R` to execute PCA on normalized SCT assays, calculate variance explained relative to z-scored residuals, test cell score correlations with technical covariates, find PC selection ranges (conservative, recommended, maximum) using a geometric elbow knee-point detector, and save detailed results.
  - Developed `scripts/python/generate_m5_report.py` to aggregate dataset-specific metrics, generate comparison tables, and write `reports/PCA_RECOMMENDATIONS.tsv` and `reports/milestones/M5_REPORT.md`.
  - Created validation unit tests in `tests/unit/test_pca.R` to verify PCA reductions, coordinate dimensions, loadings, reports, and provenance.
- **Workflow Integration**: Extended `workflow/Snakefile` with rules `run_pca_and_evaluation`, `test_pca`, and `generate_m5_report`. Updated the figure index merging rule.
- **Synthetic Validation**: Validated local Snakemake execution on synthetic datasets, checking test parameter parsing.
- **SLURM Production Run**: Submitted workflow execution to SLURM (JobID `19176730`) on partition `ihc` node `ihc-grid-1-1-1`.
- **HPC Execution Metrics**:
  - State: COMPLETED (ExitCode 0:0)
  - Elapsed: 00:03:51
  - MaxRSS: 11410872K (~10.88 GB)
- **Key Scientific Findings**:
  - PCA geometric elbows identified optimal recommended PC cutoffs: `PC8` (MPNST_1), `PC6` (MPNST_2), `PC9` (MPNST_3), and `PC5` (MPNST_4).
  - SCTransform normalization decoupled sequencing depth covariance, resulting in low correlations (R < 0.2) in leading PCs.
  - Leading PCs are heavily dominated by biological programs: antigen presentation/immune response (CD74, HLA-DRA, HLA-DRB1) and extracellular matrix structure/remodeling (COL1A1, COL1A2, COL3A1, DCN, SFRP2), representing core MPNST biology.
  - JackStraw analysis was omitted due to lack of statistical validity on regularized negative binomial SCT z-scored residuals and high computational footprint.
- **Artifacts Generated**:
  - `results/datasets/MPNST_*/MPNST_*_pca.rds` (PCA-embedded Seurat objects)
  - `reports/datasets/MPNST_*/PCA_REPORT.md` (Dataset-specific PCA reports)
  - `reports/datasets/MPNST_*/top_loading_genes.tsv` (Top loading gene lists)
  - `reports/datasets/MPNST_*/pc_technical_correlations.tsv` (Correlation metrics)
  - `reports/datasets/MPNST_*/pca_variance_explained.tsv` (Raw variance scores)
  - `reports/datasets/MPNST_*/` diagnostic plots (elbow, cumulative variance, PC loadings, PC heatmaps, and technical correlations in PDF and PNG)
  - `reports/PCA_RECOMMENDATIONS.tsv` (Machine-readable recommendations)
  - `reports/milestones/M5_REPORT.md` (Consolidated Milestone 5 report)
  - `reports/FIGURE_INDEX.tsv` (Updated global figure index with all 40 PCA plots)



