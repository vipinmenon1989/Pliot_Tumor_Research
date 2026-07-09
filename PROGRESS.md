# MPNST Phase 1 Project Progress

This document tracks the progress of the Phase 1 independent dataset analysis, reproducible workflow engineering, marker discovery, and pre-integration handoff.

| Milestone Status Summary

| Milestone | Description | Status | Target Date | Completion Date |
| --- | --- | --- | --- | --- |
| **M0** | Infrastructure, Environment, and SLURM Safety | **Completed** | 2026-07-06 | 2026-07-06 |
| **M1** | Real-Data Audit | **Completed** | 2026-07-07 | 2026-07-07 |
| **M2** | Dataset Extraction and Pre-Filter QC | **Completed** | 2026-07-08 | 2026-07-08 |
| **M3** | QC Filtering and Doublet Assessment | **Completed** | 2026-07-09 | 2026-07-09 |
| **M4** | Normalization and Variable Features | Not Started | - | - |
| **M5** | PCA and PC Evaluation | Not Started | - | - |
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
  - M0 Cleanup and README Documentation Commit: [PENDING_README_HASH]

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


