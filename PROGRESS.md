# MPNST Phase 1 Project Progress

This document tracks the progress of the Phase 1 independent dataset analysis, reproducible workflow engineering, marker discovery, and pre-integration handoff.

| Milestone Status Summary

| Milestone | Description | Status | Target Date | Completion Date |
| --- | --- | --- | --- | --- |
| **M0** | Infrastructure, Environment, and SLURM Safety | **Completed** | 2026-07-06 | 2026-07-06 |
| **M1** | Real-Data Audit | **Completed** | 2026-07-07 | 2026-07-07 |
| **M2** | Dataset Extraction and Pre-Filter QC | Not Started | - | - |
| **M3** | QC Filtering and Doublet Assessment | Not Started | - | - |
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

