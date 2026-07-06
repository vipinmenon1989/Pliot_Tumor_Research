# Milestone 0 (M0) Execution Report
**Infrastructure, Environment, and SLURM Safety**

- **Date**: 2026-07-06
- **Status**: Completed
- **Implementation Agent**: Antigravity
- **Initial Commit Hash**: [PENDING_COMMIT_HASH]
- **Documentation Commit Hash**: [PENDING_COMMIT_HASH]

---

## 1. Executive Summary

Milestone 0 has been completed successfully in full compliance with the repository architecture and safety restrictions outlined in `PROJECT.md`. The environment was successfully audited, Snakemake was safely installed, a synthetic Seurat dataset was generated, unit smoke tests were written and executed, and a trivial SLURM smoke job was run to verify cluster integration. The real input file (`processed_mpnst.rds`) was **not** loaded, inspected, deserialized, or otherwise accessed. All staged files have passed whitespace and safety verification scripts.

---

## 2. Environment Baseline & Snakemake Installation

### 2.1 Pre-installation Snapshot
- **Active Environment Location**: `/local/projects-t3/lilab/vmenon/anaconda3/envs/R_env`
- **R Version**: 4.4.3
- **Python Version**: 3.11.13
- **Critical Single-Cell Packages**:
  - `Seurat` (5.4.0)
  - `SeuratObject` (5.3.0)
  - `Matrix` (1.7.4)
  - `sctransform` (0.4.3)
  - `SingleCellExperiment` (1.28.1)
  - `SummarizedExperiment` (1.36.0)
  - `BiocGenerics` (0.52.0)
  - `future` (1.69.0)
  - `future.apply` (1.20.1)
  - `ggplot2` (4.0.1)
  - `patchwork` (1.3.2)
  - `data.table` (1.18.0)
  - `harmony` (1.2.4)
  - `DoubletFinder` (2.0.4)
  - `scDblFinder` (1.20.2)
- Saved Conda baseline specifications to [workflow/envs/R_env_baseline.yaml](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/workflow/envs/R_env_baseline.yaml) and [workflow/envs/R_env_baseline.explicit.txt](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/workflow/envs/R_env_baseline.explicit.txt).

### 2.2 Dry-Run Solve Analysis
- The dependency solve proposed installing `snakemake` (9.23.1) and python packages (`pulp`, `pydantic`, etc.).
- There were **no** proposed upgrades, downgrades, removals, or replacements of any of the critical R packages or single-cell packages listed above. The solve was determined to be safe.

### 2.3 Snakemake Installation & Verification
- **Command**: `conda install -y snakemake`
- **Verification Details**:
  - `which snakemake`: `/local/projects-t3/lilab/vmenon/anaconda3/envs/R_env/bin/snakemake`
  - `snakemake --version`: `9.23.1`
  - R startup and critical package loads were verified programmatically; all packages loaded successfully without regression.
- Saved post-installation specifications to [workflow/envs/R_env_final.yaml](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/workflow/envs/R_env_final.yaml), [workflow/envs/R_env_final.explicit.txt](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/workflow/envs/R_env_final.explicit.txt), and [workflow/envs/R_env_final_sessionInfo.txt](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/workflow/envs/R_env_final_sessionInfo.txt).

### 2.4 Conda Environment Diff Summary
A machine-readable comparison of the Conda environment before and after Snakemake installation was generated:
- **Packages Added**: 81 (including `snakemake-minimal`, `pulp`, `pydantic`, `jsonschema`)
- **Packages Removed**: 0
- **Packages Upgraded**: 2 (`ca-certificates` from 2026.1.4 to 2026.6.17, `openssl` from 3.6.0 to 3.6.3)
- **Packages Downgraded**: 0
- **Critical Packages**: Verified **unchanged** (all R and single-cell package versions were preserved).
- The detailed comparison reports are saved at [reports/audits/conda_env_diff.json](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/audits/conda_env_diff.json) and [reports/audits/conda_env_diff.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/audits/conda_env_diff.md).

### 2.5 Portable Environment Specification
A portable Conda specification was created by stripping the machine-specific prefix from the final environment yaml.
- Portable Specification: [workflow/envs/R_env_portable.yaml](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/workflow/envs/R_env_portable.yaml)
- Explanation of distinction between raw snapshots and portable specifications can be found in [workflow/envs/README.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/workflow/envs/README.md).

---

## 3. Workflow & Test Setup

- **Directory Skeleton**: Created config, workflow, scripts, tests, data, results, reports, logs, and benchmarks subdirectories.
- **Configurations**:
  - [config/config.yaml](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/config/config.yaml): Production settings using the real RDS.
  - [config/config.test.yaml](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/config/config.test.yaml): Test settings pointing to synthetic data.
  - [config/schemas/config.schema.yaml](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/config/schemas/config.schema.yaml): Configuration structure validator.
- **Logging & Provenance**: Implemented [scripts/R/logging_utils.R](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/scripts/R/logging_utils.R) and [scripts/R/provenance_utils.R](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/scripts/R/provenance_utils.R) to enforce structured stdout/stderr logging and compute file MD5 checksums for reproducible data lineage.
- **Deterministic Synthetic Data**:
  - Script: [scripts/R/generate_synthetic_data.R](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/scripts/R/generate_synthetic_data.R)
  - Outputs: Created [data/synthetic/synthetic_mpnst.rds](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/data/synthetic/synthetic_mpnst.rds) (120 genes, 600 cells, split across 4 samples, 2 patients, 2 conditions, and 3 cell types with simulated technical depth/MT effects and known markers).
  - Provenance saved to [reports/audits/synthetic_data_provenance.json](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/audits/synthetic_data_provenance.json).
- **Unit/Smoke Tests**:
  - Script: [tests/unit/test_smoke.R](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/tests/unit/test_smoke.R). Verified configuration parser, loaded Seurat, and confirmed synthetic object metadata. All tests passed.
- **Snakemake Workflow**:
  - Snakefile: [workflow/Snakefile](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/workflow/Snakefile).
  - Validated local execution of `validate_config` and `smoke_test` rules via `snakemake --cores 1 --configfile config/config.test.yaml`. All rules executed successfully.
  - Successfully verified synthetic-data reproducibility by deleting generated outputs (`data/synthetic/synthetic_mpnst.rds` and completion flags) and regenerating them cleanly through the production Snakemake workflow.

---

## 4. SLURM Smoke Job Submission

A trivial SLURM smoke job was submitted to verify scheduler integration and safety:
- **Job Script**: [tests/integration/submit_smoke.sh](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/tests/integration/submit_smoke.sh)
- **JobID**: `19160128`
- **Partition**: `ihc` (account: `ihc`, node: `ihc-grid-1-1-1`)
- **Requested Resources**: 1 CPU, 4GB RAM, 10 min walltime
- **Actual Resource Usage**:
  - **State**: COMPLETED
  - **ExitCode**: 0:0
  - **Elapsed Time**: 20 seconds
  - **MaxRSS**: 0 (lightweight job finished quickly)
- **Log Verification**:
  - Output log [logs/slurm/smoke_job_19160128.out](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/logs/slurm/smoke_job_19160128.out) confirms R, python, and Snakemake paths are correctly set on the compute node and that `Seurat` and `SeuratObject` load successfully.
  - Error log [logs/slurm/smoke_job_19160128.err](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/logs/slurm/smoke_job_19160128.err) contains only standard Seurat startup attachments.

---

## 5. Summary of Workspace Changes

### 5.1 Files Created
1. [.gitignore](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/.gitignore) - Excludes large objects, logs, and figures from Git.
2. [PROGRESS.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/PROGRESS.md) - Milestone tracking sheet.
3. [config/config.yaml](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/config/config.yaml) - Production parameters configuration.
4. [config/config.test.yaml](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/config/config.test.yaml) - Test parameters configuration.
5. [config/schemas/config.schema.yaml](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/config/schemas/config.schema.yaml) - JSON schema for config validation.
6. [workflow/Snakefile](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/workflow/Snakefile) - Snakemake workflow rules.
7. [workflow/profiles/slurm/config.yaml](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/workflow/profiles/slurm/config.yaml) - Snakemake SLURM execution profile.
8. [scripts/R/logging_utils.R](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/scripts/R/logging_utils.R) - Standardized logging helper.
9. [scripts/R/provenance_utils.R](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/scripts/R/provenance_utils.R) - Reproducibility and checksum tracker.
10. [scripts/R/generate_synthetic_data.R](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/scripts/R/generate_synthetic_data.R) - Deterministic synthetic Seurat data generator.
11. [scripts/python/conda_diff.py](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/scripts/python/conda_diff.py) - Conda comparison generator.
12. [scripts/shell/verify_milestone.sh](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/scripts/shell/verify_milestone.sh) - Milestone verification utility.
13. [tests/unit/test_smoke.R](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/tests/unit/test_smoke.R) - Environment unit test script.
14. [tests/integration/submit_smoke.sh](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/tests/integration/submit_smoke.sh) - SLURM test script.
15. [reports/milestones/M0_REPORT.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/milestones/M0_REPORT.md) - This report.

### 5.2 Commands Executed
- Directory skeleton creation: `mkdir -p ...`
- Snapshot captures: `conda env export`, `conda list --explicit`, `conda info`, `conda list`
- Dry-run solver: `conda install --dry-run snakemake`
- Installation: `conda install -y snakemake`
- Verify Snakemake and R: `which snakemake`, `snakemake --version`, `Rscript -e '...'`
- Execute Synthetic generator: `Rscript scripts/R/generate_synthetic_data.R`
- Run unit test: `Rscript tests/unit/test_smoke.R`
- Snakemake dry run: `snakemake -n --configfile config/config.test.yaml`
- Snakemake run: `snakemake --cores 1 --configfile config/config.test.yaml`
- Submit SLURM: `sbatch tests/integration/submit_smoke.sh`
- Query accounting: `squeue -j ...`, `sacct -j ...`

### 5.3 Commit-Ready Changes
- All files listed in Section 5.1 (except those ignored by `.gitignore`) are staged/ready for commit.
- Ignored paths: `processed_mpnst.rds`, `results/`, `data/`, `logs/`, `benchmarks/`, `.snakemake/`

---

## 6. Unresolved Issues and Decisions for Researcher Review

1. **Constituent Dataset Identifier**: We have set `"sample_id"` as the default `dataset_id_column` in the config, and configured `"orig.ident"` as an alternative. The cardinality, cell distribution, and patient mappings must be audited during Milestone 1 (M1) to finalize which metadata column correctly isolates constituent datasets.
