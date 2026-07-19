# MPNST Phase 1 Analysis Workflow

Reproducible, configuration-driven single-cell RNA-seq preprocessing and dataset auditing workflow for MPNST datasets.

---

## 1. Project Purpose & Scientific Objective
The primary objective of this project is to construct a modular, reproducible, configuration-driven, and HPC-ready Phase 1 workflow to preprocess and analyze constituent MPNST (malignant peripheral nerve sheath tumor) datasets independently. Crucially, all datasets are processed and evaluated in isolation before any batch correction or multi-dataset integration is performed. This prevents technical confounders from distorting early-stage QC and dimensionality reduction.

---

## 2. Current Project Status
- **Current Phase**: Phase 1 (Independent Dataset Processing)
- **Status**: **Completed: Milestones M0–M6**
- **Next Step**: **Milestone M7 — Marker Discovery and Dataset Recommendations**

> [!IMPORTANT]
> **Strict Phase 1 Scope Constraint**:
> This workflow strictly limits processing to independent dataset analysis. Downstream operations such as integration methods (Harmony, Seurat CCA/RPCA, FastMNN, scVI, BBKNN, etc.), automated cell-type annotation, condition-level differential expression, pathway analysis, trajectory analysis, cell-cell communication, and predictive modeling are strictly prohibited in Phase 1.

---

## 3. Workflow Overview

```text
M0 Infrastructure (Environment, synthetic data generation, SLURM verification)
    ↓
M1 Real-data audit (Authoritative object structure and inventory checks)
    ↓
M2 Dataset extraction + pre-filter QC (Deterministic extraction & diagnostic plotting)
    ↓
M3 QC filtering + doublet assessment (Doublet detection & dataset-specific thresholds)
    ↓
M4 Normalization + variable features (SCTransform v2 variance stabilization & HVF selection)
    ↓
M5 PCA + PC evaluation (Independent PCA & geometric elbow selection)
    ↓
M6 Clustering resolution sweep (SNN graph construction & sweeps 0.1–1.0)
    ↓
M7 Marker discovery + recommendations (Resolution-specific markers & selection)  <-- [NEXT STEP]
    ↓
M8 Combined pre-integration baseline (Consolidated baseline & integration prep)
    ↓
M9 Workflow hardening + CI/CD + Phase 1 freeze
```

---

## 4. Repository Structure

The repository is organized into distinct directories to support reproducible execution and auditing:

```text
Pilot_tumor/
├── PROJECT.md          # Authoritative project specification
├── PROGRESS.md         # Milestone execution status tracking
├── CHANGELOG.md        # Technical changes & scientific decisions log
├── README.md           # This project overview and documentation
├── .gitignore          # File exclusions for Git tracking
├── config/             # YAML configurations and JSON schemas
│   ├── config.yaml     # Production real-data configurations
│   ├── config.test.yaml# Synthetic test configuration
│   └── schemas/        # Schema files validating config structures
├── workflow/           # Snakemake orchestrations
│   ├── Snakefile       # Core Snakemake execution file
│   └── envs/           # Conda environment specifications
├── scripts/            # Modular R, Python, and Shell helper scripts
│   ├── R/              # Core Seurat/R data processing scripts
│   ├── python/         # Python summary report generation scripts
│   └── shell/          # HPC execution and milestone verification scripts
├── tests/              # Verification unit tests
│   └── unit/           # Script-specific unit tests run on synthetic data
├── reports/            # Markdown reports, TSV recommendations, & diagnostic figures
│   ├── milestones/     # Consolidated reports for M0-M6
│   ├── datasets/       # Dataset-specific QC and PCA reports/figures
│   ├── qc_optimization/# QC threshold optimization sensitivity plots
│   └── qc_comparison/  # Comparative analysis of QC strategies
├── results/            # Computed Seurat RDS objects (Git-ignored)
└── logs/               # Run logs and SLURM outputs (Git-ignored)
```

---

## 5. Reproducibility Framework

This repository enforces strict reproducibility across all milestones via:
1. **Snakemake Orchestration**: The workflow is fully managed by Snakemake, defining clear rule-level dependencies, inputs, outputs, and compute resource limits.
2. **Conda Environments**: R/Python package dependencies are locked using environment specifications.
3. **Synthetic Testing**: Synthetic data is used to validate workflow changes without executing expensive real-data jobs.
4. **SLURM Integration**: Jobs are configured with specific resource envelopes (CPUs, Memory, Runtime) and executed via SLURM on HPC clusters.
5. **Provenance Tracking**: Execution runs generate JSON provenance logs recording file checksums, R session metadata, and runtime parameters.
6. **Milestone Tracking**: Comprehensive milestone logs in `PROGRESS.md` and `CHANGELOG.md` track the historical evolution of scientific and technical decisions.

---

## 6. Environment Setup

Depending on the environment, the R and Snakemake runtime can be initialized as follows:

### Developer/HPC Environment (Pre-configured)
In the primary HPC system, a pre-compiled environment is loaded using the following commands:
```bash
# Source Conda profile
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh

# Activate the pre-configured environment
conda activate R_env
```

### Portable User Installation (New Environments)
For users running this workflow on a new system or node, the environment can be constructed using the provided portable YAML configuration:
```bash
# Create the conda environment from the portable specification
conda env create -f workflow/envs/R_env_portable.yaml -n R_env

# Activate the new environment
conda activate R_env
```

---

## 7. Running the Workflow

> [!CAUTION]
> **Real-Data Compute Safety Constraint**:
> Never load, deserialize, inspect, subset, or analyze the real raw/processed data (`processed_mpnst.rds` or objects in `results/`) on a login node.
> All production calculations must run inside a SLURM job allocation.

### Environment & Safety Verification
Run the verification script to confirm repository sanity, clean whitespaces, check paths, and ensure no direct login node RDS reads:
```bash
./scripts/shell/verify_milestone.sh
```

### Synthetic/Test Configuration
Verify the Snakemake setup and execute a smoke test run using the synthetic dataset:
```bash
# 1. Perform a dry-run to verify rule dependencies
snakemake -n --configfile config/config.test.yaml

# 2. Execute the synthetic workflow locally
snakemake --cores 4 --configfile config/config.test.yaml
```

### Real-Data Workflow Execution
To run the production workflow on real datasets:
```bash
# 1. Perform a production dry-run
snakemake -n --configfile config/config.yaml

# 2. Submit the workflow job to the SLURM partition 'ihc'
sbatch scripts/shell/run_m6_workflow.sh
```

---

## 8. Current Outputs (Through Milestone 6)

Successful execution of Milestones M0–M6 yields the following major artifacts:

- **Seurat RDS Objects** (stored in `results/datasets/{ds}/`):
  - `{ds}_raw.rds`: Raw extracted datasets split by `sample_id`.
  - `{ds}_filtered_specific.rds`: Filtered single-cell objects after scDblFinder doublet removal and dataset-specific QC thresholds.
  - `{ds}_normalized.rds`: SCTransform-normalized and variance-stabilized Seurat objects with 3,000 highly variable features.
  - `{ds}_pca.rds`: PCA-embedded Seurat objects computed on variable features.
  - `{ds}_clustered.rds`: Clustered Seurat objects containing sweep resolution metadata and active identity set to the recommended resolution.
- **Recommendations & Indexes** (stored in `reports/`):
  - [reports/PCA_RECOMMENDATIONS.tsv](file://reports/PCA_RECOMMENDATIONS.tsv): Machine-readable Recommended, Conservative, and Maximum PC counts for downstream clustering.
  - [reports/CLUSTERING_RECOMMENDATIONS.tsv](file://reports/CLUSTERING_RECOMMENDATIONS.tsv): Machine-readable selected resolutions and alternative recommendations for each dataset.
  - [reports/CLUSTERING_SWEEP_SUMMARY.tsv](file://reports/CLUSTERING_SWEEP_SUMMARY.tsv): Comprehensive metrics summary (cluster sizes, stability ARI, technical covariate correlation) across all resolution sweep values.
  - [reports/FIGURE_INDEX.tsv](file://reports/FIGURE_INDEX.tsv): Consolidated index of all 60 diagnostic plots generated during QC, normalization, PCA, and clustering.
- **Consolidated Milestone Reports** (stored in `reports/milestones/`):
  - [M0_REPORT.md](reports/milestones/M0_REPORT.md): Infrastructure, Environment, & SLURM Safety.
  - [M1_REPORT.md](reports/milestones/M1_REPORT.md): Real-Data Object & Layers Inventory.
  - [M2_REPORT.md](reports/milestones/M2_REPORT.md): Extraction & Pre-Filter metrics.
  - [M3_REPORT.md](reports/milestones/M3_REPORT.md): QC Filtering & doublet validation.
  - [M4_REPORT.md](reports/milestones/M4_REPORT.md): Normalization & HVF selection.
  - [M5_REPORT.md](reports/milestones/M5_REPORT.md): Principal Component Analysis & evaluation.
  - [M6_REPORT.md](reports/milestones/M6_REPORT.md): Clustering Resolution Sweep & Selection.

### Summary of Milestone 6 Computational Recommendations
> [!NOTE]
> These resolutions are **computational recommendations** based on stability metrics and covariate correlation heuristics. The final researcher-approved resolutions have **not** yet been selected.

- **MPNST_1**: Recommended resolution **0.6** resolving **18** clusters (Bootstrap Stability ARI: `0.919`, no technical concern).
- **MPNST_2**: Recommended resolution **0.3** resolving **9** clusters (Bootstrap Stability ARI: `0.933`, no technical concern).
- **MPNST_3**: Recommended resolution **0.6** resolving **13** clusters (Bootstrap Stability ARI: `0.910`, no technical concern).
- **MPNST_4**: Recommended resolution **0.7** resolving **14** clusters (Bootstrap Stability ARI: `0.740`, technical concern: `MT_Bias` with percent.mt $R^2 = 0.47$).

---

## 9. Reference Documentation
For detailed progress, requirements, and historical records:
- [PROJECT.md](PROJECT.md) - Authoritative project specification.
- [PROGRESS.md](PROGRESS.md) - Project milestone history and log.
- [CHANGELOG.md](CHANGELOG.md) - Technical changes & decisions log.
- [reports/milestones/](reports/milestones/) - Directory containing all milestone reports.
