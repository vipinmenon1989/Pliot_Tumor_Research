# MPNST Phase 1 Analysis Workflow

Reproducible, configuration-driven single-cell RNA-seq preprocessing and dataset auditing workflow for MPNST datasets.

---

## 1. Project Purpose & Scientific Objective
The primary objective of this project is to construct a modular, reproducible, configuration-driven, and HPC-ready Phase 1 workflow to preprocess and analyze constituent MPNST (malignant peripheral nerve sheath tumor) datasets independently. Crucially, all datasets are processed and evaluated in isolation before any batch correction or multi-dataset integration is performed. This prevents technical confounders from distorting early-stage QC and dimensionality reduction.

---

## 2. Current Project Status
- **Current Phase**: Phase 1 (Independent Dataset Processing)
- **Status**: **Completed: Milestones M0–M7**
- **Latest Completed Milestone**: **M7 — Marker Discovery and Dataset-Specific Recommendations**
- **Next Required Step**: **Post-M7 Independent Audit Gate**
- **Next Computational Milestone (After Approval)**: **Milestone M8 — Combined Pre-Integration Baseline**

> [!IMPORTANT]
> **Strict Phase 1 Scope Constraint**:
> This workflow strictly limits processing to independent dataset analysis. Downstream operations such as integration methods (Harmony, Seurat CCA/RPCA, FastMNN, scVI, BBKNN, etc.), automated cell-type annotation, condition-level differential expression, pathway analysis, trajectory analysis, cell-cell communication, and predictive modeling are strictly prohibited in Phase 1.

---

## 3. Workflow Overview

```text
M0 Infrastructure (Environment, synthetic data generation, SLURM verification) ✓
    ↓
M1 Real-data audit (Authoritative object structure and inventory checks) ✓
    ↓
M2 Dataset extraction + pre-filter QC (Deterministic extraction & diagnostic plotting) ✓
    ↓
M3 QC filtering + doublet assessment (Doublet detection & dataset-specific thresholds) ✓
    ↓
M4 Normalization + variable features (SCTransform v2 variance stabilization & HVF selection) ✓
    ↓
M5 PCA + PC evaluation (Independent PCA & geometric elbow selection) ✓
    ↓
M6 Clustering resolution sweep (SNN graph construction & sweeps 0.1–1.0) ✓
    ↓
M7 Marker discovery + recommendations (Resolution-specific markers & selection) ✓
    ↓
POST-M7 INDEPENDENT AUDIT GATE  <-- CURRENT GATE
    ↓
M8 Combined pre-integration baseline (NOT STARTED)
    ↓
M9 Workflow hardening + CI/CD + Phase 1 freeze (NOT STARTED)
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
│   ├── milestones/     # Consolidated reports for M0-M7
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
To run the production workflow on the real datasets:
```bash
# 1. Perform a production dry-run
snakemake -n --configfile config/config.yaml

# 2. Execute production steps via SLURM
# Note: A generic end-to-end production entry point is not yet implemented.
# Production execution currently uses milestone-specific SLURM wrapper scripts
# to manage computing resources safely. Complete workflow generalization is scheduled for M9.
sbatch scripts/shell/run_m7_workflow.sh
```

---

## 8. Current Outputs (Through Milestone 7)

Successful execution of Milestones M0–M7 yields the following major artifacts:

### Seurat RDS Objects (stored in `results/datasets/{ds}/`)
- `{ds}_raw.rds` (M2): Raw extracted datasets split by `sample_id`.
- `{ds}_filtered_specific.rds` (M3): Filtered single-cell objects after scDblFinder doublet removal and dataset-specific QC thresholds.
- `{ds}_normalized.rds` (M4): SCTransform-normalized and variance-stabilized Seurat objects with 3,000 highly variable features.
- `{ds}_pca.rds` (M5): PCA-embedded Seurat objects computed on variable features.
- `{ds}_clustered.rds` (M6): Clustered Seurat objects containing sweep resolution metadata (0.1 to 1.0) and active identity set to the recommended resolution.

### Recommendations & Reports (stored in `reports/`)
- [reports/PCA_RECOMMENDATIONS.tsv](file://reports/PCA_RECOMMENDATIONS.tsv) (M5): Machine-readable Recommended, Conservative, and Maximum PC counts.
- [reports/CLUSTERING_RECOMMENDATIONS.tsv](file://reports/CLUSTERING_RECOMMENDATIONS.tsv) (M6): Machine-readable selected resolutions and alternative recommendations for each dataset.
- [reports/CLUSTERING_SWEEP_SUMMARY.tsv](file://reports/CLUSTERING_SWEEP_SUMMARY.tsv) (M6): Summary of metrics (stability ARI, cluster sizes, covariate correlation) across sweep values.
- [reports/datasets/{ds}/ANALYSIS_RECOMMENDATION.md](file://reports/datasets/) (M7): Dataset-specific reports detailing the scientific rationale for resolution selection.
- [reports/datasets/{ds}/markers/resolution_{res}/](file://reports/datasets/) (M7): Marker tables for all 40 combinations:
  - `markers_all.tsv`: All detected markers.
  - `markers_filtered.tsv`: Markers filtered by significance (adj. p-value < 0.05) and log2 fold change (> 0.25).
  - `top_markers.tsv`: Top 20 ranked markers per cluster.
  - `marker_summary.tsv`: Global summary of cluster stats.
- [reports/datasets/{ds}/markers/resolution_{res}/figures/](file://reports/datasets/) (M7): Selected top-marker visualizations including downsampled (max 100 cells/cluster) heatmaps, dot plots, and FeaturePlots in PDF and PNG.
- [reports/FIGURE_INDEX.tsv](file://reports/FIGURE_INDEX.tsv): Consolidated index of all diagnostic plots generated across all milestones.

### Consolidated Milestone Reports (stored in `reports/milestones/`)
- [M0_REPORT.md](reports/milestones/M0_REPORT.md): Infrastructure, Environment, & SLURM Safety.
- [M1_REPORT.md](reports/milestones/M1_REPORT.md): Real-Data Object & Layers Inventory.
- [M2_REPORT.md](reports/milestones/M2_REPORT.md): Extraction & Pre-Filter metrics.
- [M3_REPORT.md](reports/milestones/M3_REPORT.md): QC Filtering & doublet validation.
- [M4_REPORT.md](reports/milestones/M4_REPORT.md): Normalization & HVF selection.
- [M5_REPORT.md](reports/milestones/M5_REPORT.md): Principal Component Analysis & evaluation.
- [M6_REPORT.md](reports/milestones/M6_REPORT.md): Clustering Resolution Sweep & Selection.
- [M7_REPORT.md](reports/milestones/M7_REPORT.md): Marker Discovery & Dataset Recommendations.

---

## 9. Current Scientific Handoff

As of the completion of Milestone 7, the Phase 1 pre-integration baseline is frozen with the following dataset recommendations:

| Dataset | Recommended PCs | Recommended Resolution | Alternative Resolution | Status | Resolved Clusters | Covariate Concerns / Limitations |
| --- | :---: | :---: | :---: | :---: | :---: | :--- |
| **MPNST_1** | 1–8 | **0.6** | 0.3 | CONFIRMED | 18 | None |
| **MPNST_2** | 1–6 | **0.3** | 0.5 | CONFIRMED | 9 | None |
| **MPNST_3** | 1–9 | **0.6** | 0.3 | CONFIRMED | 13 | None |
| **MPNST_4** | 1–5 | **0.7** | 0.5 | CONFIRMED | 14 | Technical `MT_Bias` ($R^2 = 0.47$ at res 0.7); alternative res 0.5 recommended to reduce MT correlation ($R^2 = 0.24$) |

### Important Scientific Context
- **Pre-Integration Freeze**: The datasets are processed and stored independently. **No integration or batch correction has been executed** (e.g. no Harmony, CCA, RPCA, MNN, scVI).
- **Multi-Resolution Preserved**: All clustering resolutions (0.1 to 1.0) and SNN graphs are preserved in the Seurat objects, with marker tables pre-computed for every resolution.
- **Pre-Annotation Status**: No final biological cell-type annotations have been assigned to clusters; clusters are currently defined by their numerical partitions and associated marker signatures.
- **Downstream Transition**: The workflow is currently halted at the **Post-M7 Independent Audit Gate**. Milestone M8 (Combined Pre-Integration Baseline) has **NOT** started.

---

## 10. Reference Documentation
For detailed progress, requirements, and historical records:
- [PROJECT.md](PROJECT.md) - Authoritative project specification.
- [PROGRESS.md](PROGRESS.md) - Project milestone history and log.
- [CHANGELOG.md](CHANGELOG.md) - Technical changes & decisions log.
- [reports/milestones/](reports/milestones/) - Directory containing all milestone reports.
