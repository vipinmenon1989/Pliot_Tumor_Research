# MPNST Phase 1 Analysis Workflow

Reproducible, configuration-driven single-cell RNA-seq preprocessing and dataset auditing workflow for MPNST datasets.

## 1. Project Purpose & Scientific Objective
The primary objective of this project is to build a modular, reproducible, configuration-driven, and HPC-ready Phase 1 workflow to preprocess and analyze constituent datasets independently before any batch correction or integration is performed.

## 2. Phase 1 Scope & Strict Prohibitions
Phase 1 covers:
- Dataset auditing, quality control (QC) filtering, doublet assessment, normalization, variable feature selection, PCA, clustering sweeps (0.1–1.0), marker gene discovery, and pre-integration diagnostics.
- **Strict Prohibitions**: Integration methods (Harmony, Seurat CCA/RPCA, FastMNN, scVI, BBKNN, etc.), automated cell-type annotation, condition-level differential expression, pathway analysis, trajectory analysis, cell-cell communication, and predictive modeling are strictly prohibited in Phase 1.

## 3. Repository Structure
```text
Pilot_tumor/
├── PROJECT.md          # Authoritative project specification
├── PROGRESS.md         # Milestone execution status tracking
├── README.md           # This project overview and documentation
├── .gitignore          # File exclusions for Git tracking
├── config/             # YAML configurations and JSON schemas
├── workflow/           # Snakemake workflows, environments, profiles
├── scripts/            # Modular R, Python, and Shell helper scripts
├── tests/              # Smoke and unit tests
├── data/               # Input data (ignored by Git, except synthetic)
├── results/            # Outputs (ignored by Git)
├── reports/            # Markdown reports and audits
└── logs/               # Log output directories
```

## 4. Environment & Execution

### Conda Environment
- Required environment: `R_env`
- Sourced from: `/local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh`
- Active env location: `/local/projects-t3/lilab/vmenon/anaconda3/envs/R_env`
- Snakemake version: `9.23.1` (installed directly into `R_env` without modifying critical packages)

### Activating the Environment
```bash
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
```

### Running the Synthetic Workflow Locally
```bash
snakemake --cores 1 --configfile config/config.test.yaml
```

### Running a Snakemake Dry-Run
```bash
snakemake -n --configfile config/config.test.yaml
```

---

## 5. Real-Data Compute Safety Rules

> [!IMPORTANT]
> **Strict Safety Constraint**: Never load, deserialize, inspect, subset, or analyze `processed_mpnst.rds` directly on the login node.
> Any command loading the real RDS dataset must run inside a SLURM allocation. Real-data computation must be submitted via `sbatch`.

### SLURM Configuration
- **Account**: `ihc`
- **Partition**: `ihc`
- **Preferred Node**: `ihc-grid-1-1-1`
- **Maximum Resource Envelope**: 32 CPUs, 450GB RAM, 72 hours walltime

---

## 6. Milestone Model & Roles
- **Execution Model**: Incremental development with mandatory STOP gates at each milestone boundary. Progression requires explicit researcher approval.
- **Roles**:
  - **Antigravity**: Primary implementation agent.
  - **Cursor**: Independent audit agent.

### Current Project Status
- **Milestone 0**: Completed / Frozen.
- **Milestone 1**: Completed / Frozen.
- **Milestone 2**: Completed (Awaiting researcher approval).

---

## 7. Reference Documentation
- [PROJECT.md](PROJECT.md) - Authoritative project specification.
- [PROGRESS.md](PROGRESS.md) - Project milestone history and log.
- [M0_REPORT.md](reports/milestones/M0_REPORT.md) - Detailed Milestone 0 report.
- [M1_REPORT.md](reports/milestones/M1_REPORT.md) - Detailed Milestone 1 report.
- [M2_REPORT.md](reports/milestones/M2_REPORT.md) - Detailed Milestone 2 report.

---

## 8. Reproducibility
- **Environment Snapshots**: Pre- and post-installation environments are stored in `workflow/envs/`.
- **Portable Specification**: [R_env_portable.yaml](workflow/envs/R_env_portable.yaml) represents the machine-agnostic environment configuration.
- **Synthetic Tests**: Exercises production workflows using a clean generated Seurat dataset.
- **Logging & Provenance**: Standardized structured logs and MD5 checksum tracking for data lineage.
- **Milestone History**: Git commits track clean boundaries for audit validation.
