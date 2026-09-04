# MPNST Phase 1 Analysis Workflow
[![Phase 1 CI Workflow](https://github.com/vipinmenon1989/Pliot_Tumor_Research/actions/workflows/ci.yml/badge.svg)](https://github.com/vipinmenon1989/Pliot_Tumor_Research/actions/workflows/ci.yml)

Reproducible, configuration-driven single-cell RNA-seq preprocessing and dataset auditing workflow for MPNST datasets.

---

## 1. Project Purpose & Scientific Objective
The primary objective of this project is to construct a modular, reproducible, configuration-driven, and HPC-ready Phase 1 workflow to preprocess and analyze constituent MPNST (malignant peripheral nerve sheath tumor) datasets independently. Crucially, all datasets are processed and evaluated in isolation before any batch correction or multi-dataset integration is performed. This prevents technical confounders from distorting early-stage QC and dimensionality reduction.

---

## 2. Current Project Status
- **Current Phase**: **Phase 6 — complete and frozen** (2026-09-04). Phases 1–5 are complete and
  frozen and were **not modified**; Phase 5 read the frozen Phase 4 object read-only, Phase 6 read
  the frozen Phase 5 object, and all upstream checksums were re-verified after every milestone.
  **Phase 7 is NOT initiated and requires separate authorization.**
- **Phase 5 headline — a negative result, and a clean one.** Phase 4 found 0 of 8 *discrete*
  malignant states recurrent across patients. Phase 5 asked the *continuous* version on the same
  6,434 malignant cells with the same recurrence rule shape, using consensus matrix factorization
  (cNMF, K = 8 by a rule declared before any factor was seen). **The answer is also 0: 0 recurrent,
  1 shared-limited, 7 patient-private — and the one shared-limited program is technical
  (ribosomal), so not a single biologically interpretable program is carried by even two patients.**
  This survives the most aggressive check available: a full re-run with ribosomal, pseudogene and
  ambient-myeloid genes removed across the whole K grid still gives **0 recurrent programs at every
  K from 5 to 15**. The one program that comes close is the **ECM/mesenchymal program P3**, which
  reaches three patients only when the activity threshold is relaxed from 20% to 10%.
  The five patient-private `Mesenchymal_ECM-*` states are **not** one shared program — they map to
  three, though MPNST_4's three ECM states do collapse onto one. Phase 5 also derived a transparent
  **92-gene signature separating malignant ECM-like MPNST cells from genuine fibroblasts**
  (evaluable in 2 of 4 patients); the fibroblast side is marked by **CDH19, APOD, SCN7A and
  ABCA6/8/9/10** — nerve-associated / endoneurial stroma. Handoff
  [reports/phase5/PHASE5_HANDOFF.md](reports/phase5/PHASE5_HANDOFF.md).
- **Phase 6 headline — transcriptional phenotype is largely decoupled from clone structure.**
  Combining Phase 5 programs with the frozen Phase 4 SCEVAN clones (MPNST_1/2/4 only; **MPNST_3's
  clone structure is excluded**, it failed the Phase 4 immune sanity gate): **median η² = 0.059, so
  roughly 94% of each program's variance sits WITHIN clones, and between-clone divergence is only
  0.7% of within-clone dispersion** — a patient's CNA-defined clones are transcriptionally
  near-interchangeable. **Model C (mixed architecture) was selected** against pre-specified numeric
  criteria. Within-clone program diversity is real but patient-specific (MPNST_1 clones span
  1.69–2.83 effective programs; MPNST_2's span exactly 1.00). Regulatory layer: **E2F4/E2F1/MYC for
  Cycling and HIF1A/ATF4/HSF1 for the translation-stress program, recovered without being imposed**,
  with PROGENy and Hallmark independently supporting labels derived from the programs' own genes.
  Handoff [reports/phase6/PHASE6_HANDOFF.md](reports/phase6/PHASE6_HANDOFF.md).
- **What Phase 5/6 cannot do.** n = 4 patients and `sample_id` = patient = dataset, so a
  patient-private *biological* program cannot be separated from patient-level *technical* structure;
  clone conclusions rest on 3 patients; SCEVAN infers copy number from expression, so the
  CNA→expression analysis is an internal consistency check, **not** validation; and within-clone
  program diversity is consistent with plasticity but does **not** demonstrate a state transition.
- **Phase 4 — complete and frozen** (2026-09-03)
- **Phase 1**: complete and frozen (M0–M9) · **Phase 2**: complete and frozen (M10–M17 + M15A)
  · **Phase 3**: complete and frozen (M18–M27) · **Phase 4**: complete and frozen (M28–M35,
  plus **M35A** — SCEVAN figure consolidation, visualization only)
  · **Phase 5**: complete and frozen (M36–M41) · **Phase 6**: complete and frozen (M42–M50)
- **Phase 4 headline**: the malignant compartment is **1.88× larger than Phase 2 estimated**
  (3,420 → **6,434 cells**, 17.35% → **32.63%**), and **4,036 of 5,064 cells Phase 2 called
  `Fibroblast` carry inferred copy-number alterations** — majority-malignant in 3 of 4 patients.
  The Phase 3 CCC architecture **survives** the relabelling (high-concordance interactions
  −2.6%, median change across 23 named axes −3.1%, no axis lost); what changed is
  **attribution**, with 3,051 sender reassignments concentrated in the ECM→integrin axis.
  Final object `results/phase4/phase4_final_object.rds` (md5 `e85ba848…`, 23/23 validation
  checks). Handoff [reports/phase4/PHASE4_HANDOFF.md](reports/phase4/PHASE4_HANDOFF.md).
- **Phase 4 caveats that travel with those numbers**: the 32.63% figure is **not
  patient-robust** (dropping MPNST_4 returns it to the Phase 2 value); **no malignant
  transcriptional state is recurrent** across patients (0 of 8), so a tumour-state → TME model
  could not be built; and one sanity-gate amendment was made **after** seeing which sample it
  excluded — documented as such in
  [reports/phase4/MALIGNANCY_DECISION_RULES.md](reports/phase4/MALIGNANCY_DECISION_RULES.md).
- **M35A — the SCEVAN evidence, made visible** (2026-09-03, **visualization only; nothing was
  re-run and no call, rule, clone, state or threshold changed; the final object was neither loaded
  nor modified**). The native SCEVAN CNA heatmaps are surfaced **verbatim** into
  `results/phase4/figures/final/` as `17`–`20`, each beside a companion panel built from the same
  frozen clone assignments, and eight custom figures (`21`–`28`) carry the argument:
  **clone composition is the strongest evidence** — SCEVAN builds subclones from copy number
  alone and never sees a Phase 2 label, and **4 of 4 MPNST_2 clones and 7 of 8 MPNST_4 clones are
  fibroblast-dominated** while all 7 MPNST_1 clones mix the disputed identities. Malignant-called
  fibroblast-labelled cells separate from non-malignant ones with **Cliff's delta 0.78–0.92** and
  their genome-wide CNA profile correlates **0.934 / 0.972** with the malignant compartment against
  **0.184** with non-malignant fibroblasts — **in the 2 of 4 patients where that contrast can be
  made at all**, which the figures state rather than hide. **Fibroblast → Malignant is 4,036 at
  every tested threshold**; the 32.63% cohort fraction is not. MPNST_3 is shown failing, not
  omitted. Audit
  [reports/phase4/M35A_SCEVAN_FIGURE_AUDIT.md](reports/phase4/M35A_SCEVAN_FIGURE_AUDIT.md);
  handoff §28. **SCEVAN infers copy number from expression and provides no DNA-level proof.**
- **Phase 3 headline**: 36,486 supported ligand–receptor interactions across LIANA, CellChat and
  CellPhoneDB; **547 tumour-centric interactions with ≥3 independent evidence streams**. Primary
  deliverable `results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv`.

### Phase 4 — what it is for

Phase 2 defined `MPNST-Tumor` **conservatively**: 3,420 cells (17.35%), only where positive
MPNST/Schwann/neural-crest marker evidence existed. That was the right call on marker evidence
alone, but it leaves 5,064 `Fibroblast`, 1,231 `Candidate-Malignant-Unresolved`, 720
`Uncertain` and 435 `Pericyte-VSMC` cells whose malignant status marker expression cannot
settle — and that uncertainty propagates into every tumour-centric Phase 3 interaction.

**SCEVAN** (De Falco *et al.*, Nat Commun 14:1074, 2023; doi:10.1038/s41467-023-36790-9)
supplies the orthogonal evidence: inferred large-scale copy number from raw counts, with
automatic confident-normal detection and subclonal resolution. Phase 4 uses it to refine
malignant identity, resolve tumour states and clones, and then run a **targeted sensitivity
analysis** on the Phase 3 conclusions.

**Phase 4 does not invalidate Phase 2 or Phase 3.** Phase 2 annotation is historical
biological annotation; Phase 4 malignancy is additional orthogonal evidence. No Phase 2 or
Phase 3 label is overwritten — Phase 4 only adds metadata fields. See `PROJECT.md` §P4.

- **Phase 4 directories**: scripts `scripts/phase4/{scevan,malignancy,tumor_states,ccc_refinement,utils}`
  · SLURM `scripts/shell/phase4/` · results `results/phase4/` · reports `reports/phase4/`
  · logs `logs/phase4/`
- **Phase 4 environment**: `R_env`, **unchanged** — SCEVAN 1.0.3 and yaGST 2017.8.25 were
  already installed, so no dependency moved. Phase 3's CCC stack is reused at identical
  versions so the Phase 3 → Phase 4 comparison is not confounded by tool drift.
- **Phase 4 directories**: results `results/phase4/` (figures `figures/final/`, tables
  `tables/final/`, manifest `phase4_manifest.json`) · reports `reports/phase4/` · scripts
  `scripts/phase4/` · SLURM `scripts/shell/phase4/` · logs `logs/phase4/`
- **Phase 4 environment**: `R_env` **unchanged** — SCEVAN 1.0.3 and yaGST were already
  installed, so **no dependency moved**. One addition anywhere: Python umap-learn 0.5.12 in an
  isolated `p4_umap_env`, required by SCEVAN's subclone stage.
- **Phase 5**: **not begun; requires separate authorization.**
- **Status**: **Completed: Milestones M0–M9** (Phase 1) and **M10–M17** (Phase 2)
- **Latest Completed Milestone**: **M17 — Phase 2 Validation, Freeze and Handoff** (re-frozen 2026-09-03 after the M15A CCC-annotation amendment)
- **Phase Status**: **PHASE 1 = FROZEN** | **PHASE 2 = COMPLETE**
- **Final Phase 2 object**: `results/phase2/phase2_final_object.rds`
  (19,716 cells; md5 `153d5f6acc70f9c05aa48cabc4f4ac2d`) — carries **both** the detailed
  literature-supported annotation and the CCC-oriented `annotation_ccc` layer
  (`MPNST-Tumor` = 3,420 cells, 17.35%)
- **Phase 2 handoff**: [reports/phase2/PHASE2_HANDOFF.md](reports/phase2/PHASE2_HANDOFF.md)
  · CCC readiness [reports/phase2/CCC_READINESS.md](reports/phase2/CCC_READINESS.md)
  · manifest `results/phase2/phase2_manifest.json`
- **Phase 3 architecture**: sample-aware CCC → multi-method inference
  (LIANA · CellChat · CellPhoneDB) → concordance → tumour↔immune prioritisation →
  NicheNet receiver-response → LochNESS receiver-state → integrated interaction map.
  See `PROJECT.md` §P3 and [reports/phase3/PHASE3_METHOD_PLAN.md](reports/phase3/PHASE3_METHOD_PLAN.md).
- **Phase 3 environment**: `R_env` (unchanged Phase 2 stack) plus liana/CellChat/nichenetr;
  CellPhoneDB in an isolated `cpdb_env`. See
  [reports/phase3/environment/PHASE3_INSTALL_LOG.md](reports/phase3/environment/PHASE3_INSTALL_LOG.md).
- **Where things live**: final figures `results/phase3/figures/final/` · final tables
  `results/phase3/tables/` and `results/phase3/ccc/prioritized/` · milestone reports
  `reports/phase3/milestones/` · Phase 3 manifest `results/phase3/phase3_manifest.json`.

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
M8 Combined pre-integration baseline (SCT, shared PCA, shared UMAP, neighbor mixing metrics) ✓
    ↓
M9 Workflow hardening + CI/CD + Phase 1 freeze (reconciled and tested) ✓
    ↓
=== PHASE 1 FROZEN / PHASE 2 NOT STARTED ===
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
Verify system packages, configuration compatibility, directory permissions, and dataset visibility before running execution:
```bash
# Run preflight checker in synthetic mode (for testing/CI)
python scripts/python/preflight_checker.py --mode synthetic

# Run preflight checker in real/HPC mode (for real cluster execution)
python scripts/python/preflight_checker.py --mode real
```

### Synthetic/Test Configuration (Quick-Start)
Verify the Snakemake setup and execute a smoke test run using the synthetic dataset:
```bash
# 1. Run the preflight checker in synthetic mode
python scripts/python/preflight_checker.py --mode synthetic

# 2. Perform a dry-run to verify rule dependencies
snakemake -n --configfile config/config.test.yaml

# 3. Execute the synthetic workflow locally
snakemake --cores 4 --configfile config/config.test.yaml
```

### Real-Data/HPC Workflow Execution
To run the production workflow on the real datasets:
```bash
# 1. Run the preflight checker in real/HPC mode
python scripts/python/preflight_checker.py --mode real

# 2. Perform a production dry-run
snakemake -n --configfile config/config.yaml
```

---

## 8. Current Outputs (Through Milestone 9)

### Seurat RDS & Manifest Objects (stored in `results/datasets/{ds}/` or `results/`)
- `{ds}_raw.rds` (M2): Raw extracted datasets split by `sample_id`.
- `{ds}_filtered_specific.rds` (M3): Filtered single-cell objects after scDblFinder doublet removal and dataset-specific QC thresholds.
- `{ds}_normalized.rds` (M4): SCTransform-normalized and variance-stabilized Seurat objects with 3,000 highly variable features.
- `{ds}_pca.rds` (M5): PCA-embedded Seurat objects computed on variable features.
- `{ds}_clustered.rds` (M6): Clustered Seurat objects containing sweep resolution metadata (0.1 to 1.0) and active identity set to the recommended resolution.
- `results/combined/pre_integration/combined_preintegration.rds` (M8): Combined pre-integration Seurat object containing all 19,716 cells in a shared non-integrated expression space.
- `results/phase1_manifest.json` (M9): Programmatic, machine-readable validation manifest contract for Phase 2.

### Recommendations & Reports (stored in `reports/`)
- [reports/PCA_RECOMMENDATIONS.tsv](file://reports/PCA_RECOMMENDATIONS.tsv) (M5): Machine-readable Recommended, Conservative, and Maximum PC counts.
- [reports/CLUSTERING_RECOMMENDATIONS.tsv](file://reports/CLUSTERING_RECOMMENDATIONS.tsv) (M6): Machine-readable selected resolutions and alternative recommendations for each dataset.
- [reports/CLUSTERING_SWEEP_SUMMARY.tsv](file://reports/CLUSTERING_SWEEP_SUMMARY.tsv) (M6): Summary of metrics (stability ARI, cluster sizes, covariate correlation) across sweep values.
- [reports/datasets/{ds}/ANALYSIS_RECOMMENDATION.md](file://reports/datasets/) (M7): Dataset-specific reports detailing the scientific rationale for resolution selection.
- [reports/datasets/{ds}/markers/resolution_{res}/](file://reports/datasets/) (M7): Marker tables for all 40 combinations.
- [reports/PRE_INTEGRATION_ASSESSMENT.md](file://reports/PRE_INTEGRATION_ASSESSMENT.md) (M8): Pre-integration baseline assessment report.
- [reports/INTEGRATION_PREPARATION.md](file://reports/INTEGRATION_PREPARATION.md) (M8): Integration preparation design report.
- [reports/combined/pre_integration/](file://reports/combined/pre_integration/) (M8): Pre-integration inventories, dictionaries, composition summaries, confounding analyses, and neighborhood mixing diagnostic tables.
- [reports/audits/M0_M8_RECONCILIATION.md](file://reports/audits/M0_M8_RECONCILIATION.md) (M9): Three-layer verification audit report.
- [reports/audits/M0_M8_RECONCILIATION.tsv](file://reports/audits/M0_M8_RECONCILIATION.tsv) (M9): Machine-readable audit TSV.
- [reports/PHASE1_HANDOFF.md](file://reports/PHASE1_HANDOFF.md) (M9): Human-readable final Phase 1 handoff documentation.
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
- [M8_REPORT.md](reports/milestones/M8_REPORT.md): Combined Pre-Integration Baseline.
- [M9_REPORT.md](reports/milestones/M9_REPORT.md): Workflow Hardening, CI/CD, Provenance, and Phase 1 Freeze.

---

## 9. Current Scientific Handoff

As of the completion of Milestone 9, the Phase 1 pre-integration baseline is frozen with the following dataset recommendations:

| Dataset | Recommended PCs | Recommended Resolution | Alternative Resolution | Status | Resolved Clusters | Covariate Concerns / Limitations |
| --- | :---: | :---: | :---: | :---: | :---: | :--- |
| **MPNST_1** | 1–8 | **0.6** | 0.3 | CONFIRMED | 18 | None |
| **MPNST_2** | 1–6 | **0.3** | 0.5 | CONFIRMED | 9 | None |
| **MPNST_3** | 1–9 | **0.6** | 0.3 | CONFIRMED | 13 | None |
| **MPNST_4** | 1–5 | **0.7** | 0.5 | CONFIRMED | 14 | Technical MT correlation ($R^2 = 0.47$ at res 0.7); alternative res 0.5 recommended to reduce MT correlation ($R^2 = 0.24$) |

### Important Scientific Context
- **Pre-Integration Freeze**: The datasets are combined but remain non-integrated. **No integration or batch correction has been executed** (e.g. no Harmony, CCA, RPCA, MNN, scVI).
- **Multi-Resolution Preserved**: All clustering resolutions (0.1 to 1.0) and SNN graphs are preserved in the Seurat objects, with marker tables pre-computed for every resolution.
- **Pre-Annotation Status**: No final biological cell-type annotations have been assigned to clusters; clusters are currently defined by their numerical partitions and associated marker signatures.
- **Downstream Transition**: The workflow has completed Milestone M9 (Workflow Hardening, CI/CD, and Phase 1 Freeze). **Phase 1 is FROZEN. Phase 2 is NOT STARTED.** No integration or batch correction has yet been performed.

---

## 10. Reference Documentation
For detailed progress, requirements, and historical records:
- [PROJECT.md](PROJECT.md) - Authoritative project specification.
- [PROGRESS.md](PROGRESS.md) - Project milestone history and log.
- [CHANGELOG.md](CHANGELOG.md) - Technical changes & decisions log.
- [reports/milestones/](reports/milestones/) - Directory containing all milestone reports.
- [reports/PHASE1_HANDOFF.md](reports/PHASE1_HANDOFF.md) - Final Phase 1 handoff documentation.
- [results/phase1_manifest.json](results/phase1_manifest.json) - Final Phase 1 manifest contract.

