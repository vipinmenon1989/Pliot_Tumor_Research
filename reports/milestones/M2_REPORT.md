# Milestone 2 (M2) Execution Report — Dataset Extraction & Pre-Filter QC

*Generated on: 2026-07-08 19:40:00 EDT*

## 1. Execution Summary
- **Authorized Milestone**: Milestone 2 (M2) — Dataset Extraction and Pre-Filter QC
- **Input Object Audited & Extracted**: `/local/projects-t3/lilab/vmenon/Pilot_tumor/processed_mpnst.rds`
- **Output Objects Created**: 4 dataset-specific, immutable working Seurat objects
- **HPC Execution Details**:
  - **SLURM Job ID**: 19161116
  - **Node**: `ihc-grid-1-1-1`
  - **State**: COMPLETED (ExitCode 0:0)
  - **Elapsed Time**: 6 minutes 35 seconds
  - **Requested Resources**: 8 CPUs, 64 GB RAM, 12 hours
  - **Peak Memory (MaxRSS)**: 66,935,856K (~63.83 GB) — Representing ~99.7% resource efficiency.
- **Git Commit Hash**: `472429e70d564c230ce08a2fb9d36e52ccadf0f0`

---

## 2. Extracted Constituents Inventory

The following table summarizes the four constituent datasets extracted from the original Seurat object by splitting on `sample_id`:

| Dataset ID | Cells | Genes | Default Assay | Size in Memory | MD5 Checksum | Output Path |
| --- | --- | --- | --- | --- | --- | --- |
| **MPNST_1** | 8,338 | 29,708 | SCT | 3.83 GB (4,108,548,872 B) | `f523e61322c585b16e632e800dc12bde` | `results/datasets/MPNST_1/MPNST_1_raw.rds` |
| **MPNST_2** | 2,830 | 29,708 | SCT | 1.19 GB (1,275,558,992 B) | `cd15a88b7f70a6492a1f05c94725cd8c` | `results/datasets/MPNST_2/MPNST_2_raw.rds` |
| **MPNST_3** | 3,682 | 29,708 | SCT | 1.42 GB (1,526,715,264 B) | `ecaca4b46e0a23f7bb8798cdb6adf280` | `results/datasets/MPNST_3/MPNST_3_raw.rds` |
| **MPNST_4** | 7,811 | 29,708 | SCT | 3.15 GB (3,380,746,280 B) | `f32811cf69595dcab9c5f5067efccbaa` | `results/datasets/MPNST_4/MPNST_4_raw.rds` |

### Assay and Layer Isolation
- All extracted objects preserve the assay structure (`RNA`, `SCT`, and `mnn.reconstructed`).
- For the `RNA` assay (class `Assay5`), the specific raw counts layers (`counts.MPNST_1` through `counts.MPNST_4`) and normalized layers (`data.MPNST_1` through `data.MPNST_4`) have been correctly retained.
- All dimensional reductions and legacy graphs/neighborhood tables have been dropped from the working objects to prevent downstream leakage or contamination.
- The `percent.ribo` column has been calculated for all objects using the standard ribosomal gene prefix pattern `^RP[SL]`.

---

## 3. Scientific Observations and QC recommendation Summary

For every constituent dataset, a detailed descriptive summary was generated in its respective `QC_RECOMMENDATION.md` report. Below is a comparison table of observed distributions and expected filtering impacts:

### Observed Median Values
| Dataset ID | Median nCount_RNA | Median nFeature_RNA | Median percent.mt | Median percent.ribo |
| --- | --- | --- | --- | --- |
| **MPNST_1** | 8,145 | 3,337 | 0.00% | 4.22% |
| **MPNST_2** | 7,657 | 2,419 | 7.10% | 15.67% |
| **MPNST_3** | 4,137 | 1,558 | 2.89% | 16.77% |
| **MPNST_4** | 5,244 | 2,179 | 8.09% | 15.52% |

### Expected Excluded Cells (Default Configured Thresholds)
*Configured thresholds: min_features=200, max_features=6000, min_counts=500, max_counts=50000, max_percent_mt=15.0%, max_percent_ribo=20.0%*

| Dataset ID | Cells Before | Expected Excluded (All Filters) | Expected Retained | Expected Removed % | Key Removal Driver |
| --- | --- | --- | --- | --- | --- |
| **MPNST_1** | 8,338 | 1,352 | 6,986 | 16.21% | `max_features` (> 6,000 features, 15.9%) |
| **MPNST_2** | 2,830 | 227 | 2,603 | 8.02% | `max_features` (> 6,000, 4.4%) & `percent.mt` (> 15%, 2.7%) |
| **MPNST_3** | 3,682 | 200 | 3,482 | 5.43% | `max_features` (> 6,000, 3.7%) & `percent.mt` (> 15%, 1.7%) |
| **MPNST_4** | 7,811 | 318 | 7,493 | 4.07% | `percent.mt` (> 15%, 2.2%) & `max_features` (> 6,000, 1.8%) |

### Key Biological Findings
1. **Mitochondrial Expression Anomalies**: In `MPNST_1`, the mitochondrial fraction (`percent.mt`) is exactly 0.00% across all cells. This indicates that either the library was sequenced/mapped without mitochondrial genes or they were excluded prior to object creation. In contrast, `MPNST_2`, `MPNST_3`, and `MPNST_4` have expected mitochondrial fractions (medians ~7.10%, 2.89%, and 8.09% respectively).
2. **High Complexity Cells in MPNST_1**: Nearly 16% of the cells in `MPNST_1` exceed 6,000 genes. This represents either a highly complex tumor cell population or potential doublets. The researcher must confirm whether a uniform 6,000 gene limit is appropriate or if it should be customized for `MPNST_1`.

---

## 4. Generated Artifacts List
The following files were successfully generated or updated during this milestone:

### Reports & Documentation
- [M2_REPORT.md](file://reports/milestones/M2_REPORT.md) (This file)
- [FIGURE_INDEX.tsv](file://reports/FIGURE_INDEX.tsv) (Combined figure index for publication plots)
- [CHANGELOG.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/CHANGELOG.md) (Updated with M2 details)
- [PROGRESS.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/PROGRESS.md) (Updated to complete M2)

### Dataset-Specific Subdirectories
For each dataset, the following files exist in `results/datasets/<id>` and `reports/datasets/<id>`:
- **`MPNST_1`**:
  - `results/datasets/MPNST_1/MPNST_1_raw.rds`
  - `results/datasets/MPNST_1/extraction_provenance.json`
  - `reports/datasets/MPNST_1/manifest.json`
  - `reports/datasets/MPNST_1/DATASET_VALIDATION.md`
  - `reports/datasets/MPNST_1/QC_RECOMMENDATION.md`
  - `reports/datasets/MPNST_1/qc_violins.pdf` / `qc_violins.png`
  - `reports/datasets/MPNST_1/qc_scatter.pdf` / `qc_scatter.png`
  - `reports/datasets/MPNST_1/qc_histograms.pdf` / `qc_histograms.png`
  - `reports/datasets/MPNST_1/qc_metrics.pdf` / `qc_metrics.png`
  - `reports/datasets/MPNST_1/qc_metrics_summary.tsv`
  - `reports/datasets/MPNST_1/qc_plots_provenance.json`
  - `reports/datasets/MPNST_1/qc_rec_provenance.json`
- **`MPNST_2`**, **`MPNST_3`**, and **`MPNST_4`**:
  - Equivalent validation reports, QC recommendations, plots, metrics, manifest, and provenance logs are in place.

---

## 5. Verification Status
- **Snakemake validation**: Successful
- **Unit Tests**:
  - `tests/unit/test_smoke.R`: Passed
  - `tests/unit/test_audit.R`: Passed
  - `tests/unit/test_extraction.R`: Passed
- **Validation check**: Verified that output RDS files do not contain dimensional reductions, preventing technical data leakage.

---

## 6. Researcher Decisions Required (Gate to Milestone 3)
Before proceeding to QC filtering in Milestone 3 (M3), the researcher must approve:
1. **The 6,000 feature limit on `MPNST_1`**: Since this removes 15.9% of all cells in `MPNST_1`, we should confirm whether these are true tumor cells with high transcriptional activity or doublets.
2. **The 15% mitochondrial limit on `MPNST_2`, `MPNST_3`, and `MPNST_4`**.
3. **The 20% ribosomal limit across all datasets**.

**DO NOT PROCEED TO M3 OR FILTER ANY CELLS UNTIL APPROVAL IS GRANTED.**
