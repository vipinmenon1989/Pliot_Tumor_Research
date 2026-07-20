# Milestone 0 to Milestone 8 Reconciliation Audit Report
*Generated on: 2026-07-20 12:00:00*

---

## 1. Audit Summary
This report performs a comprehensive, three-layer verification pass across all previously completed milestones (Milestones M0 through M8). Each milestone was evaluated for:
1. **Execution Evidence**: Verifying code, Snakemake rules, tests, logs, SLURM job execution records, output artifacts, and benchmark TSVs.
2. **Scientific/Technical Reporting**: Verifying consolidated milestone reports, dataset-specific markdown files, summary TSVs/JSONs, execution provenance records, and `reports/FIGURE_INDEX.tsv` registrations.
3. **Public Project Documentation**: Verifying consistency across `README.md`, `PROGRESS.md`, `CHANGELOG.md`, and `PROJECT.md`.

---

## 2. Reconciliation Matrix

The table below summarizes the results of the three-layer verification audit. (The machine-readable version is available at [reports/audits/M0_M8_RECONCILIATION.tsv](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/audits/M0_M8_RECONCILIATION.tsv)):

| Milestone | Execution Verified | Required Outputs Verified | Report Exists | Report Complete (PROJECT.md) | PROGRESS Consistent | CHANGELOG Consistent | README Consistent | FIGURE_INDEX Consistent | Provenance Consistent | Unresolved Issues | Severity | Resolution |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- | :---: | :--- |
| **M0** | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | None | None | None |
| **M1** | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | None | None | None |
| **M2** | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | None | None | None |
| **M3** | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | None | None | None |
| **M4** | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | None | None | None |
| **M5** | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | None | None | None |
| **M6** | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | None | None | None |
| **M7** | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | None | None | None |
| **M8** | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | None | None | None |

---

## 3. Layer-by-Layer Verification Details

### M0 — Infrastructure, Environment, and SLURM Safety
- **Execution Evidence**: Checked Snakemake execution engine, Conda portable specification `workflow/envs/R_env_portable.yaml`, smoke tests, and the initial SLURM test log.
- **Reporting**: Verified [M0_REPORT.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/milestones/M0_REPORT.md) exists and is complete.
- **Documentation**: Verified M0 is marked as completed in all docs.

### M1 — Real-Data Audit
- **Execution Evidence**: Verified read-only audit script `scripts/R/audit_object.R` and unit tests in `tests/unit/test_audit.R`.
- **Reporting**: Verified [M1_REPORT.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/milestones/M1_REPORT.md), [DATA_AUDIT.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/DATA_AUDIT.md), and [object_inventory.tsv](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/object_inventory.tsv).
- **Documentation**: Checked consistency across all documentation.

### M2 — Dataset Extraction and Pre-Filter QC
- **Execution Evidence**: Checked extraction script `scripts/R/extract_dataset.R` and generated raw RDS files `results/datasets/MPNST_*/MPNST_*_raw.rds`.
- **Reporting**: Checked dataset-specific validation reports [DATASET_VALIDATION.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/datasets/MPNST_1/DATASET_VALIDATION.md) and QC recommendations.
- **Documentation**: Verified entries in `FIGURE_INDEX.tsv`.

### M3 — QC Filtering and Doublet Assessment
- **Execution Evidence**: Verified scDblFinder doublet removal execution and dataset-specific threshold filtering in `scripts/R/filter_and_detect_doublets.R`.
- **Reporting**: Checked [M3_REPORT.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/milestones/M3_REPORT.md), doublet reports, and the strategy comparison report [QC_COMPARISON_REPORT.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/QC_COMPARISON_REPORT.md).
- **Documentation**: Reconciled all QC figures in the global figure index.

### M4 — Normalization and Variable Features
- **Execution Evidence**: Verified SCTransform v2 normalization and 3,000 highly variable features selection.
- **Reporting**: Checked consolidated report [M4_REPORT.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/milestones/M4_REPORT.md) and dataset-specific reports.
- **Documentation**: Checked project logs.

### M5 — PCA and PC Evaluation
- **Execution Evidence**: Verified independent PCA run on variable features and elbow determination.
- **Reporting**: Checked PCA recommendations table [PCA_RECOMMENDATIONS.tsv](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/PCA_RECOMMENDATIONS.tsv) and consolidated [M5_REPORT.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/milestones/M5_REPORT.md).
- **Documentation**: Verified PC numbers matched recommendations.

### M6 — Neighbors, UMAP, and Clustering Resolution Sweep
- **Execution Evidence**: Checked Louvain clustering sweeps from 0.1 to 1.0 (in increments of 0.1) on dataset-specific PCA spaces.
- **Reporting**: Checked consolidated report [M6_REPORT.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/milestones/M6_REPORT.md) and recommendation summaries.
- **Documentation**: Reconciled UMAP plots in the global index.

### M7 — Marker Discovery and Dataset Recommendations
- **Execution Evidence**: Verified Wilcoxon rank-sum marker discovery across all 40 combinations and downsampled heatmap/dotplot/FeaturePlot visualizations.
- **Reporting**: Checked consolidated [M7_REPORT.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/milestones/M7_REPORT.md) and recommendations [ANALYSIS_RECOMMENDATION.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/datasets/MPNST_1/ANALYSIS_RECOMMENDATION.md).
- **Documentation**: Reconciled the 24 recommended-resolution visualizations and their figure index entries.

### M8 — Combined Pre-Integration Baseline
- **Execution Evidence**: Checked combined pre-integration Seurat object containing 19,716 cells, shared PCA (PCs 1-30), shared UMAP (`umap_preintegration`), and block-based neighbor mixing diagnostics.
- **Reporting**: Checked expanded [M8_REPORT.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/milestones/M8_REPORT.md), assessment [PRE_INTEGRATION_ASSESSMENT.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/PRE_INTEGRATION_ASSESSMENT.md), and design [INTEGRATION_PREPARATION.md](file:///local/projects-t3/lilab/vmenon/Pilot_tumor/reports/INTEGRATION_PREPARATION.md).
- **Documentation**: Verified that all stale references, unsupported scientific patient-specific claims, and trailing whitespaces were corrected.

---

## 4. Verification Conclusion
- **Critical or High Issues**: **None**.
- **Conclusion**: All previous milestones (M0–M8) are fully completed, technically validated, scientifically reportable, and consistently documented. The project is ready to proceed to the final freeze.
