[M0]

Infrastructure

* Environment creation
* Snakemake installation
* Synthetic workflow
* SLURM validation
* Logging
* Provenance

[M1]

Real Data Audit

* Object audit
* Metadata audit
* Object inventory
* Raw count verification
* Dataset identifier validation
* Legacy analysis inventory

[M2]

Dataset Extraction and Pre-Filter QC

* Deterministic dataset extraction
* Immutable dataset-specific Seurat objects
* Dataset validation reports (DATASET_VALIDATION.md)
* Pre-filter QC figures and tables (nFeature_RNA, nCount_RNA, percent.mt, percent.ribo, etc.)
* QC recommendation reports (QC_RECOMMENDATION.md)
* Updated FIGURE_INDEX.tsv
* Updated Snakemake workflow and unit/integration tests

[M3]

QC Filtering and Doublet Assessment

* Dynamic command-line configuration for QC filters
* Reproducible doublet detection with scDblFinder and full fallbacks
* QC filtering applying researcher-approved thresholds (200 features, 500 counts, 10% MT, 20% ribosomal)
* Retention of high-feature singlets to preserve transcriptional complexity
* Post-filter QC figures (violins, scatter, density, histograms, doublet summaries, filtering summaries)
* Dataset-specific filter reports (FILTER_REPORT.md)
* Cross-dataset comparison analysis and M3 milestone report
* Updated FIGURE_INDEX.tsv
* Updated Snakemake workflow and validation unit tests
* Dataset-specific QC threshold framework (Strategy B) in addition to global Strategy A
* Dynamic R script support for custom file suffixes to allow parallel runs
* Threshold sensitivity simulations and grid intersection overlap plots (qc_optimization_review)
* Quantitative comparative analysis of strategies (compare_qc_strategies)
* Strategy comparison report (reports/QC_COMPARISON_REPORT.md)
* Comprehensive validation tests for dataset-specific filtering (test_filtering_specific)
* Fully consolidated global figure index (FIGURE_INDEX.tsv) indexing all 32 diagnostic plots

[M4]

Normalization and Variable Features

* **Added**:
  * Configurable normalization and variable feature selection framework supporting `SCTransform` v2 and `LogNormalize`.
  * Highly variable feature (HVF) selection sorting and metrics reporting (`reports/datasets/MPNST_*/variable_features.tsv`).
  * Variable feature diagnostic plots (scatter, metric distribution, and top 6 expression violins in PDF/PNG) under `reports/datasets/MPNST_*/`.
* **Changed**:
  * Upgraded Snakemake workflow to include normalization and variable feature discovery rules.
* **Scientific decisions**:
  * Selected `SCTransform` v2 as the primary normalization strategy to handle sequencing depth variations and stabilize variance in MPNST datasets.
  * Dynamically handle libraries lacking mitochondrial variance (e.g., `MPNST_1`) by omitting mitochondrial regression to prevent numerical singularity errors.
* **Workflow changes**:
  * Added portable Snakemake rules: `normalize_and_find_features`, `test_normalization`, and `generate_m4_report`.
* **Validation**:
  * Created automated validation unit test suite `tests/unit/test_normalization.R` to check normalised object structure and diagnostic reports.
  * Added Snakemake verification rule `test_normalization` running on synthetic data.
* **Outputs**:
  * Seurat objects: `results/datasets/MPNST_*/MPNST_*_normalized.rds`.
  * Reports: `reports/datasets/MPNST_*/NORM_REPORT.md` and consolidated `reports/milestones/M4_REPORT.md`.
  * Diagnostic figures: `var_features_scatter.png`, `var_features_distribution.png`, and `top_features_violins.png`.

[M5]

Principal Component Analysis and Evaluation

* **Added**:
  * Configurable PCA and PC evaluation rule supporting SCTransform-derived scale.data.
  * Non-arbitrary, quantitative PC range selection using a geometric elbow detector (recommended range), marginal variance drop < 1.0% (conservative range), and variance noise floor < 0.3% (maximum range).
  * Correlation testing of leading PC scores against technical covariates (`nCount_RNA`, `nFeature_RNA`, `percent.mt`, `percent.ribo`) reporting Pearson/Spearman coefficients and statistical significance.
  * Premium quality publication figure generation (Elbow plots, cumulative variance plots, PC loading barplots, PC DimHeatmaps, and technical correlation heatmaps) in PDF and PNG.
  * Machine-readable recommendation files (`reports/PCA_RECOMMENDATIONS.tsv`).
* **Changed**:
  * Integrated Snakemake rules: `run_pca_and_evaluation`, `test_pca`, and `generate_m5_report`.
  * Updated global figure index (`reports/FIGURE_INDEX.tsv`) capturing all 40 newly generated PCA figures.
* **Scientific decisions**:
  * Performed PCA independently per dataset to avoid batch mixing before clustering evaluation.
  * Omitted JackStraw permutation checks due to statistical incompatibility with regularized negative binomial SCT z-scored residuals and high computational resource footprint.
* **PCA methodology**:
  * Computed PCA on the top 3,000 variable features using Seurat's RunPCA on the SCT assay scaled residuals.
* **PC-selection methodology**:
  * Used geometric elbow/knee-point detection to define the Recommended PC limits: `PC8` (MPNST_1), `PC6` (MPNST_2), `PC9` (MPNST_3), and `PC5` (MPNST_4).
* **Validation**:
  * Automated unit test suite `tests/unit/test_pca.R` validating PC coordinates, loadings, reports, and provenance.
* **Outputs**:
  * Seurat objects: `results/datasets/MPNST_*/MPNST_*_pca.rds`.
  * Reports: `reports/datasets/MPNST_*/PCA_REPORT.md` and consolidated `reports/milestones/M5_REPORT.md`.
  * Tables: `reports/datasets/MPNST_*/top_loading_genes.tsv`, `reports/datasets/MPNST_*/pc_technical_correlations.tsv`, and `reports/datasets/MPNST_*/pca_variance_explained.tsv`.
  * Diagnostic figures: `pca_elbow.png`, `pca_cumulative_variance.png`, `pca_loadings.png`, `pca_heatmaps.png`, and `pca_correlations.png`.
* **HPC execution**:
  * SLURM production job ID `19176730` completed with ExitCode 0:0, elapsed time 00:03:51, and MaxRSS ~10.88 GB.
