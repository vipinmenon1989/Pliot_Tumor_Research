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

* Configurable normalization and variable feature selection framework supporting SCTransform v2 and LogNormalize
* Dynamic regression handling for libraries lacking mitochondrial variance (e.g. MPNST_1) to prevent numerical failures
* Highly variable feature (HVF) selection sorting and metrics reporting (variable_features.tsv)
* Variable feature diagnostic plots (scatter, metric distribution, and top 6 expression violins in PDF/PNG)
* Portable Snakemake rules (`normalize_and_find_features`, `test_normalization`, `generate_m4_report`)
* Consolidated Milestone 4 execution report (M4_REPORT.md)
* Automated validation unit test suite (test_normalization.R)

[M5]

Principal Component Analysis and Evaluation

* Configurable PCA and PC evaluation rule supporting SCTransform-derived scale.data.
* Non-arbitrary, quantitative PC range selection using a geometric elbow detector (recommended range), marginal variance drop < 1.0% (conservative range), and variance noise floor < 0.3% (maximum range).
* Correlation testing of leading PC scores against technical covariates (nCount_RNA, nFeature_RNA, percent.mt, percent.ribo) reporting Pearson/Spearman coefficients and statistical significance.
* Omission of JackStraw permutation checks with documented R-based and single-cell regression justifications.
* Premium quality publication figure generation (Elbow plots, cumulative variance plots, PC loading barplots, PC DimHeatmaps, and technical correlation heatmaps) in PDF and PNG.
* Machine-readable recommendation files (PCA_RECOMMENDATIONS.tsv) and dataset-specific PCA reports (PCA_REPORT.md).
* Integrated Snakemake rules (`run_pca_and_evaluation`, `test_pca`, `generate_m5_report`) and automated unit test suite (test_pca.R) passing active dataset parameters.
* Updated global figure index (FIGURE_INDEX.tsv) capturing all 40 newly generated figures.

