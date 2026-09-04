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

[M6]

Clustering Resolution Sweep and Dataset-Specific Cluster Selection

* **Added**:
  * Configurable clustering resolution sweep rule running Louvain clustering (Algorithm 1) from resolution 0.1 to 1.0 in steps of 0.1 on independent SNN graphs.
  * Graph construction using recommended PC dimensions inherited from Milestone 5: `PC1:8` (MPNST_1), `PC1:6` (MPNST_2), `PC1:9` (MPNST_3), and `PC1:5` (MPNST_4).
  * Subsampling bootstrap stability analysis (5 rounds of 80% cells) using corrected cell-order alignment to calculate true Adjusted Rand Index (ARI) metrics.
  * Robust technical covariate correlation check via linear model $R^2$, including safety checks to prevent `NaN` crashes on constant (zero-variance) covariates (e.g. `percent.mt` in `MPNST_1`).
  * Trade-off recommendation heuristic targeting standard single-cell resolutions (0.3 to 0.8), filtering for min cluster size $\ge 5$ and stability ARI $\ge 0.60$, and optimizing for high stability and low technical covariate correlation.
  * Multi-resolution UMAP grids, recommended resolution UMAP projection, sweep metric plots, stability/covariate correlation plots, and ggplot2 transition trees (dendrogram-like cluster splits) in PDF and PNG.
  * Machine-readable recommendation files (`reports/CLUSTERING_RECOMMENDATIONS.tsv` and `reports/CLUSTERING_SWEEP_SUMMARY.tsv`).
* **Changed**:
  * Integrated Snakemake rules: `run_clustering_sweep`, `test_clustering`, and `generate_m6_report`.
  * Updated global figure index (`reports/FIGURE_INDEX.tsv`) capturing all 20 newly generated clustering sweep figures.
* **Scientific decisions**:
  * Performed graph construction and clustering independently on each dataset to prevent technical integration artifacts.
  * Adopted the Louvain algorithm as a robust modularity-maximization method due to the lack of Python's `leidenalg` in the R HPC environment.
  * Recommended biologically-defensible, high-stability resolutions for downstream marker discovery: `0.6` (MPNST_1, 18 clusters, stability ARI = 0.919), `0.3` (MPNST_2, 9 clusters, stability ARI = 0.933), `0.6` (MPNST_3, 13 clusters, stability ARI = 0.910), and `0.7` (MPNST_4, 14 clusters, stability ARI = 0.740).
  * Flagged a technical MT percentage covariate correlation bias ($R^2 = 0.47$) in `MPNST_4` at the recommended resolution of 0.7 for monitoring in downstream analyses.
* **Validation**:
  * Automated unit test suite `tests/unit/test_clustering.R` validating metadata columns, cell counts, coordinates, graphs, and reports.
  * Added Snakemake verification rule `test_clustering` running on synthetic data.
* **Outputs**:
  * Seurat objects: `results/datasets/MPNST_*/MPNST_*_clustered.rds`.
  * Reports: `reports/datasets/MPNST_*/CLUSTERING_REPORT.md` and consolidated `reports/milestones/M6_REPORT.md`.
  * Tables: `reports/datasets/MPNST_*/clustering_sweep_stats.tsv`, `reports/datasets/MPNST_*/clustering_stability_metrics.tsv`, and `reports/datasets/MPNST_*/clustering_recommendation_summary.tsv`.
  * Diagnostic figures: `pca_umap_grid.png`, `umap_recommended.png`, `clustering_metrics.png`, `clustering_stability.png`, and `clustering_tree.png`.
* **HPC execution**:
  * SLURM production job ID `19399140` completed with ExitCode 0:0, elapsed time 00:07:58, and MaxRSS ~28.22 GB.

[M7]

Marker Discovery and Dataset-Specific Recommendations

* **Added**:
  * Wilcoxon rank-sum marker discovery scripts running independently across all 4 datasets and 10 resolutions (40 combinations).
  * FeaturePlot generators and dot plot generators for recommended resolutions, downsampling cell count to a maximum of 100 cells per cluster to construct clean, publication-ready heatmaps.
  * Publication-quality UMAP FeaturePlot grid panel (2x3 format) for 5-6 top representative markers per dataset at its recommended resolution.
  * Standardized local figure index files (`figure_index_m7.tsv`) for recommended resolution visualizations.
  * Dynamically queries active assay (`SCT` vs standard `RNA`) and executes `PrepSCTFindMarkers` if SCT is active.
  * Validation unit tests (`tests/unit/test_markers.R`) verifying marker presence, specificity metrics, column ranges, and checking all 24 required visual files (heatmap PDF/PNG, dotplot PDF/PNG, FeaturePlot PDF/PNG across 4 datasets) and their global figure index registration.
  * Global consolidated report (`reports/milestones/M7_REPORT.md`) and dataset-specific recommendation reports (`reports/datasets/{ds}/ANALYSIS_RECOMMENDATION.md`).
* **Changed**:
  * Integrated rules `discover_markers`, `visualize_markers`, `generate_m7_report`, and `test_markers` into `workflow/Snakefile`.
  * Expanded Snakemake `visualize_markers` rule to explicitly list all 6 figure paths (dotplot, heatmap, representative FeaturePlots in PDF/PNG) plus the local index as target outputs to guarantee they are tracked correctly.
  * Added `reports/FIGURE_INDEX.tsv` as an explicit input to the Snakemake `test_markers` rule to enforce a strict DAG dependency and eliminate parallel race conditions.
  * Consolidated global `reports/FIGURE_INDEX.tsv` to register all 12 recommended-resolution marker figures with metadata and git commits.
* **Scientific decisions**:
  * Assessed 40 dataset-resolution combinations to evaluate marker quality (median markers per cluster) and check for weak cluster support (fewer than 5 distinct markers) and small clusters (fewer than 10 cells).
  * Confirmed that recommended resolutions provide robust, biologically-relevant marker support with zero weak or small clusters.
  * Addressed mitochondrial bias in MPNST_4 (correlation with `percent.mt` $R^2 = 0.47$ at recommended resolution 0.7) and proposed resolution 0.5 (reducing correlation to $R^2 = 0.24$) as the primary alternative baseline.
* **HPC execution**:
  * SLURM production job ID `19399784` completed with Snakemake execution of clustering sweep, marker sweep, visualization, reporting, and validation in 11m 8s, utilizing ~28.43 GB of memory.
  * SLURM production job ID `19399799` executed the standardized visualization, local index generation, and global index merge in 2m 21s, using MaxRSS 5.22 GB.

[M8]

Combined Pre-Integration Baseline and Integration-Readiness Assessment

* **Added**:
  * Unified combination script `scripts/R/combine_pre_integration.R` that merges the four constituent datasets, renames cells, namespaces resolution columns (`preint_MPNST_{ds}_res_{resolution}`), checks cell and feature counts, and outputs metadata inventories and dictionaries.
  * Unified pre-integration analysis script `scripts/R/analyze_pre_integration.R` that executes global SCTransform normalization, PCA, neighbor graph construction, UMAP calculation (`umap_preintegration`), neighborhood dataset-mixing diagnostics, composition analysis, and generates reports and visualizations.
  * Custom block-based matrix distance neighbor mixing diagnostic algorithm in pure R, calculating same-dataset neighbor fraction and Shannon entropy for each cell using a memory-efficient matrix approach without external package dependencies.
  * Validation unit test suite `tests/unit/test_preintegration.R` asserting cell counts per dataset, uniqueness of cell names, namespaced resolution columns with proper `NA` mappings, dimensional reductions, absence of integrated embeddings, reports existence, and global figure index registration.
  * SLURM execution script `scripts/shell/run_m8_workflow.sh` configured for partition `ihc` node `ihc-grid-1-1-1`.
  * Consolidated global reports: `reports/PRE_INTEGRATION_ASSESSMENT.md`, `reports/INTEGRATION_PREPARATION.md`, and consolidated milestone report `reports/milestones/M8_REPORT.md`.
  * Local figure index `reports/combined/pre_integration/figure_index_m8.tsv` and 10 pre-integration visualizations (PDF and PNG formats, total 20 visual files) including elbow, dataset PCA, dataset/sample/QC UMAPs, faceted views, recommended clusters UMAP, and neighbor mixing boxplots.
* **Changed**:
  * Integrated rules `combine_pre_integration`, `run_pre_integration_analysis`, and `test_preintegration` into `workflow/Snakefile`.
  * Updated global `reports/FIGURE_INDEX.tsv` to register all 10 Milestone 8 pre-integration figures.
* **Scientific decisions**:
  * Established a mathematically valid shared non-integrated baseline space for Phase 1 datasets by running a single unified SCTransform on the merged raw counts.
  * Prohibited dataset integration (Harmony, CCA, RPCA, fastMNN, scVI, etc.) or dimensional reduction concatenation to preserve a true baseline.
  * Audited clinical covariates, documenting that patient/sample identity is 100% confounded with dataset identity due to missing patient demographics.
  * Quantified dataset segregation via neighborhood diagnostics, showing an extremely high mean same-dataset neighbor fraction (>98%) across all datasets.
  * Provided a scientific justification and proposed design for Phase 2 integration (Harmony, CCA, and RPCA benchmarking).
* **HPC execution**:
  * SLURM production job ID `19403199` completed with Snakemake execution of combination, normalization, embedding, mixing diagnostics, visualizations, reports, and validation tests in 17m 17s, utilizing ~32.00 GB of memory.

[M9]

Workflow Hardening, CI/CD, Provenance Finalization, Phase 1 Handoff, and Phase 1 Freeze

* **Added**:
  - Top-level canonical Snakemake target `phase1_complete` which pulls in the Phase 1 manifest, handoff report, milestone report, audits, and consistency logs.
  - Portable preflight environment validator script `scripts/python/preflight_checker.py` checking packages, configuration, and permissions.
  - Python project-state consistency checker `validate_project_state.py` validating PROGRESS, CHANGELOG, README, and milestone report file synchronizations.
  - Programmatic Phase 1 manifest generator `scripts/python/generate_phase1_manifest.py` generating SHA256 checksums and parameters for all RDS objects.
  - GitHub Actions CI workflow configuration `.github/workflows/ci.yml` performing static, config, linting, and clean-room synthetic workflow testing.
  - Human-readable final handoff documentation `reports/PHASE1_HANDOFF.md`.
  - Consolidated Milestone 9 execution report `reports/milestones/M9_REPORT.md`.
  - Three-layer verification audit report `reports/audits/M0_M8_RECONCILIATION.md` and machine-readable `M0_M8_RECONCILIATION.tsv`.
* **Changed**:
  - Updated `workflow/Snakefile` adding rules `generate_phase1_manifest`, `validate_project_state`, and `phase1_complete`, and linking them to rule `all`.
  - Fully synchronized `README.md` and `PROGRESS.md` to reflect M0-M9 complete, Phase 1 frozen, and Phase 2 not started status.
* **Scientific decisions**:
  - Declared Phase 1 frozen and established the entry contract for Phase 2, defining the immutable `preint_*` namespace and requiring Phase 2 coordinates to use `postint_*`.
* **HPC execution**:
  - SLURM validation job ID `19403619` executed the Snakemake verification tests on the compute node in 1m 15s, using MaxRSS 16.00 GB.



[M10]

Phase 2 Reconstruction and Integration Planning

*2026-09-02 — Phase 2 begins. No Harmony executed.*

* **Added**:
  - Phase 2 directory skeleton: `results/phase2/{harmony,clustering,markers,annotation,composition,handoff}`, `reports/phase2/{milestones,figures}`, `scripts/R/phase2/`, `scripts/shell/phase2/`, `logs/phase2/{slurm,milestones}`, `benchmarks/phase2/`.
  - Read-only Phase 1 handoff validation script `scripts/R/phase2/inspect_phase1_handoff.R`, which computes md5 and sha256 of the frozen object, enumerates assays/layers/SCT models/reductions/graphs/metadata, cross-tabulates candidate batch variables, and summarises per-sample QC — without modifying any Phase 1 artefact.
  - SLURM launcher `scripts/shell/phase2/run_m10_inspect.sh` (ihc / ihc-grid-1-1-1, 4 CPUs, 96G, 02:00:00).
  - Modular Phase 2 plan `reports/phase2/PHASE2_PLAN.md`: input contract, `postint_*` namespace contract, per-milestone script/SLURM module map, native mixing-metric design, carried-forward prohibitions.
  - Harmony grouping-variable decision document `reports/phase2/HARMONY_VARIABLE_DECISION.md`, including candidate screening of all 93 metadata columns, cross-tabulations, over/under-correction risk analysis, the exact proposed `RunHarmony()` call, an explicit-vs-default parameter table, and a researcher sign-off block.
  - Milestone report `reports/phase2/milestones/M10_REPORT.md`.
  - Environment records `reports/phase2/PHASE2_ENVIRONMENT.tsv` and `reports/phase2/PHASE2_SESSIONINFO.txt`.
  - Eight machine-readable inspection outputs under `results/phase2/handoff/`.
* **Changed**:
  - `PROGRESS.md` restructured: Phase 2 section added ahead of the frozen Phase 1 historical record; status line changed from `PHASE 2 = NOT STARTED` to `PHASE 2 = IN PROGRESS`.
  - `README.md` Phase 2 status synchronised.
  - `scripts/python/validate_project_state.py`: the Phase 2 status check now requires an *explicit* `PHASE 2 = ...` declaration (`NOT STARTED` | `IN PROGRESS` | `FROZEN` | `COMPLETE`) in `PROGRESS.md` and `README.md`, instead of the hard-coded `NOT STARTED` that became false when Phase 2 began.
* **Workflow decisions**:
  - **Snakemake is not used in Phase 2.** Each milestone is a standalone argument-driven R script under `scripts/R/phase2/` invoked by an explicit `sbatch` script under `scripts/shell/phase2/`.
  - Phase 1's `postint_*` entry contract is reconciled with the Phase 2 specification's reduction/metadata names by prefixing: `postint_harmony`, `postint_umap_harmony`, `postint_harmony_clusters_res_0.1`…`_res_1.0`, `postint_celltype_level1/2/3`, `postint_annotation_confidence`.
  - `lisi`, `kBET` and `clustree` are **not installed**; the M12 mixing metrics and the M13 transition diagram are implemented natively, extending the block-wise kNN routine validated in Phase 1 `analyze_pre_integration.R`. No package was installed, upgraded or removed.
* **Scientific decisions**:
  - Phase 2 input confirmed as `results/combined/pre_integration/combined_preintegration.rds`; both Phase 1 checksums (md5 `88a442688f912d882f6c6da01820e329`, sha256 `c3fdce8b…dc66`) reproduce exactly.
  - Proposed Harmony grouping variable `sample_id`; proposed dimensions `1:30` of `pca`; all Harmony scientific parameters left at `harmony 1.2.4` defaults. **Proposed, not executed** — Phase 2 §8 requires researcher approval because dataset = sample = patient = presumed technical batch is a single inseparable 4-level variable.
  - The frozen Phase 1 `pca` and `umap_preintegration` reductions are adopted as the pre-Harmony baseline and preserved unmodified; no duplicate copy is created.
  - Legacy `orig.anno`, legacy clusters and legacy Harmony/CCA/RPCA/MNN embeddings inventoried for provenance only and excluded from all parameter selection.
* **Findings reported, not repaired**:
  - The Phase 1 `SCT` assay contains **four** SCTransform models, not the single global model described in `PRE_INTEGRATION_ASSESSMENT.md` §3.2 and `PHASE1_HANDOFF.md` §3.1 — Seurat v5 ran SCTransform once per merged `RNA` layer. Integration remains valid; `PrepSCTFindMarkers()` becomes mandatory in M14.
  - `percent.mt` is identically zero across all 7,615 MPNST_1 cells while Phase 1 M8 passed `vars.to.regress = "percent.mt"` unconditionally. A non-finite-value integrity check is added to the front of M11.
  - Nine Phase 1 prose-vs-record documentation defects (D1–D9) catalogued in `M10_REPORT.md` §7, including `config/config.yaml: input_rds` pointing at the deleted `/local/projects-t3/lilab/vmenon/Pilot_tumor/` project root.
* **HPC execution**:
  - SLURM job `19886411` (`p2_m10_inspect`, ihc-grid-1-1-1, 4 CPUs, 96G, 02:00:00): **COMPLETED**, ExitCode 0:0, Elapsed 00:02:17, MaxRSS 4,959,508K (4.73 GiB), TotalCPU 00:01:30, empty stderr. Memory efficiency 4.9% — the M11 request (8 CPUs / 64G / 04:00:00) is sized from this observation and the Phase 1 M8 benchmark rather than from the resource envelope maximum.
* **STOP**: Phase 2 halted at the M10 gate pending researcher approval of the Harmony grouping variable, dimensions and naming contract.

[M11]

Default Harmony Integration

*2026-09-02 — Harmony executed on the approved configuration. Integration quality NOT yet accepted.*

* **Added**:
  - `scripts/R/phase2/run_harmony_integration.R` — modular, argument-driven M11 Harmony script. Enforces six pre-run guards (output-path safety against every protected Phase 1 location; md5 + sha256 verification of the input against the M10 record; cell count and duplicate-barcode check; grouping-variable level and per-level count check against M10; reduction/dimension availability and row-name ordering; numerical-integrity precheck for non-finite values). After Harmony it proves Phase 1 preservation by comparing pre/post digests of the `pca` and `umap_preintegration` embeddings, the metadata frame, column names and cell names, then re-reads the saved RDS from disk and re-validates 17 structural properties. Aborts with a non-zero exit status on any failure. Includes a `--validation-mode` flag that relaxes only the real-data-specific assertions, for synthetic smoke testing.
  - `scripts/shell/phase2/run_m11_harmony.sh` — SLURM launcher (account/partition `ihc`, node `ihc-grid-1-1-1`, 8 CPUs, 64G, 04:00:00), sized from M10 accounting rather than the resource envelope maximum.
  - `results/phase2/harmony/phase2_harmony_integrated.rds` — the Phase 2 object (new file; the Phase 1 object was never opened for writing).
  - `results/phase2/harmony/harmony_parameters.json` and `harmony_parameters.tsv` — full parameter record separating explicitly supplied arguments from package defaults.
  - `results/phase2/harmony/harmony_convergence.tsv` and `harmony_kmeans_objective.tsv` — Harmony objective per iteration and per clustering step.
  - `results/phase2/harmony/harmony_embedding_dimension_summary.tsv` and `harmony_grouping_composition.tsv` — machine-readable embedding and grouping summaries.
  - `results/phase2/harmony/prov_m11_harmony.json` — provenance written through the existing Phase 1 `scripts/R/provenance_utils.R` helpers rather than a parallel system.
  - `reports/phase2/milestones/M11_REPORT.md`.
* **Changed**:
  - `PROGRESS.md` — M11 section added; Phase 2 status advanced to "M11 complete, STOPPED awaiting authorization before M12".
  - `CHANGELOG.md` — this entry.
  - `reports/FIGURE_INDEX.tsv` — **deliberately unchanged**; M11 generated no figures.
* **Scientific decisions**:
  - Harmony executed exactly as approved at the M10 gate: `harmony::RunHarmony(object = obj, group.by.vars = "sample_id", reduction.use = "pca", dims.use = 1:30, reduction.save = "postint_harmony", verbose = TRUE)` with `set.seed(42)`.
  - **Every scientific parameter left at harmony 1.2.4 defaults** — `theta` 2, `sigma` 0.1, `lambda` 1, `nclust` 100 (the automatic `min(round(N/30), 100)` cap), `max_iter` 10, `early_stop` TRUE, `ncores` 1, `project.dim` TRUE, and `harmony_options()` untouched. Nothing was tuned, and no parameter was chosen on the basis of appearance.
  - Harmony was applied at the embedding level only. No RNA or SCT expression value was modified.
  - A UMAP on the Harmony reduction was **deliberately deferred to M12/M13**: it is not required to validate that Harmony executed successfully, and the milestone scope restricts M11 to diagnostics.
  - Convergence provenance: because `RunHarmony.Seurat` discards Harmony's internal state, a second matrix-level call with identical input, seed and defaults was used to recover the objective history, and its embedding was asserted **bit-identical** to the primary result (max abs difference 0.000e+00). The Harmony result is confirmed deterministic under seed 42.
* **Results**:
  - Harmony **converged after 9 of a maximum 10 iterations** (`early_stop` triggered); objective 753.089 → 353.656 across eight iterations, +0.35% at iteration 9.
  - Output object: 19,716 cells (in = out), reductions `pca` + `umap_preintegration` + `postint_harmony` (19,716 × 30, key `postintharmony_`, assay `SCT`), graphs `SCT_nn`/`SCT_snn`, 93 metadata columns, sha256 `6085976f788633b68c07c1ddff49e20e6f6f6727579cb7900854febd2ae6e69d`.
  - **Phase 1 state proven preserved byte-for-byte**, both in memory and after a disk round-trip.
  - **M10 Finding 2 CLOSED**: the zero-variance `percent.mt` regression for MPNST_1 produced 0 non-finite values in `pca` (all 50 dims), 0 in `SCT@scale.data` (5,192 × 19,716; none in any sample) and 0 in the Harmony embedding.
  - Incidental finding for M14: `SCT@scale.data` spans **5,192** genes rather than the 3,000 requested in Phase 1 M8, because the four SCT models contribute a union of variable features.
* **Scientific status**: Harmony execution has completed, but **integration quality has not yet been scientifically accepted**. A falling objective function is an optimiser meeting its own criterion, not evidence of correct integration. Over-correction of the patient-private malignant compartment and of the sample-restricted populations flagged at M10 remains untested. Pre/post-Harmony evaluation is reserved for M12.
* **HPC execution**:
  - SLURM job `19886486` (`p2_m11_harmony`, ihc-grid-1-1-1, 8 CPUs, 64G, 04:00:00): **COMPLETED**, ExitCode 0:0, Elapsed 00:08:24, MaxRSS 8,788,416K (8.38 GiB, 13.1% of request), TotalCPU 00:07:58. Zero R warnings (`options(warn = 1)`), zero R errors; stderr contains only Harmony's own progress messages. Harmony itself accounted for 10.1 s of the runtime; `saveRDS` (256.7 s) and checksumming (110 s) of the 5.6 GB object dominated.
  - Before submission the script was validated end-to-end on a synthetic structural mirror of the handoff object (400 cells, `pca` with 50 dims, `umap_preintegration`, `SCT_nn`/`SCT_snn`, 4-level `sample_id`), and the abort guards were confirmed to fire for a protected output path, an existing output file, and a wrong cell count. No real object was loaded on the login node.
* **STOP**: Phase 2 halted at the M11 gate pending researcher authorization for M12.

[M12]

Pre/Post Harmony Evaluation

*2026-09-02 — Default Harmony assessed. Recommendation: ACCEPT DEFAULT HARMONY WITH CAVEATS.*

* **Added**:
  - `scripts/R/phase2/evaluate_harmony.R` — M12 evaluation script. Verifies the M11 object's md5/sha256, builds a post-Harmony UMAP with parameters matched exactly to the frozen Phase 1 UMAP, re-checks that the Phase 1 PCA/UMAP and the M11 Harmony embedding are unchanged, then releases the Seurat object and computes all metrics on extracted matrices: exact-kNN neighbourhood composition (same-sample fraction, Shannon entropy, inverse Simpson, richness, dominance ratio vs global composition) at k = 15 and k = 50 in both spaces; global and within-sample kNN retention; coherence with the independent Phase 1 per-sample clusters; silhouette widths for a technical and a biological label; a sample-restricted-population over-correction check; 12 figures and 12 machine-readable outputs. Includes `--validation-mode` for synthetic smoke testing.
  - `scripts/shell/phase2/run_m12_evaluate.sh` — SLURM launcher (8 CPUs, 96G, 04:00:00, `ihc` / `ihc-grid-1-1-1`), sized from M11 accounting plus the two n x n distance matrices.
  - `reports/phase2/HARMONY_ASSESSMENT.md` — the M12 scientific assessment, sections A–J.
  - `reports/phase2/milestones/M12_REPORT.md`.
  - `results/phase2/harmony/evaluation/` — 12 outputs including `pre_post_mixing_summary.tsv`, `neighborhood_mixing_metrics_by_sample.tsv`, `technical_silhouette_summary.tsv`, `biological_preservation_summary.tsv`, `sample_restricted_population_check.tsv`, `biological_program_gene_sets.tsv`, `program_by_sample_composition.tsv`, `harmony_evaluation_parameters.tsv`, `m12_headline_metrics.json`, `prov_m12_evaluation.json`, `umap_harmony_m12_embedding.{tsv,rds}`.
  - `reports/phase2/figures/m12/` — 12 figures in PDF and PNG.
* **Changed**:
  - `reports/FIGURE_INDEX.tsv` — schema extended with `phase`, `milestone` and `slurm_job_id`; all 260 Phase 1 rows back-filled as `phase1` / `NA_backfilled_phase1`; 24 M12 rows appended (284 total). `snakemake_rule` retained for backwards compatibility, with Phase 2 rows carrying `NA_no_snakemake_in_phase2`.
  - `PROGRESS.md`, `CHANGELOG.md`.
* **Methodological decisions**:
  - **Fair comparison enforced**: identical cells (all 19,716, no subsampling), identical 30 dimensions, identical metadata, identical seed 42, and a post-Harmony UMAP built with the same `RunUMAP` defaults as the frozen Phase 1 UMAP with only the input reduction changed. No UMAP or neighbour parameter was optimised for appearance.
  - **Harmony was neither re-run nor re-tuned.** The M11 default result was evaluated as produced.
  - **`lisi`, `kBET` and `clustree` remain uninstalled.** Inverse Simpson was computed natively on exact kNN (the same estimator `lisi` uses) and a dominance ratio against global composition replaced kBET. Nothing was installed, upgraded or removed.
  - **Composition-adjusted references** are used throughout: with unequal sample sizes, perfect mixing means a same-sample fraction of 0.3065, entropy of 1.2683 nats and inverse Simpson of 3.2627 — not 0, log K and 4.
  - **Legacy `orig.anno` was deliberately excluded** as legacy integration-derived annotation. Two non-legacy biological references were used instead: the 54 independently derived Phase 1 per-sample clusters, and ten canonical broad-lineage program scores (a documented sanity check, explicitly not an annotation, never written into any saved object).
  - The Harmony UMAP is persisted as a 2-column embedding only; no third multi-gigabyte object was created.
* **Results**:
  - **Pipeline validated**: M12 PRE-Harmony per-sample metrics reproduce Phase 1 M8 to 15 significant figures. Phase 1's headline 0.9613 is the unweighted mean of four per-sample means (cell-weighted 0.9728); Phase 1's entropy 0.0966 is in log2 (0.0490 nats).
  - **Technical mixing improved**: same-sample neighbour fraction 0.9728 → 0.7784 at k = 15 (29.2% of the achievable gap) and 0.9443 → 0.7053 at k = 50 (37.5%); entropy 0.0490 → 0.3565 nats; inverse Simpson 1.0549 → 1.4408; dominance ratio 3.30 → 2.53. Every sample improved.
  - **Mixing is compartment-specific — the decisive result**: Panleukocyte 0.933 → 0.527, Myeloid 0.958 → 0.611, T/NK 0.953 → 0.689, Endothelial 0.962 → 0.709, Mural 0.951 → 0.724, versus Schwann/neural-crest 0.990 → 0.900 and B/plasma 0.986 → 0.918. Harmony aligned the shared immune and vascular compartments while leaving the presumptive malignant Schwann-lineage compartment (86.5% MPNST_1) patient-private, which is the desired behaviour and explains the modest global figure.
  - **Biology preserved**: within-sample kNN retention 0.8206; coherence with the integration-free Phase 1 clusters 0.8644 → 0.8427 (−2.5% relative); overall canonical-program silhouette 0.0627 → 0.0729 (improved), with six of ten programs more cohesive.
  - **Over-correction present but bounded**: B/plasma cohesion −44% (0.437 → 0.243; 93.5% MPNST_3) and fibroblast cohesion −44% (0.174 → 0.097; 70.6% MPNST_4); MPNST_4 shows the weakest within-sample preservation (0.781).
  - **Under-correction real but largely appropriate**: only 29.2% of the mixing gap closed and dominance ratios remain above 2.2 (highest for the small samples, MPNST_2 4.35 and MPNST_3 4.95, reflecting uniform default `theta = 2` under a 3.3x size imbalance) — but the residual sits in compartments where patient-private structure is expected.
  - **A metric that failed**: global silhouette of `sample_id` was −0.0104 *pre*-Harmony despite 97% same-sample neighbours, because each sample spans malignant, immune and stromal states so within- and between-sample mean distances nearly cancel. Recorded as a methodological caveat; it contributes nothing to the conclusion.
* **Recommendation**: **`ACCEPT DEFAULT HARMONY WITH CAVEATS`**. Two caveats carried into M13–M17: (C1) any downstream B/plasma cluster may be distorted and must be cross-checked against the preserved non-integrated baseline before annotation; (C2) the same for fibroblast/stromal clusters and MPNST_4-derived structure. A sensitivity analysis is **not** recommended and was not run — raising `theta` would press hardest on the already-fragile Schwann/neural-crest and B/plasma populations, and lowering it would undo the immune/vascular alignment that justifies integrating; two candidate experiments are documented in `HARMONY_ASSESSMENT.md` §J for the researcher to authorise if wanted.
* **Confounding**: `sample_id` remains simultaneously dataset, patient and the only batch proxy. Residual structure cannot be apportioned between uncorrected batch effect and real between-patient tumour biology; the compartment-resolved result is a consistency argument, not proof. The integrated embedding cannot support any between-tumour, between-condition or differential-abundance claim.
* **HPC execution**:
  - SLURM job `19886628` (`p2_m12_evaluate`, ihc-grid-1-1-1, 8 CPUs, 96G, 04:00:00): **COMPLETED**, ExitCode 0:0, Elapsed 00:04:09, MaxRSS 8,338M (8.14 GiB, 8.5% of request), TotalCPU 00:03:35. Zero R warnings, zero R errors; stderr holds one informational Seurat message about the `RunUMAP` backend. Validated first on an 800-cell synthetic structural mirror carrying a real Harmony reduction; no real object was loaded on the login node.
* **STOP**: Phase 2 halted at the M12 gate pending researcher authorization for M13.

[M13]

Post-Harmony Neighbours, UMAP and Clustering Sweep

*2026-09-02 — continuous execution M12→M17 authorized. Primary resolution 1.0 (26 clusters), alternative 0.7 (21 clusters).*

* **Added**:
  - `scripts/R/phase2/cluster_sweep_harmony.R` — builds the Harmony neighbour graph and UMAP, sweeps Louvain resolutions 0.1–1.0, preserves every solution, computes per-resolution diagnostics (sizes, tiny-cluster and sample-dominance flags, silhouette over a single shared distance matrix, canonical-program coverage and purity, adjacent-resolution ARI/NMI), applies a documented composite selection rule, carries the M12 caveats forward as machine-readable per-cluster flags, writes 22 figures and 10 tables, and validates the saved object by disk round-trip.
  - `scripts/shell/phase2/run_m13_clustering.sh` (8 CPUs, 128G, 08:00:00).
  - `scripts/R/phase2/m12_supplement_figures.R` and `scripts/shell/phase2/run_m12_supplement.sh` — the individually-named full-page M12 pre/post panels requested in the completion authorization.
  - `reports/phase2/CLUSTERING_ASSESSMENT.md`, `reports/phase2/milestones/M13_REPORT.md`.
  - `results/phase2/clustering/phase2_harmony_clustered.rds` (5.64 GB, md5 `b91f0eecd73e202804ac7c6859db4c46`) plus 10 tables; `results/phase2/figures/M13/` (44 files); `results/phase2/figures/M12/` (16 files + a not-applicable note).
* **Changed**: `PROGRESS.md`, `CHANGELOG.md`.
* **Scientific decisions**:
  - Neighbours and clustering use `k.param = 20` and Louvain, matching Phase 1 `config/config.yaml`, so the post-integration sweep is comparable to the Phase 1 per-sample sweeps.
  - All ten resolutions preserved as `postint_harmony_clusters_res_*`; nothing overwritten.
  - Resolution selected by a coded composite rule with fixed weights; **UMAP appearance is explicitly not a criterion**.
  - Resolution 1.0 sits at the sweep boundary, so its stability term is one-sided and inflated — stated as a caveat, with the interior, two-sided-validated resolution 0.7 offered as the conservative alternative.
  - `FindClusters()` overwrites the legacy `seurat_clusters` column and the active identities as a side effect; both are snapshotted before the sweep and restored verbatim, and a per-column digest check proves every one of the 93 original metadata columns is unchanged in value and R type.
* **Results**: 13 → 26 clusters across the sweep; **no tiny cluster at any resolution** (smallest 141 cells); every adjacent-resolution ARI ≥ 0.89; 11 of 26 primary clusters >60% one sample (expected, per M12); 3 clusters flagged for the M12 B/plasma caveat and 6 for the fibroblast caveat; strongest Schwann-lineage/tumour-like candidates C9, C14, C8 (all MPNST_1-dominated).
* **HPC execution**: three documented failures — 19886683 (preservation guard fired; the whole-frame metadata digest was replaced with a stricter per-column check), 19886685 (unnamed-vector indexing), 19886687 (`sprintf("%d", median())`) — then **19886690 COMPLETED**, ExitCode 0:0, Elapsed 00:09:37, MaxRSS 12,679,688K (12.09 GiB), zero R warnings. All failures were code defects fixed at source; no resource was increased in response to any failure. M12 supplement: 19886675 FAILED (`dpi = NA` rejected by ggplot2 4.0.1), 19886682 COMPLETED 00:01:32.

[M14]

Marker Discovery

*2026-09-02 — cluster-characterisation markers on the M13 primary clustering (resolution 1.0, 26 clusters).*

* **Added**:
  - `scripts/R/phase2/discover_markers_harmony.R` — inspects Seurat v5 assay/layer structure explicitly, runs `PrepSCTFindMarkers()` when the SCT assay carries more than one model, runs `FindAllMarkers()` on both the primary and the alternative resolution, writes complete/filtered/top-10/20/50 marker tables plus a per-cluster summary, and produces the marker figure suite.
  - `scripts/shell/phase2/run_m14_markers.sh` (12 CPUs, 250G, 12:00:00).
  - `reports/phase2/MARKER_REPORT.md`, `reports/phase2/milestones/M14_REPORT.md`.
  - `results/phase2/markers/` — 11 tables; `results/phase2/figures/M14/` — 6 figures × PDF + PNG.
* **Changed**: `PROGRESS.md`, `CHANGELOG.md`.
* **Scientific decisions**:
  - Seurat v5 layer handling was **inspected, not assumed**: the `RNA` assay's counts/data layers are split four ways and the `SCT` assay carries **four SCTransform models**, so `PrepSCTFindMarkers()` is mandatory and was applied before any test.
  - Testing used the `SCT` `data` layer. **Harmony coordinates were never used for marker testing** — Harmony is an embedding-level method and carries no expression values.
  - Marker parameters match Phase 1 `config/config.yaml: markers` (`wilcox`, `min.pct = 0.25`, `logfc.threshold = 0.25`, `only.pos = TRUE`), keeping M14 comparable to the Phase 1 per-sample runs.
  - Every marker table and report carries an explicit statement that these are cluster-characterisation markers and **not** condition-level differential expression.
* **Results**: 35,437 marker rows over 26 clusters, 31,774 significant; every cluster returned markers (212–2,901). Key evidence for M15: C9 retains MPZ with GFRA3/ABCB5; C8 expresses L1CAM; C14 expresses SHH with GAL3ST1/KLK6; C21 markers (LILRA4/CLEC4C/SPIB/GZMB) are decisive for pDC and override the coarse M12 programme score; C15 is mitochondrial-dominated (technical) and C23 ribosomal-pseudogene-dominated; C12 is genuinely ambiguous between a hypoxia programme and the reported perineurial GLUT1/ITGB4 signature.
* **HPC execution**: SLURM job `19886699` COMPLETED, ExitCode 0:0, Elapsed 00:04:26, MaxRSS 17,193,456K (16.40 GiB, 6.6% of the 250G request), zero R warnings, no retries.

[M15]

Literature-Grounded Cell-Type Annotation

*2026-09-02 — hierarchical, conservative, CNV-free annotation of the 26 M13 primary clusters.*

* **Added**:
  - `config/phase2/annotation_map_M15.tsv` — the auditable annotation map. Every cluster's Level 1/2/3 label, positive markers, negative markers, supporting pathways, conflicting evidence, confidence, source and DOI/URL live in one reviewable TSV.
  - `scripts/R/phase2/annotate_celltypes.R` — applies the map, builds `ANNOTATION_EVIDENCE.tsv` (merging the curated evidence with data-derived top markers and per-sample counts), and produces the annotation figure suite.
  - `scripts/shell/phase2/run_m15_annotation.sh` (8 CPUs, 128G, 06:00:00).
  - `reports/phase2/ANNOTATION_REPORT.md`, `reports/phase2/milestones/M15_REPORT.md`.
  - `results/phase2/annotation/phase2_harmony_annotated.rds` plus `ANNOTATION_EVIDENCE.tsv`, `cluster_to_annotation_map.tsv`, records and provenance; `results/phase2/figures/M15/` (PDF + PNG).
* **Changed**: `PROGRESS.md`, `CHANGELOG.md`.
* **Scientific decisions**:
  - **Annotation is data, not code.** The script applies the map and infers no label, keeping scientific judgement visible and reviewable.
  - **No CNV inference** — outside the authorised Phase 2 scope. Malignant calls integrate Schwann/neural-crest lineage markers, MPNST-specific literature, absence of convincing immune/stromal/endothelial identity, and sample provenance. Provenance is never used alone and the inference "non-immune ⇒ tumour" is never made.
  - **Graded, conservative malignant labels**: MPNST-like malignant (SCP-like) 1,680 cells — C9 retains MPZ with GFRA3/CRYAB/NOV/ABCB5 and C8 expresses L1CAM; MPNST-like malignant (NC-like) 440 — C14 expresses SHH with GAL3ST1/KLK6; Schwann-lineage tumour-like 1,011; Candidate malignant 1,231 at Low confidence; Cycling tumour-like 289.
  - **The fibroblast-vs-Mes-NC-like ambiguity is stated rather than hidden**: the malignant fraction is 23.6% conservatively and could be ~49% if the four fibroblast clusters are Mes-NC-like malignant. CDKN2A expression argues against malignancy in C25; the PI16+ epineurial signature argues for genuine fibroblast in C3.
  - **Uncertain labels retained, not forced**: C12 (hypoxia programme vs the reported perineurial GLUT1/ITGB4 signature) and C23 (ribosomal pseudogenes); C15 labelled technical (mitochondrial-high), matching the Phase 1 MPNST_4 flag.
  - **Marker evidence overrode a programme-score prior**: C21 was called B/plasma by the coarse M12 score because pDCs share IGJ and MZB1; LILRA4/CLEC4C/SPIB/GZMB are decisive for plasmacytoid dendritic cells and the annotation follows the markers.
  - Ten primary sources with DOIs/PMIDs were triangulated (MPNST spatial transcriptomics, MPNST single-cell multiomics, NF1-PNST immunotyping, peripheral-nerve fibroblast/perineurial biology, canonical immunology). No single paper was copied.
* **Results**: Level 1 — Fibroblast/Stromal 5,499 (27.9%), Malignant/tumour 4,651 (23.6%), Myeloid 4,154 (21.1%), B/Plasma 1,942 (9.8%), T/NK 1,352 (6.9%), Endothelial 960 (4.9%), Uncertain 720 (3.7%), Other/technical 438 (2.2%). 17 Level 2 cell types. Confidence High 51.3% / Moderate 38.8% / Low 9.9%.
* **HPC execution**: SLURM job `19886713` COMPLETED, ExitCode 0:0, Elapsed 00:08:35, MaxRSS 8,963,852K (8.55 GiB, 6.7% of request), zero R warnings, no retries; all 7 round-trip validation checks TRUE.

[M16]

Annotation Refinement and Composition

*2026-09-02 — computational annotation review plus descriptive composition. Zero labels changed.*

* **Added**:
  - `scripts/R/phase2/refine_and_compose.R` — scores six canonical compartment panels per cell, averages them per cluster, compares each cluster's assigned Level 1 with its own panel ranking under four documented rules, writes `*_initial` and `*_refined` annotation columns plus the review rule/outcome, and builds the descriptive composition tables and figures.
  - `scripts/shell/phase2/run_m16_composition.sh` (8 CPUs, 128G, 06:00:00).
  - `reports/phase2/milestones/M16_REPORT.md`.
  - `results/phase2/composition/phase2_harmony_refined.rds` plus 12 tables; `results/phase2/figures/M16/` — 10 figures × PDF + PNG.
* **Changed**: `PROGRESS.md`, `CHANGELOG.md`.
* **Scientific decisions**:
  - **The panel review adjusts confidence, never labels.** M15 labels rest on specific marker genes and cited literature, which is stronger evidence than a coarse six-panel mean; disagreement therefore lowers confidence and is recorded, rather than overriding a literature-grounded call. Zero labels were changed.
  - **Previous annotation preserved, not replaced**: `postint_celltype_level1/2/3_initial` and `postint_annotation_confidence_initial` hold the M15 values alongside the `*_refined` versions, with `postint_annotation_review_rule` and `postint_annotation_review_outcome` recording the reason.
  - **Composition is descriptive only** — no inferential condition-level or differential-abundance test was performed, and every composition table carries that statement in a `note` column. Cells are not independent biological replicates.
  - `sample_id` is simultaneously dataset and patient, so those three composition groupings are one table; all three filenames are written for convenience and are identical by construction. No condition table was fabricated — a `cell_proportions_by_condition_NOT_APPLICABLE.txt` note records that no such field exists.
* **Results**: 18 clusters confirmed, 5 confidence-downgraded (C7, C13, C18, C20 malignant calls unsupported by the classical Schwann panel — consistent with reported MPNST dedifferentiation but honestly downgraded; C21 pDC downgraded as a conservative artefact of IGJ/MZB1 panel overlap), 3 deliberately unscored. Confidence after refinement: High 9,828 (49.8%), Moderate 6,637 (33.7%), Low 3,251 (16.5%). **4 of 17 cell types are driven >80% by a single sample/patient** and are flagged as such.
* **HPC execution**: SLURM job `19886728` COMPLETED, ExitCode 0:0, Elapsed 00:07:34, MaxRSS 8,890,820K (8.48 GiB, 6.6% of request), zero R warnings, no retries; all 6 round-trip checks TRUE.

[M17]

Phase 2 Validation, Freeze and Handoff

*2026-09-03 — PHASE 2 COMPLETE AND FROZEN.*

* **Added**:
  - `scripts/R/phase2/validate_and_freeze.R` — runs 26 pre-save validation checks, verifies Phase 1 immutability by md5 and mtime, writes the deterministic final object, **reloads it from disk** and runs 14 post-reload checks, assembles the final table suite by copy, and builds the Phase 2 manifest.
  - `scripts/shell/phase2/run_m17_freeze.sh` (8 CPUs, 128G, 06:00:00).
  - `results/phase2/phase2_final_object.rds` — 6,055,929,223 bytes, md5 `63146e84e43d8f036cca6fd6d1ef99b3`, sha256 `51f0833f9b1e76a83c8a94046fd8bffcf647d842c59dae8f7b0705b0e9d9d233`, 19,716 cells.
  - `results/phase2/phase2_manifest.json`, `results/phase2/tables/final/` (17 files), `results/phase2/figures/final/` (28 files + a not-applicable note).
  - `reports/phase2/PHASE2_HANDOFF.md` (19 sections, self-contained), `reports/phase2/milestones/M17_REPORT.md`.
* **Changed**:
  - `reports/FIGURE_INDEX.tsv` — rebuilt to 438 rows covering 260 back-filled Phase 1 figures and all 178 Phase 2 figure files across M12–M16 and the final suite.
  - `PROGRESS.md` status advanced to **PHASE 2 = COMPLETE**; `CHANGELOG.md`.
* **Scientific decisions**:
  - **No new exploratory biology in M17** — the milestone validates, freezes and documents.
  - **`saveRDS()` success is not treated as proof of validity**: the final object is written, memory freed, re-read from disk and re-verified against in-memory digests.
  - **Scientific-scope guards are machine-checked**: the final object is asserted to contain no column matching `cnv|infercnv|copykat`, `pseudobulk|deseq|edger|condition_de` or `pseudotime|velocity|monocle|slingshot`, proving no prohibited Phase 3 analysis leaked into Phase 2.
  - **Phase 1 immutability re-verified at the freeze** by md5 and mtime on `processed_mpnst.rds` and the Phase 1 handoff object.
  - Requested condition figures and tables were **not fabricated**; each absence is recorded in a `*_NOT_APPLICABLE.txt` note explaining that no biological-condition field exists in the data.
* **Results**: all 26 pre-save and 14 post-reload checks TRUE. Phase 2 completion criteria all met. Total successful Phase 2 compute across nine jobs: **54 minutes**, peak memory **16.40 GiB** (3.6% of the 450G envelope). Four documented job failures, all code defects fixed at source, with no resource ever increased in response to a failure.
* **HPC execution**: SLURM job `19886743` COMPLETED, ExitCode 0:0, Elapsed 00:07:39, MaxRSS 8,520,016K (8.13 GiB, 6.4% of request), zero R warnings.
* **STOP**: **Phase 2 is frozen. Phase 3 was not begun and requires separate authorization.**

[AMENDMENT — CCC-Oriented Annotation Layer]

*2026-09-03 — Phase 2 expanded, after the M17 freeze, to add a cell–cell-communication-oriented annotation layer alongside the detailed annotation. Prior history is unchanged; this entry is appended.*

* **Authorized change**: preserve the detailed literature-supported annotation as the source of truth, and additionally create a simplified `annotation_ccc` layer in which the evidence-supported malignant MPNST populations are collapsed into a single `MPNST-Tumor` population while immune, stromal and endothelial populations retain biologically meaningful identities. The simplified layer is intended for downstream Phase 3 tumour–microenvironment communication analysis and **does not replace or destroy the detailed annotation**.
* **Rationale**: the detailed annotation resolves the malignant compartment into five states across eight clusters. That is scientifically correct but fragments the tumour compartment into populations too small and too patient-specific to act as a communication source or receiver, so a collapsed tumour identity is needed for CCC while immune/stromal/endothelial resolution must be retained to evaluate tumour↔immune, tumour↔stromal and tumour↔endothelial axes.
* **Scope discipline**: applied as a new milestone **M15A** consuming the validated M16 object, followed by a re-run of **M17**. **Harmony (M11), integration assessment (M12), clustering (M13) and marker discovery (M14) were NOT rerun** — their outputs are valid and were reused, per amendment §36.
* **Collapse rule — explicitly not "non-immune = tumour"**: applied to the detailed Level 2 label. `MPNST-like malignant (SCP-like)`, `MPNST-like malignant (NC-like)`, `Schwann-lineage tumour-like` and `Cycling tumour-like` collapse into `MPNST-Tumor`. `Candidate malignant` does **not** — it is retained separately as `Candidate-Malignant-Unresolved` because the label is provisional, its confidence is Low, and its only evidence is patient-private neural expression. Fibroblast, endothelial, pericyte, uncertain and low-quality populations are never absorbed. Every decision including every exclusion is recorded with rationale in `results/phase2/annotation/CCC_ANNOTATION_MAPPING.tsv`.
* **Tumour heterogeneity preserved**: the detailed malignant state labels remain in the final object alongside `annotation_ccc`, so tumour-state-specific communication analysis remains possible later without re-running annotation.
* **Documentation updated first**: `PROJECT.md` gained a full Phase 2 architecture section (P2.1–P2.7) describing the hierarchical annotation strategy, the `MPNST-Tumor` definition and prohibition, tumour-heterogeneity preservation, the sample-aware statistical rule for Phase 3, the annotation metadata fields, and the prohibitions still in force. `PROGRESS.md` records the amendment with rationale, planned fields, collapse rule and expected outputs.

[M15A]

CCC-Oriented Annotation Layer

*2026-09-03 — `annotation_ccc` built from the detailed annotation; Phase 2 re-frozen.*

* **Added**:
  - `config/phase2/ccc_annotation_map.tsv` — the auditable collapse rule, one row per detailed Level 2 label, carrying `annotation_ccc`, `collapsed_to_mpnst_tumor`, `ccc_compartment`, `ccc_ready` and a written `mapping_rationale` for every decision **including every exclusion**.
  - `scripts/R/phase2/build_ccc_annotation.R` — applies the map, sub-resolves the single T/NK cluster by canonical marker gating, and writes the mapping table, summary, population-size audit, MPNST-Tumor composition, CCC marker summary, by-sample/by-patient composition tables and the 14-figure CCC suite. Validates that every detailed annotation column is byte-identical and that no prohibited population entered `MPNST-Tumor`.
  - `scripts/shell/phase2/run_m15a_ccc_annotation.sh` (8 CPUs, 128G, 06:00:00).
  - `results/phase2/annotation/phase2_harmony_ccc.rds`; `CCC_ANNOTATION_MAPPING.tsv`, `CCC_ANNOTATION_SUMMARY.tsv`, `CCC_POPULATION_SIZE_AUDIT.tsv`, `MPNST_TUMOR_COMPOSITION.tsv`, `CCC_MARKER_SUMMARY.tsv`, `m15a_ccc_annotation_record.json`, `prov_m15a_ccc.json`; `ccc_counts_by_{sample,patient}.tsv` and `ccc_proportions_by_{sample,patient}.tsv`; `results/phase2/figures/CCC_annotation/` (14 figures × PDF + PNG).
  - `reports/phase2/CCC_READINESS.md` — Phase 3 readiness audit, verdict **READY WITH CAVEATS**.
* **Changed**:
  - `scripts/R/phase2/validate_and_freeze.R` — now consumes the CCC object and adds eleven CCC-specific pre-save guards plus four post-reload guards (37 and 18 checks in total).
  - `results/phase2/phase2_final_object.rds` **re-frozen** — md5 `153d5f6acc70f9c05aa48cabc4f4ac2d`, sha256 `62e97524836309170c38be9335a10a4c76baa48a246b2963d2d440cebead0cf3` (was `63146e84…`), now carrying both annotation layers.
  - `results/phase2/phase2_manifest.json` — new `ccc_annotation` section with field names, mapping table and config, tumour-collapse criteria, per-population counts and proportions, MPNST-Tumor totals, readiness report, figures, script and record, plus three additional known limitations and the M15A/re-freeze JobIDs.
  - `results/phase2/tables/final/` 17 → **27 files**; `results/phase2/figures/final/` renumbered to the amendment's scheme (**34 files**); `reports/FIGURE_INDEX.tsv` 438 → **472 rows** with `figure_type`, `annotation_field`, `input_object` and `notes` columns.
  - `reports/phase2/PHASE2_HANDOFF.md` — new section 8A (CCC-Oriented Annotation) and refreshed final-object checksums.
  - `reports/phase2/milestones/M15_REPORT.md`, `M16_REPORT.md`, `M17_REPORT.md` — amendment addenda appended; prior content unchanged.
  - `PROJECT.md`, `PROGRESS.md`, `README.md`.
* **Results**: 17 CCC identities. **`MPNST-Tumor` = 3,420 cells (17.35%)** collapsed from clusters C7, C8, C9, C14, C20 (four detailed malignant states, all retained). Immune 7,448 (37.8%) across nine identities; stromal/endothelial 6,459 (32.8%); uncertain/excluded 2,389 (12.1%). The single T/NK cluster was sub-resolved into CD8-T 382, CD4-T 606, NK 106, T-cell-other 258. **14 of 17 identities are CCC-ready and all 14 are present in all four samples.**
* **Scientific decisions**: the prohibited inference `if (!immune) -> MPNST-Tumor` was never used, and M17 machine-verifies it; `Candidate malignant` (C13, C17, C18; 1,231 cells) was deliberately **not** collapsed and is retained as `Candidate-Malignant-Unresolved`; fibroblasts were **not** absorbed despite the unresolved Mes-NC-like hypothesis; the immune compartment was **not** collapsed into a single label; `NK` (106 cells) was flagged rather than merged; tumour heterogeneity is fully recoverable; **no marker rediscovery** — the M14 tables were reused and Harmony coordinates were never treated as expression.
* **HPC execution**: SLURM `19893067` (M15A) COMPLETED, ExitCode 0:0, Elapsed 00:07:33, MaxRSS 8,853,320K (8.44 GiB), 1 R warning (sparse CD4 detection softens the CD4-T boundary). SLURM `19893267` (M17 re-freeze) COMPLETED, ExitCode 0:0, Elapsed 00:07:59, MaxRSS 8,521,956K (8.13 GiB), 0 warnings, no check FALSE. No failures or retries.
* **STOP**: Phase 2 complete and re-frozen with both annotation layers. Phase 3 not begun; requires separate authorization.

[PHASE 3 — M18 to M27]

MPNST Tumour–Immune Cell–Cell Communication, Multi-Method Concordance, Receiver-State Analysis and LochNESS

*2026-09-03 — Phase 3 executed continuously M18→M27 and frozen.*

* **M18 — reconstruction and method audit**: audited the perturb-seq pipeline (33 components: 9 ADAPT, 11 REFERENCE ONLY, 13 NOT APPLICABLE; energy distance/MMD and ORA rejected on **scientific**, not technical, grounds) and the LochNESS implementation across three sources (publication, official MMCA R, our Python port). **MMCA is a repository of scripts, not an installable R package**, so `remotes::install_github()` was not attempted; it was vendored read-only to `external/MMCA_ref/` at commit `af629c49`. Six divergences identified between the official R implementation and our Python port; the two that matter are the **absent same-sample exclusion** and a fixed rather than `round(0.5·√N)` k. Compute-node network access verified before installing anything inside SLURM.
* **Environment**: installed liana 0.1.14, CellChat 2.2.0.9001, nichenetr 2.2.1.1, OmnipathR 3.14.0, decoupleR 2.12.0, NMF 0.28, FNN 1.1.4.1, systemfonts 1.3.2, ggraph 2.2.2, ggpubr 1.0.0; CellPhoneDB 5.0.1 in a **deliberately isolated** `cpdb_env`. **Seurat 5.4.0, SeuratObject 5.3.0, harmony 1.2.4, Matrix 1.7.4 and R 4.4.3 verified unchanged after every stage** (`r-base` explicitly pinned). One dependency change recorded honestly: **igraph 2.2.1 → 2.1.4**, downgraded by the conda solver; Phase 2 outputs are frozen artefacts and unaffected. Three install attempts were needed; each failure was root-caused (source compilation of `systemfonts`/`Deriv` cascading to svglite/car/rstatix/ggpubr; then `systemfonts` 1.2.3 < the 1.3.0 svglite requires; then nichenetr's `shadowtext`/`units`/`sf` system-library chain) and fixed with prebuilt conda-forge binaries rather than by weakening the environment.
* **M19 — sample-aware CCC input preparation**: `scripts/phase3/ccc/prepare_ccc_inputs.R`. Expression basis is the **RNA assay, joined + LogNormalize**, not SCT, because the SCT assay carries four models and is not on a common footing across samples; **Harmony coordinates are never used as expression**. Minimum-cell policy stated rather than silently chosen (≥10 cells per population per sample, CellChat's `min.cells` default), with the explicit rule that a population below the minimum is **NOT EVALUABLE, never "no signalling"**. Result: 14 populations, 17,327 cells, **all populations and all 182 directed pairs evaluable in all four samples**.
* **M20 — multi-method execution**: LIANA, CellChat and CellPhoneDB each run **per sample**. Three documented failures fixed at source: liana 0.1.14 calls the **defunct** `GetAssayData(slot=)` (fixed by routing through `SingleCellExperiment` rather than downgrading SeuratObject); CellPhoneDB 5.0.1's `download_database()` keyword is `cpdb_version`; and the serial CellChat job was **cancelled at 44% of sample 1** when it projected a ~5 h runtime and re-run as four concurrent per-sample jobs, finishing in 44 min with identical results.
* **M21 — concordance**: 502,846 distinct keys, 36,486 supported. **Raw scores never averaged** — each framework contributes a boolean flag. Every row carries `testable_*` alongside `supported_*`, because 97% of keys exist in only one framework's resource, so "single-method" usually means resource non-overlap rather than disagreement. High concordance 5,848 · Moderate 499 · Single-method 24,654 · **Discordant 5,485 (preserved, not averaged away)**. Among the 1,347 jointly testable interactions, 763 (57%) were supported by all three. **LIANA/CellPhoneDB partial non-independence declared.**
* **M22 — prioritisation**: lexicographic tiers on concordance, sample recurrence and expression support; **no weighted composite invented**. 12 prioritised axes curated with DOIs in `CCC_LITERATURE_EVIDENCE.tsv`; **no novelty claimed from a failed search**.
* **M23 — receiver response (NicheNet)**: gene set of interest = receiver marker genes (cluster characterisation, **not** condition DE). Myeloid receivers converge on **CSF1**, lymphoid on **IL15**, endothelial on TGFB1/VEGFA. **Fibroblast is a clear negative (best AUPR 0.022).** Reported divergence: **APP is not among the top NicheNet ligands for macrophages** despite being the strongest LR finding — bounding that claim to "predicted engagement", not "drives the macrophage state".
* **M24 — LochNESS**: design written **before** execution. Official MMCA formulation in R with same-sample exclusion, `k = round(0.5·√N)`, L2 Harmony space, exact `FNN::get.knnx`, and two nulls. **Implementation comparison: R vs Python with mathematically equivalent inputs gives Pearson 1.0000 and max abs diff 0**; vs the perturb-seq formulation as written, 0.3764, because that formulation omits same-sample exclusion. **Biological result is a negative**: no lineage shows receiver-state structure associated with the tumour-derived APP context beyond chance (best descriptive p = 0.333 against a hard 0.167 floor at n = 4), and the score is unstable across embeddings (ρ 0.02–0.42) and reference samples (ρ 0.31–0.50). LochNESS was never treated as a ligand-receptor method.
* **M25 — integration**: an **evidence matrix, not a composite score**. 325 interactions reach all four streams; 547 tumour-centric reach ≥3. LochNESS enters only as receiver-lineage context.
* **M26 — robustness**: leave-one-patient-out retention 71–79%, so no single patient dominates; but **53.6% of supported interactions rest on one patient** and are flagged, not promoted. Rare populations have a median of 1–2 supporting patients and were flagged rather than merged.
* **M27 — freeze**: `results/phase3/phase3_manifest.json` (31 sections, 11 caveats, every JobID including every failure and root cause), `reports/phase3/PHASE3_HANDOFF.md` (26 sections), 36 final figures, 29 final tables, 10 reports, 10 milestone reports, `reports/FIGURE_INDEX.tsv` extended to 578 rows.
* **Scientific findings**: the MPNST tumour compartment is predicted to act principally on **myeloid cells and the vasculature** — a reproducible myeloid-directed set (APP→CD74 across five receivers, CD99→PILRA, ANXA1→FPR1, HLA-F→LILRB1/2) with the myeloid programme independently best explained by **CSF1**; canonical angiogenesis (VEGFA→KDR/FLT1/NRP1); and a **reciprocal perivascular Notch circuit**, in a tumour type where Notch is already implicated in Schwann-cell transformation. The lymphoid picture is **mixed, not uniformly suppressive**.
* **Prohibitions respected**: no survival analysis, treatment-response modelling, deep learning, Transformers, trajectory, RNA velocity, CNV inference, spatial inference, unrelated pseudobulk condition DE, large unrelated pathway screens, or arbitrary composite scores. No Snakemake.
* **STOP**: Phase 3 complete and frozen. Phase 4 (spatial validation) not begun; requires separate authorization.

[PHASE 4 — M28 to M35]

MPNST SCEVAN Malignancy Refinement, Tumour-State Resolution and Targeted CCC Reassessment

*2026-09-03 — Phase 4 authorized and begun.*

* **Scope**: resolve the malignant compartment using inferred copy number, then reassess only
  the tumour-dependent Phase 3 conclusions. Phase 2's conservative `MPNST-Tumor` definition
  (3,420 cells, 17.35%) was correct given marker evidence alone, but leaves 5,064 `Fibroblast`,
  1,231 `Candidate-Malignant-Unresolved`, 720 `Uncertain` and 435 `Pericyte-VSMC` cells
  unresolved — and that uncertainty propagates into every tumour-centric interaction.
* **Method**: SCEVAN (De Falco *et al.*, Nat Commun 14:1074, 2023; doi:10.1038/s41467-023-36790-9).
  Chosen for raw-count input, automatic confident-normal detection, large-scale CNA profiles,
  subclonal capability, multi-sample clonal comparison and low computational burden in a native
  R workflow. **Not claimed to be universally superior to inferCNV**; no second CNV method is
  added for benchmarking.
* **Environment — no change**: SCEVAN 1.0.3 and yaGST 2017.8.25 were **already installed** in
  `R_env`, so no installation was attempted and no dependency moved. R 4.4.3, Seurat 5.4.0,
  SeuratObject 5.3.0, harmony 1.2.4, Matrix 1.7.4 remain exactly as in Phase 2 and Phase 3.
  Pre-Phase-4 environment captured to `reports/phase4/environment/R_env_PRE_PHASE4.yml`.
* **Anti-circularity rule adopted before any run**: the primary SCEVAN run passes
  `norm_cell = NULL` so confident-normal detection sees **no Phase 2 label at all**. Reading
  the SCEVAN source established that `FIXED_NORMAL_CELLS = TRUE` executes
  `cellType_pred[!cellType_pred %in% norm_cell_names] <- "malignant"` — every non-reference
  cell forced to malignant. That is the non-immune-equals-tumour inference Phase 2 banned,
  arriving through a function argument, so **`FIXED_NORMAL_CELLS = TRUE` is prohibited in this
  project**. `Fibroblast`, `Pericyte-VSMC`, `Candidate-Malignant-Unresolved` and `Uncertain`
  are never used as fixed normal references — they are the question.
* **Documentation updated before production jobs** (§4): `PROJECT.md` §P4.1–P4.12,
  `README.md`, `PROGRESS.md`, `CHANGELOG.md`.
* **New directories**: `scripts/phase4/{scevan,malignancy,tumor_states,ccc_refinement,utils}`,
  `scripts/shell/phase4/`, `results/phase4/{scevan/by_sample,malignancy,tumor_states,ccc_refinement,figures/final,tables/final}`,
  `reports/phase4/{milestones,environment}`, `logs/phase4/`.

* **M28 — reconstruction and SCEVAN feasibility** (2026-09-03; SLURM 19896576 F, 19896623 F,
  19896654 F, 19896672, **19896711**): reconstructed Phase 2/3 from the repository rather than
  from memory, then established the exact SCEVAN input basis. **19,716 cells; RNA assay, four
  per-sample `counts` layers, all integer-verified; gene identifiers are symbols; 28,340 of
  31,764 features (89.2%) map to SCEVAN's chr1–22 annotation.** Cell IDs round-trip 1:1 with
  zero duplicates. Per-sample sparse counts written to disk so the four M29 jobs never reload
  the 6 GB Phase 2 object, and `as.matrix()` was never applied to the full dataset — SCEVAN's
  own internal densify was measured instead (largest 1.71 GiB) and used to size M29 at 64 G
  rather than the 450 G maximum. The Phase 2 default `SCT` assay is deliberately unused: it
  carries four SCTransform models, the same reason Phase 3 chose RNA + LogNormalize.
* **M28 defects found and fixed without weakening anything**: (1) SCEVAN 1.0.3 accepts
  `output_dir` but hardcodes `path = "./output"` in `getScevanCNV`, `getScevanCNVfinal`,
  `plotAllClonalCN`, `plotAllSubclonalCN`, `plotConsensusCNA` and `analyzeSegm2`, and
  `plotCNclonal()` does not forward it — classification succeeded and then plotting died with
  `cannot open the connection`; fixed by giving each run its own working directory so the
  parameterised and hardcoded paths coincide, **without patching the package**. (2)
  `subcloneAnalysisPipeline` → `plotTSNE` requires **Python** umap-learn via reticulate, and
  reading the source showed `plotTSNE` runs *before* the line that writes subclone labels into
  `classDf` — so it is not an optional cosmetic step and skipping it would lose the clone
  assignments; fixed by installing umap-learn 0.5.12 into an **isolated `p4_umap_env`** reached
  through `RETICULATE_PYTHON`, with every `R_env` version printed before and after to prove
  nothing moved. (3) Our own defect: `set -u` versus conda's `qt-main` activation hook
  (`QT_XCB_GL_INTEGRATION: unbound variable`), fixed by dropping `-u` per conda's guidance.
  **No resource request was raised in response to any failure**; 96 G was set from the Phase 3
  measured load and peaked at 26.90 GiB.
* **M28 anti-circularity decisions, recorded before any result existed**: the primary SCEVAN
  run passes `norm_cell = NULL`; the sensitivity run supplies only the 7,448 high-confidence
  immune cells with `FIXED_NORMAL_CELLS = FALSE`; the 11,308 disputed cells are never a fixed
  reference; **Endothelial is also excluded from the reference set** so its malignancy call
  stays independent of the Phase 3 VEGFA and JAG/NOTCH findings that concern it. And
  `reports/phase4/MALIGNANCY_DECISION_RULES.md` — 20 ordered rules with every threshold fixed
  a priori — was written and committed **before the first SCEVAN classification was inspected**,
  so no rule could be tuned to the answer.
* **M28 limitations recorded**: SCEVAN's annotation covers chromosomes 1–22, so X/Y events are
  not assessable; and SCEVAN removes cell-cycle genes **and all `HLA-*` genes** before
  inference, so the CNV analysis is structurally blind to the HLA-E/HLA-F loci Phase 3
  highlighted. Both stated now rather than discovered later.
* **Smoke test** (technical only, never propagated): 600 random MPNST_2 cells through the full
  pipeline — 363 tumour / 214 normal / 23 filtered, 5 subclones, 74 output files, 3.86 min.
* **M29 — per-patient SCEVAN** (2026-09-03; SLURM array **19896712_[1-4]**, all COMPLETED, no
  failures and no retries): SCEVAN run **independently for each patient** because
  `sample_id` = patient = dataset, and **twice per patient** — primary with `norm_cell = NULL`
  (non-circular by construction) and sensitivity with the 7,448 high-confidence immune cells
  and `FIXED_NORMAL_CELLS = FALSE`. **18,449 of 19,716 cells assessed** (1,267 SCEVAN-filtered,
  recorded as *not assessed*, never as non-malignant); **8,794 SCEVAN-malignant**; 22 subclones.
  Runtime 8–42 min per patient; MaxRSS 6.2–**60.9 GiB against a 64 G request** — the M28-derived
  estimate was closer than intended and that is recorded rather than hidden. Primary/sensitivity
  agreement **0.9951 / 0.9663 / 0.0864 / 0.9974**: running both strategies is the only reason the
  one failure below was detectable at all.
* **M29 — the fibroblast result, and two sanity failures reported as such.** Malignant fraction
  of assessed cells for Phase 2 `Fibroblast`: **0.628 (MPNST_1), 0.559 (MPNST_2), 1.000
  (MPNST_4)** — replicated in three patients. Clone composition is the most direct evidence:
  **all four MPNST_2 subclones and seven of eight MPNST_4 subclones are fibroblast-dominated**.
  Conversely SCEVAN corroborates only **36.7%** of the cells Phase 2 called `MPNST-Tumor`, which
  is reported rather than buried — §65C anticipated copy-number-quiet malignant cells, which is
  why those cells became `Ambiguous` rather than `Non-malignant`. **MPNST_3 failed an immune
  sanity check outright**: its primary run calls the lymphoid compartment malignant and the
  myeloid compartment normal while the sensitivity run inverts it, 6 of 9 immune populations
  individually exceed 25% malignant, only 25 confident normal cells were found, it has the
  lowest depth of the four, and its three "clones" are pure immune lineages. **MPNST_4's high
  overall immune rate (0.293) is a different phenomenon** — 98.2% of its plasma cells, an
  immunoglobulin-locus (IGH 14q32 / IGK 2p11 / IGL 22q11) artefact documented for this method
  family, with 0 of 9 populations failing and the highest agreement of the four.
* **M29 warning**: every run emitted a `plotCloneTree` error caught by SCEVAN's own `tryCatch` —
  ggtree calls `ggplot2:::is.waive()`, removed in ggplot2 4.x. Only the clone **phylogeny plot**
  is lost; all clone assignments, `.seg` profiles, CNA matrices and onco-heatmaps are intact.
  **ggplot2 was deliberately not downgraded**, because that would destabilise the Phase 3 figure
  suite to gain a dendrogram that is not a Phase 4 deliverable.
* **M30 — malignancy integration** (2026-09-03; SLURM 19897096 superseded, **19897148**
  reported, 00:01:00, MaxRSS 1.38 GiB): **Phase 2 3,420 (17.35%) → SCEVAN 8,794 of 18,449
  assessed (47.7%) → refined 6,434 (32.63%)**, of which 3,261 High confidence; plus 9,078
  Non-malignant, **3,766 Ambiguous (19.10%)** and 438 Excluded. **R99 fail-safe = 0 cells**, so
  rule coverage is complete. The answer to "was 17.35% too conservative?" is **yes, but
  two-directionally**: Phase 4 adds 4,036 fibroblasts, 836 candidate-malignant, 92 uncertain and
  65 pericytes while withdrawing confidence from 2,015 of the 3,420 cells Phase 2 called
  malignant. Phase 2 did not so much undercount as look in the wrong place.
* **M30 — Amendments A1 and A2, both post-hoc and both labelled as such.** A1 added a per-sample
  immune sanity gate that disables malignant **promotions** from a failing sample, deliberately
  asymmetric (an unreliable CNV run may still support a `Non-malignant` call, because those
  always require concordant lineage evidence too). A2 was then written **after observing that
  the flat gate excluded MPNST_4, the sample carrying the strongest evidence** — the rules
  document states that ordering explicitly, because it is exactly the circumstance in which a
  reader should be suspicious. A2 excludes plasma cells from the sanity denominator on
  mechanistic grounds, adds a **breadth** criterion so a global inversion cannot dilute past the
  gate (which independently condemns MPNST_3), and adds rule **R7p** so 1,125 plasma cells are
  called `Non-malignant` with the artefact named rather than dropped into `Ambiguous`. Both
  gates are published for all four samples with a `flat_gate_would_have_failed` column so the
  stricter version can be applied by anyone who prefers it.
* **M30 — threshold sensitivity published, not discovered later.** The a priori
  `pop_frac_low = 0.25` was recomputed across 0.15–0.40. The **fibroblast finding is completely
  threshold-independent** (4,036 malignant at every grid point), while **`MPNST-Tumor`
  retention is entirely threshold-dependent** (3,266 versus 1,405 across the 0.20/0.25 boundary,
  because MPNST_1's pop_frac is 0.222). The latter is named as the most fragile number in
  Phase 4.
* **M30 — recurrent broad CNA events** in the three reliable patients: **chr18 loss 3/3, chr2
  gain 3/3, chr7 gain 3/3**, with chr3/chr5/chr6/chr11/chr15/chr19/chr22 events in 2/3. MPNST_4
  carries a 38.1 Mb CN=1 segment on chr17 whose **interval contains** the NF1 locus (17q11.2),
  and MPNST_1 and MPNST_4 both carry a 34.2 Mb CN=1 segment on chr22 containing NF2 (22q12).
  Reported as segmental events containing those loci, **not** as gene-level deletions (§51).
* **M30 — other decisions**: marker evidence may only **reduce** confidence, never raise it,
  since letting the markers Phase 2 already used raise confidence in a Phase 4 call would
  reintroduce the circularity §16 forbids; `Ambiguous` is a terminal outcome and never becomes
  tumour; and MPNST_4's endothelium (52.7% malignant, forming its own subclone) was **flagged**
  `Ambiguous` by rule R11 rather than promoted, endothelium having been kept out of the normal
  reference set precisely so that call would be independent.
* **Process note**: two duplicate dependency chains (19897098–19897103 and 19897104–19897111)
  were submitted simultaneously by a background waiter and a monitor both watching M29. Both
  were cancelled; neither produced scientific output.
* **`.gitignore`**: phase manifests un-ignored (`results/phase{2,3,4}/phase*_manifest.json`) as
  small JSON documentation, using the un-ignore/re-ignore/allow pattern git requires for files
  inside an ignored directory. All RDS, RData, MTX, PDF and PNG outputs remain ignored.
* **Process note — a risky command that happened to be safe.** While stopping the duplicate
  dependency chains, `for j in $(squeue -u $USER -h -o "%A"); do scancel $j; done` was used.
  That cancels **every** job belonging to the user, not only Phase 4's. It was verified
  afterwards that the queue contained only `p4_*` jobs at that instant and that the unrelated
  long-running job `3763550 cv_pert_tf_10000` was **not** affected (still RUNNING). The command
  was nonetheless wrong for the task — cancellation should be by explicit JobID — and it was not
  repeated.
* **M31 — malignant-only tumour states** (2026-09-03; SLURM 19897617 superseded, **19897639**
  reported, 00:02:22, MaxRSS 5.41 GiB): 6,434 `malignancy_refined == "Malignant"` cells subset
  from `malignancy_refined`, **not** the broad Phase 2 label. **New** reductions
  `malignant_pca` → `malignant_harmony` (batch `sample_id`) → `malignant_umap`, dims 1:30
  capturing 93.0% of variance; Phase 2's `pca`/`postint_harmony`/`postint_umap_harmony` never
  modified and M35 asserts every Phase 2 embedding is numerically identical in the final object.
  The new reduction was **justified by measurement**: kNN(k=20) overlap with the Phase 2 Harmony
  embedding is only **0.1733**. Resolution swept 0.2–1.0 with the selection rule printed before
  application; **8 states**, no cluster needing the `Uncertain` label.
* **M31 — headline negative result (§36)**: **0 of 8 malignant states are recurrent** across
  patients, 1 is shared between 2, and **97.2% of malignant cells sit in patient-private
  states**. An earlier version of this analysis classified recurrence on patient *presence* and
  therefore labelled `Mesenchymal_ECM-2` "recurrent across ≥3 patients" when 1,587 of its 1,598
  cells were MPNST_4 and the others contributed 1, 3 and 7 cells — a direct §36 violation. The
  rule was rewritten to require that no patient exceeds 80% **and** that ≥3 patients each
  contribute ≥5% and ≥10 cells, the script now prints the negative result explicitly, and M31/M33
  were re-executed. Clone counts per state range 2–19, so a state routinely spans many CNV
  clones: clone and state are **not** equivalent (§33).
* **M31 — the coherent biological finding.** Phase 2's `MPNST-Tumor` cells map to `Cycling`
  (248), `Schwann_like` (175) and `Interferon` (72), while the promoted fibroblasts map to the
  five `Mesenchymal_ECM` states (1,245 + 1,128 + 637 + 512 + 481). **Phase 2 detected the
  marker-legible malignant cells and missed the Mes-NC-like ECM ones, because an ECM programme
  is what a fibroblast looks like** — which is precisely why copy-number evidence rather than a
  better marker panel was the right instrument.
* **M32 — targeted CCC sensitivity** (2026-09-03; SLURM 19897175, 19897176, 19897177_[1-4],
  19897178, **19897179 F**, 19897616): Phase 3's `run_liana.R`, `run_cellchat.R`,
  `run_cellphonedb.py` and `build_concordance.R` reused **completely unmodified** at Phase 3
  versions, resources, thresholds, seed and expression basis — possible because Phase 3's own
  scripts read a generic `ccc_label`/`ccc_sample` contract, so only a new input-preparation
  script was needed. Refined input **15,036 cells, 14 populations** (md5 `99fd641f…`). Results:
  LIANA 31,661 supported of 164,271 rows (7,010 tumour-involving); CellPhoneDB 21,844 supported
  (6,292 tumour-involving); CellChat 18–38 min per patient, parallelised from the outset per the
  Phase 3 lesson.
* **M32 — the sensitivity answer.** supported **36,486 → 34,893 (−4.4%)**, **High concordance
  5,848 → 5,698 (−2.6%)**, Moderate 499 → 513, Discordant 5,485 → 5,380. **4,036 cells entered
  the tumour compartment and 2,015 left, and high-concordance interactions fell by only 2.6% —
  the Phase 3 architecture survives.** Change classes over 506,605 union keys: Stable 18,276 ·
  Weakened 8,446 · Newly-supported 5,421 · Lost 3,963 · **Sender-reassigned 3,051** ·
  Strengthened 2,746. **Median change across the 23 named axes −3.1%, and no axis was lost.**
  APP→CD74 −3.1%, ANXA1→FPR1 0.0%, HLA-E→KLRC1 −2.0%, HLA-F→LILRB1/2 −1.4%/−3.1%, JAG1→NOTCH3
  0.0%, JAG2→NOTCH2 0.0% all effectively unchanged; FN1→ITGAV_ITGB8 **+33.3%**, CD99→PILRA
  +14.0%, VEGFA→KDR +11.1% strengthened; COL1A1→ITGAV_ITGB8 **−44.4%**, SLIT2→ROBO1 −38.9%,
  DLL4→NOTCH2 −31.6%, COL6A2→ITGAV_ITGB8 −30.0%, VEGFA→FLT1 −28.0% genuinely weakened.
* **M32 — the fibroblast question answered (§43).** Of 9,483 Phase 3 supported interactions
  involving `Fibroblast`: **6,547 (69.0%) retained, 1,565 (16.5%) reassigned to the refined
  tumour compartment, 1,371 (14.5%) lost or no longer testable.** **All four collagen/FN1 →
  ITGAV_ITGB8 axes changed sender** — the ECM→integrin signalling Phase 3 read as
  stroma-to-tumour is substantially **tumour-autocrine**, from malignant Mes-NC-like cells to
  malignant cells. FN1 *rose* because the refined tumour compartment itself expresses it.
* **M32 failure, root-caused**: 19897179's *concordance* step exited 0 with intact outputs; only
  the sensitivity script failed, from two defects — `fib3` derived from `t3` had no `rk` column
  so `NULL %in% reasg` returned `logical(0)`, and `canon()` stripping punctuation collapsed an
  upstream resource's `"CD8 RECEPTOR"` and `"CD8_RECEPTOR"` into one key, leaving 4 duplicate
  keys per table and a many-to-many join. Fixed by computing `rk` explicitly and by **collapsing
  duplicate keys while keeping the strongest support and logging how many were affected**, rather
  than silencing the warning. **Only the failed script was rerun — the 96-minute CCC computation
  was not repeated, and no resource was increased.** Refinement also emptied populations in some
  patients (`Fibroblast` 1 cell in MPNST_4, `Pericyte-VSMC` 4), so evaluable populations fell to
  10 in MPNST_3 and 12 in MPNST_4; the change classes therefore separate `Ambiguous` (no longer
  testable) from `Lost` (testable, no longer supported).
* **M33 — tumour-state → TME model: NOT ESTABLISHABLE, and the §68 structure was not imposed.**
  Of 176 evidence rows, **12 reach ≥3 patients, 34 reach ≥2, 121 rest on one patient** — and
  **all 12 of the ≥3-patient rows belong to `Cycling`**, a 290-cell state that is itself 83% one
  patient. `Cycling` leads all six programmes purely because it is the only state present in
  enough patients to be evaluated: an **evaluability artefact, not a biological preference**. The
  five `Mesenchymal_ECM` states hold **91.3%** of the malignant compartment and have **exactly
  one evaluable patient each**. Deliberately targeted (§45): LR frameworks were **not** run per
  state; evidence is per-state, per-patient ligand (40/40 present) and receptor (39 present)
  expression at the Phase 3 detection floor, joint only when both sides pass **in the same
  patient**. **No composite score invented.** Phase 3 NicheNet and LochNESS reused verbatim, and
  §46's distinction preserved: NicheNet is not forced to support APP–CD74, and the myeloid ligand
  programme and the CSF1/IL15 receiver-state programme are carried as separate streams.
* **M34 — robustness** (SLURM 19897639): states multi-patient? **no** (0/8). CNV patterns
  patient-specific? **predominantly** — pairwise Jaccard 0.20–0.23 among reliable patients, only
  1 of 59 broad events in all four, shared core **chr18 loss / chr2 gain / chr7 gain each 3/3**.
  CCC changes driven by one patient? **the `Newly-supported` class is — 91.8% single-patient**,
  only 27 of 5,421 at ≥3 patients, whereas `Strengthened` is 10.3% single-patient with 1,031 at
  ≥3. Core interactions LOSO-stable? **yes**, retention 0.727 / 0.825 / 0.916 / 0.955.
* **M34 — the most important caveat in Phase 4.** Dropping MPNST_4 gives a refined malignant
  fraction of **0.214 against a Phase 2 fraction of 0.215** — fold change 1.00. **The headline
  17.35% → 32.63% is NOT patient-robust**: outside MPNST_4 the ~1,149 promoted fibroblasts and
  ~1,933 demoted `MPNST-Tumor` cells roughly cancel. The *fibroblast* conclusion is a separate
  claim and **is** robust — `Fibroblast` is the only disputed population majority-malignant in a
  majority of patients (0.62 / 0.55 / 0.98; 0.08 in the sanity-failed sample), while
  `Candidate-Malignant-Unresolved` rests on MPNST_1 alone and `Pericyte-VSMC` / `Uncertain` on
  MPNST_4 alone.
* **M35 — freeze** (SLURM **19897640**, COMPLETED 00:12:16, MaxRSS 13.38 GiB, three stages all
  exit 0): `results/phase4/phase4_final_object.rds` (5.64 GB, md5 `e85ba8486e456917e2483f2773bdbaf3`,
  sha256 `a658447744e5ee8621fc5a3127492ceb11bce37bac4d6f3f0f7255bdcb647ae4`), **23/23 validation
  checks passed** after saving and reloading. Preservation asserted **column by column** rather
  than by a whole-frame digest — the lesson from Phase 2 M13: all 149 pre-existing metadata
  columns byte-identical, **every Phase 2 embedding numerically identical (PCA and Harmony
  included)**, `annotation_ccc_phase3` a verbatim copy of the frozen Phase 3 layer, no
  non-Malignant cell carrying the refined `MPNST-Tumor` label, tumour states only on Malignant
  cells, 0 non-finite values across 1,656,144 embedding values, save/load md5 stable. 32 Phase 4
  metadata fields added and **nothing overwritten**. `phase4_manifest.json` (35 sections,
  `figures_missing: []`, `tables_missing: []`, 11 limitations, every JobID including failures) ·
  `PHASE4_HANDOFF.md` (27 sections) · **86 files in `figures/final/`** (16/16 required + 27
  supporting, PDF + PNG) · **47 tables in `tables/final/`** · `FIGURE_INDEX.tsv` 578 → **664
  rows**. `04_scevan_cna_heatmap` was built as a genuine ComplexHeatmap figure (genome-ordered
  chr1–22, four annotation tracks, seeded 2,000-cell subsample per patient, in-figure banner on
  sanity-failed samples) rather than substituting another figure for it.
* **Phase 4 accounting**: peak memory **60.87 GiB — 13.5% of the 450 G envelope**. Five job
  failures, every one root-caused and fixed at source (three SCEVAN 1.0.3 defects, one conda
  `set -u` interaction, two bugs of mine in the sensitivity script). **No failure was ever
  addressed by increasing RAM or walltime, and no Phase 2/3 package was downgraded to make a
  tool work.** Prohibitions respected: no spatial analysis, trajectory, RNA velocity, survival,
  treatment-response modelling, deep learning, Transformers, pan-cancer analysis, large
  condition-level DE, **inferCNV or CopyKAT**, and no second CNV method added for benchmarking.
* **STOP**: Phase 4 complete and frozen. **Phase 5 not initiated; requires separate
  authorization.**

[M35A]

Phase 4 — SCEVAN Figure Consolidation and Evidence Visualization

*Visualization milestone appended after the Phase 4 freeze. **No analysis was re-run.** SCEVAN,
`malignancy_refined`, `malignancy_confidence`, the 20 decision rules, amendments A1/A2, the SCEVAN
calls and clones, the tumour states, the CCC results and every threshold were re-used exactly as
frozen. `results/phase4/phase4_final_object.rds` was neither loaded nor modified — md5
`e85ba8486e456917e2483f2773bdbaf3` unchanged.*

* **Problem addressed**: the SCEVAN CNA evidence existed but was not visible. 267 native SCEVAN
  files sat behind `SCEVAN_NATIVE_FIGURE_INDEX.tsv`, and none of the 16 required Phase 4 final
  figures showed a native CNA heatmap. Audit: `reports/phase4/M35A_SCEVAN_FIGURE_AUDIT.md`.
* **Native figures surfaced, never redrawn** (SLURM **19899301 / 19899307 / 19899308 / 19899311 /
  19899312**, all COMPLETED ~00:04:25, ReqMem 24 G, peak MaxRSS **2.97 GiB**): the
  subclone-annotated CNA heatmap, the all-cell heatmap, the tumour-only heatmap and the consensus
  clonal CN profile for MPNST_1, MPNST_2 and MPNST_4, embedded **verbatim as rasters** so SCEVAN's
  own CNA matrix, cell ordering and subclone track are untouched — `17_scevan_native_cna_MPNST1` ·
  `18_..._MPNST2` · `19_..._MPNST4` (2 pages each) · `20_scevan_cna_reliable_patients`. The Phase 2
  identity, clone size and refined malignancy the native plot cannot carry are supplied as an
  **aligned companion panel built from the same frozen clone assignments**. Redrawing the CNA
  values would have meant reconstructing SCEVAN's internal ordering — fragile, and no scientific
  gain. **Not surfaced, with reasons stated**: the cytoband onco-heatmaps (labels unreadable at
  page scale and they invite gene-level over-reading — broad events are shown at segment resolution
  in `27` instead) and the CNA-space / expression-space UMAPs (redundant with `01/02`, `05`,
  `31_03`). **All 8 `CloneTree.png` files are blank** (pixel sd exactly 0) — the documented
  ggtree/ggplot2 4.x defect of handoff §4 and limitation K. **ggplot2 was not downgraded and clone
  phylogeny plotting was not forced.**
* **Custom figures**: `21_scevan_clone_composition_phase2` (the central evidence figure — clone ×
  Phase 2 identity, per patient, with clone size, refined-malignant fraction, fibroblast-dominance
  and disputed-triple marks) · `22_fibroblast_cna_burden_by_patient` ·
  `23_fibroblast_malignant_vs_nonmalignant_cna_profile` ·
  `24_phase2_to_phase4_malignancy_transition` · `25_MPNST3_scevan_failure_qc` ·
  `26_fibroblast_malignancy_threshold_robustness` · `27_broad_cna_recurrence` ·
  `28_phase4_scevan_evidence_summary`. PDF primary, PNG alongside for 21, 22, 24, 25, 26, 28.
* **Why clone composition is the strongest evidence**: SCEVAN builds subclones from inferred CNA
  profiles **without reference to any Phase 2 label**, so the Phase 2 composition of a CNA-defined
  clone is not circular. **MPNST_2: 4 of 4 clones fibroblast-dominated. MPNST_4: 7 of 8** (the
  eighth is 146/152 endothelial and was sent to `Ambiguous/Low` by rule R11, not promoted).
  **MPNST_1: all 7 clones mix `Candidate-Malignant-Unresolved` with `MPNST-Tumor`, and 6 of 7 also
  carry `Fibroblast` cells.**
* **Malignant vs non-malignant fibroblast-labelled cells, quantified from stored output**: Cliff's
  delta across `cnv_burden` / `cnv_frac_gain` / `cnv_frac_loss` / `cnv_mean_abs` is **0.834–0.917
  (MPNST_1, n = 512 vs 303)** and **0.781–0.890 (MPNST_2, n = 637 vs 503)**; median `cnv_burden`
  0.369 vs 0.184 and 0.311 vs 0.186. Genome-wide mean CNA profiles built from the **stored** native
  `_CNAmtx.RData` matrices give r(Fib→Mal, MPNST-Tumor→Mal) = **0.934** against r(Fib→Mal,
  Fib→Non-mal) = **0.184** in MPNST_1, and **0.972** against an immune baseline of −0.381 in
  MPNST_4. **Reported with its limit: the contrast exists in 2 of 4 patients only** — MPNST_4
  retains 1 non-malignant fibroblast and MPNST_3 contributes no malignant fibroblasts — and every
  figure prints NOT EVALUABLE rather than pooling those away. **Effect sizes and patient
  consistency are reported instead of pooled cell-level p-values**, because cells within a patient
  are not independent and a pooled test on 5,064 cells would manufacture significance from n = 4.
* **MPNST_3 shown failing rather than silently omitted**: primary↔sensitivity agreement **0.0864**
  against 0.9663–0.9974; its three "clones" are T/NK (674 cells), plasma/pDC/B (402) and
  plasma-dominated (338); six of nine canonical immune populations called ~100% malignant against
  the a priori 25% gate. Depth (median 1,594 genes/cell against 3,139 / 2,490 / 2,114) is shown as
  **context only** and the figure states that low depth alone is **not** claimed as the cause.
* **Threshold robustness separated from the fragile number**: **Fibroblast → Malignant = 4,036 at
  every tested `pop_frac_low`** (0.15, 0.20, 0.25, 0.30, 0.40), as is Candidate-Malignant-Unresolved
  at 836, while retained `MPNST-Tumor` moves 3,266 → 1,405 and the refined fraction moves 42.07% →
  32.63% — which is also not patient-robust.
* **Numerical validation**: 18 frozen-value assertions run at the head of **every** figure script
  and all passed — 19,716 cells · Fibroblast 5,064 → 4,036 / 908 / 120 ·
  Candidate-Malignant-Unresolved 836 / 395 / 0 · MPNST-Tumor 1,405 / 2,015 · refined 6,434 / 9,078 /
  3,766 / 438 · clones 7 / 4 / 3 / 8. A script that disagreed would have stopped the figure freeze,
  not plotted anyway.
* **Discrepancy found and recorded, not smoothed**: handoff §12 states "Every MPNST_1 clone mixes
  `Candidate-Malignant-Unresolved`, `Fibroblast` and `MPNST-Tumor` together". The frozen clone
  assignments show **6 of 7** — `MPNST_1_clone6` (176 cells) is 108 `Candidate-Malignant-Unresolved`
  + 68 `MPNST-Tumor` with **no `Fibroblast`**, and `MPNST_1_clone7` has 1. No call, count, clone
  assignment or conclusion changes (the two exceptions are 230 of MPNST_1's 1,952 malignant cells).
  §12 is left as written; the correction of record is handoff §28.7 and the figure plots the
  composition as it actually is.
* **A defect caught by not trusting an exit code** (SLURM **19899313**, COMPLETED, rc 0, and
  **wrong**): `jsonlite::toJSON` defaults to `digits = 4`, so re-serialising
  `phase4_manifest.json` silently rounded frozen values already in it — run agreement
  `0.99514117 → 0.9951`, refined fraction `0.32633394 → 0.3263`, every `elapsed_min`, every
  `pop_frac`, every CNA reference threshold. Caught by diffing against a pre-write backup. The
  manifest was restored and regenerated under **19899315** with `digits = NA`, and
  `m35a_finalize.R` now re-reads the file after writing and **asserts all 33 pre-existing sections
  are identical to a pre-write snapshot** — only `figures`, `tables` (appended to) and the new
  `m35a_figure_consolidation` block differ.
* **Bookkeeping**: `FIGURE_INDEX.tsv` 664 → **682 rows** (path, phase, milestone, analysis, input,
  input checksum, script, parameters, git commit, SLURM JobID, notes for every new file) ·
  `phase4_manifest.json` 35 → **36 sections** with per-figure md5 and byte size ·
  `results/phase4/tables/final/` 47 → **50 tables**.
* **Environment unchanged**: `R_env` untouched, nothing installed, upgraded or downgraded. No CNV
  method was added; inferCNV and CopyKAT were not installed. No spatial, trajectory, velocity,
  survival, treatment-response, deep-learning or pan-cancer analysis was performed.
* **STOP**: M35A complete. **Phase 4 remains complete and frozen. Phase 5 is NOT initiated and
  requires separate authorization.**
