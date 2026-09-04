# MPNST Phase 1 Project Specification

## Independent Seurat Analysis, Reproducible Workflow Engineering, Marker Discovery, and Pre-Integration Handoff

---

# 0. AUTHORITATIVE PROJECT CONFIGURATION

## 0.1 Project Root

`/local/projects-t3/lilab/vmenon/Pilot_tumor`

All project code, configuration, workflow files, tests, reports, logs, benchmarks, provenance records, and lightweight results must be organized under this project root unless explicitly configured otherwise.

Do not assume the current working directory.

Verify it before execution.

---

## 0.2 Immutable Real Input Dataset

`/local/projects-t3/lilab/vmenon/Pilot_tumor/processed_mpnst.rds`

This is a very large previously processed Seurat RDS object.

The original RDS file is immutable.

Never:

- overwrite it
- modify it in place
- move it
- rename it
- delete it
- save intermediate results over it

All derived objects must use different deterministic paths.

---

## 0.3 Known Metadata Columns

The following metadata columns are known to exist:

- `orig.ident`
- `nCount_RNA`
- `nFeature_RNA`
- `sample_id`
- `percent.mt`
- `orig.anno`
- `unintegrated_clusters`
- `seurat_clusters`
- `pca100_harmony_clusters`
- `pca100_cca_clusters`
- `pca100_mnn_clusters`
- `pca100_rpca_clusters`
- `nCount_SCT`
- `nFeature_SCT`
- `pca100.sct_harmony_clusters`
- `pca100.sct_cca_clusters`

`sample_id` is the leading candidate for defining constituent datasets/samples.

`orig.ident` must also be evaluated.

Do not assume either field is correct without auditing cardinality, cross-tabulation, patient/sample structure, cell counts, and biological meaning.

Existing cluster fields and integration-derived metadata are legacy analysis metadata.

The following must not be used to choose QC thresholds, PCs, clustering resolutions, markers, or scientific conclusions in the new Phase 1 workflow:

- `unintegrated_clusters`
- `seurat_clusters`
- `pca100_harmony_clusters`
- `pca100_cca_clusters`
- `pca100_mnn_clusters`
- `pca100_rpca_clusters`
- `pca100.sct_harmony_clusters`
- `pca100.sct_cca_clusters`
- prior integrated reductions
- prior integration-derived graphs
- prior integration-derived neighbor structures

Legacy results may be inventoried for provenance but must not influence the new Phase 1 analysis.

---

## 0.4 Conda Configuration

Conda initialization:

`source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh`

Required Conda environment:

`R_env`

Activate using:

`conda activate R_env`

All project R, Python, Snakemake, testing, workflow, and SLURM commands must use `R_env` unless GitHub Actions uses a separate CI-specific environment.

Snakemake must be installed into `R_env`.

Before installing Snakemake:

1. Activate `R_env`.
2. Record `conda info`.
3. Record `conda list`.
4. Export the environment without builds.
5. Export the explicit package specification.
6. Record R version.
7. Record `sessionInfo()`.
8. Record Python version.
9. Record critical R package versions.
10. Run a dry-run dependency solve for Snakemake.
11. Save the solver output.
12. Inspect whether the proposed transaction upgrades, downgrades, removes, or replaces critical R or single-cell packages.

Critical packages include, when installed:

- R
- Seurat
- SeuratObject
- Matrix
- sctransform
- SingleCellExperiment
- SummarizedExperiment
- BiocGenerics
- future
- future.apply
- ggplot2
- patchwork
- data.table
- harmony
- DoubletFinder
- scDblFinder

If the dry-run solve proposes material changes to critical packages:

STOP.

Report the proposed changes.

Do not install Snakemake until the researcher approves a safe installation strategy.

If the solve is safe, install Snakemake into `R_env`.

After installation:

1. Record the exact installation command.
2. Record the final package transaction.
3. Verify `which snakemake`.
4. Verify `snakemake --version`.
5. Verify `snakemake --help`.
6. Verify R starts successfully.
7. Verify Seurat loads successfully.
8. Verify SeuratObject loads successfully.
9. Verify other installed critical packages still load.
10. Record final `conda list`.
11. Export final environment YAML.
12. Export final explicit package specification.
13. Record final `sessionInfo()`.
14. Run project smoke tests.

Never silently install packages.

Every environment modification must be recorded and reflected in reproducible environment specifications.

---

## 0.5 SLURM Configuration

SLURM account:

`ihc`

SLURM partition:

`ihc`

Preferred node:

`ihc-grid-1-1-1`

Maximum CPU envelope per job:

`32 CPUs`

Maximum memory envelope per job:

`450G`

Maximum walltime envelope per job:

`72:00:00`

GPU requirement:

`None`

Project working directory:

`/local/projects-t3/lilab/vmenon/Pilot_tumor`

Standard SLURM initialization:

```bash
#!/bin/bash
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1

source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

cd /local/projects-t3/lilab/vmenon/Pilot_tumor
```

The maximum resource envelope is not the default request for every job.

Use tiered rule-level resources.

Suggested initial resource tiers:

### Tier 1: Lightweight Jobs

Examples:

- configuration validation
- report generation
- small summaries
- provenance generation
- lightweight tests

Initial request:

- CPUs: 1-2
- Memory: 4G-16G
- Walltime: 00:30:00-02:00:00

### Tier 2: Moderate Jobs

Examples:

- synthetic workflow
- dataset metadata summaries
- figure generation from small intermediate tables
- moderate QC calculations

Initial request:

- CPUs: 2-8
- Memory: 16G-64G
- Walltime: 02:00:00-12:00:00

### Tier 3: Large Dataset Jobs

Examples:

- loading large dataset-specific Seurat objects
- normalization
- SCTransform
- PCA
- neighbor construction
- UMAP
- clustering
- marker discovery

Initial request:

- CPUs: 8-24
- Memory: 64G-300G
- Walltime: 12:00:00-48:00:00

### Tier 4: Exceptional Jobs

Examples:

- loading the complete real RDS
- initial real-data audit
- combined-object construction
- operations demonstrated by benchmarks to require exceptional memory

Maximum request:

- CPUs: 32
- Memory: 450G
- Walltime: 72:00:00

Do not request Tier 4 resources without justification.

Do not repeatedly increase resources after failures without diagnosing the failure.

Record SLURM JobID, requested resources, State, ExitCode, Elapsed, AllocCPUS, MaxRSS, and relevant efficiency information when available.

---

# 1. REAL-DATA COMPUTE SAFETY — MANDATORY

The real input RDS is a very large file.

Never load, deserialize, inspect, validate, subset, summarize, or analyze the real RDS directly on the login node.

Any command that may call `readRDS()` on the real input file must execute inside a SLURM allocation.

Production real-data jobs must use `sbatch`.

`srun` may be used only for short interactive debugging after explicit CPU and memory resources have been allocated.

If execution context is uncertain, assume the current shell is a login node.

Do not load the real RDS.

The login node may be used only for:

- reading and editing text files
- reading and editing source code
- Git operations
- inspecting Git history and Git diff
- configuration validation that does not load the real RDS
- Snakemake dry runs
- lightweight synthetic-data tests
- lightweight unit tests
- inspecting small logs and reports
- `sbatch`
- `squeue`
- `sacct`
- `scontrol`
- workflow orchestration

Before submitting any real-data job:

1. Generate the required code and configuration.
2. Validate paths.
3. Validate syntax.
4. Test functions on synthetic data when applicable.
5. Run unit tests when applicable.
6. Run Snakemake dry-run when applicable.
7. Determine justified SLURM resources.
8. Submit using `sbatch`.
9. Record the JobID.
10. Monitor the job.
11. Inspect stdout.
12. Inspect stderr.
13. Inspect exit status.
14. Inspect elapsed time.
15. Inspect MaxRSS when available.
16. Inspect generated outputs.
17. Diagnose failures from evidence.
18. Update `PROGRESS.md`.

Never test whether a real-data command is "small enough" by running it on the login node.

---

# 2. AGENT INDEPENDENCE AND SOURCE OF TRUTH

This repository may be developed and reviewed using different AI coding agents.

Primary implementation agent:

Antigravity.

Independent audit agent:

Cursor.

No workflow state, scientific decision, implementation assumption, error diagnosis, or required instruction may exist only in an AI conversation.

The repository is the authoritative source of project state.

Every agent must reconstruct project state using:

- `PROJECT.md`
- `PROGRESS.md`
- Git status
- Git diff
- Git history
- configuration files
- Snakemake workflow
- source code
- tests
- CI configuration
- logs
- reports
- benchmarks
- provenance records
- phase manifests
- independent audit reports

Agents must not assume access to previous AI conversations.

Changing agents does not authorize progression to another milestone.

---

# 3. EXECUTION CONTROL — MANDATORY

This project must be implemented incrementally.

Reading `PROJECT.md` does not authorize implementation of all tasks.

The implementation agent may execute only the milestone explicitly authorized in the current researcher prompt.

Future milestones are requirements and context only.

Before modifying files for an authorized milestone:

1. Read `PROJECT.md`.
2. Read `PROGRESS.md` if it exists.
3. Inspect Git status.
4. Inspect relevant Git history.
5. Inspect existing code, configuration, tests, reports, and logs.
6. State the proposed milestone execution plan.
7. Identify scientific decisions requiring researcher approval.
8. Identify expected files to create or modify.
9. Identify expected SLURM jobs.
10. Wait for researcher approval when the current prompt explicitly requires approval before implementation.

At the end of every milestone:

1. Run required tests.
2. Run required validation.
3. Inspect all generated logs.
4. Inspect warnings.
5. Inspect errors.
6. Inspect SLURM accounting for completed real-data jobs.
7. Generate the milestone report.
8. Update `PROGRESS.md`.
9. Record files created.
10. Record files modified.
11. Record commands executed.
12. Record SLURM jobs submitted.
13. Record resource usage.
14. Record environment modifications.
15. Record unresolved issues.
16. Record scientific decisions.
17. Record recommended next milestone.
18. Summarize commit-ready changes.
19. STOP.

Do not begin the next milestone.

Successful completion of one milestone is not authorization to continue.

Wait for explicit researcher approval.

---

# 4. PRIMARY SCIENTIFIC OBJECTIVE

Build a scientifically defensible, reproducible, modular, configuration-driven, HPC-ready Phase 1 single-cell RNA-seq analysis workflow for the MPNST Seurat dataset.

The workflow must independently analyze constituent datasets before any integration or batch correction.

Phase 1 includes:

- dataset and environment auditing
- identification of constituent datasets
- independent dataset processing
- pre-filter QC
- post-filter QC
- doublet assessment
- normalization
- variable-feature selection
- scaling where appropriate
- PCA
- PC evaluation
- neighbor construction
- UMAP
- clustering resolution sweep from 0.1 through 1.0
- marker discovery for every resolution
- cluster-quality evaluation
- PC recommendation per dataset
- clustering-resolution recommendation per dataset
- combined pre-integration object construction
- pre-integration diagnostics
- integration-readiness assessment
- Snakemake orchestration
- SLURM execution
- structured logging
- synthetic-data testing
- unit testing
- integration testing
- GitHub Actions CI
- provenance capture
- frozen Phase 1 handoff

Phase 1 ends before integration.

---

# 5. STRICT PHASE 1 PROHIBITIONS

Phase 1 must not execute:

- Harmony
- Seurat CCA integration
- Seurat RPCA integration
- FastMNN
- MNN integration
- scVI integration
- BBKNN
- ComBat batch correction
- any other batch-correction method
- any integrated latent-space method
- final cell-type annotation
- automated cell-type annotation
- condition-level differential-expression inference
- pseudobulk inference
- pathway analysis
- trajectory inference
- cell-cell communication
- predictive modeling

Phase 1 may evaluate whether integration should be attempted in Phase 2.

Phase 1 must not execute integration.

---

# 6. SCIENTIFIC GUARDRAILS

## 6.1 Statistical Independence

Cells are not independent biological replicates for condition-level inference.

Phase 1 marker discovery is exploratory and intended for cluster characterization and clustering evaluation.

Do not interpret cluster markers as patient-level or condition-level differential-expression evidence.

---

## 6.2 Data Leakage

Legacy cluster labels, integration results, prior annotations, and downstream results must not influence:

- QC thresholds
- normalization selection
- variable-feature selection
- PC selection
- clustering-resolution selection
- marker-based clustering evaluation

If predictive analysis is added in a future phase, patient/sample leakage must be explicitly prevented.

---

## 6.3 Reproducibility

The workflow must:

- use deterministic random seeds where applicable
- record package versions
- record R version
- record Python version
- record Snakemake version
- record configuration
- record configuration checksum
- record Git commit
- record Git status
- record input checksum
- record intermediate object checksums
- record commands
- record SLURM JobIDs
- record logs
- record warnings
- record errors
- preserve intermediate outputs required for downstream reproducibility
- avoid undocumented manual analysis

---

# 7. REQUIRED REPOSITORY ARCHITECTURE

Create or adapt the repository toward:

```text
Pilot_tumor/
├── PROJECT.md
├── PROGRESS.md
├── README.md
├── .gitignore
├── config/
│   ├── config.yaml
│   ├── config.test.yaml
│   └── schemas/
├── workflow/
│   ├── Snakefile
│   ├── rules/
│   ├── envs/
│   ├── profiles/
│   │   └── slurm/
│   └── scripts/
├── scripts/
│   ├── R/
│   ├── python/
│   └── shell/
├── tests/
│   ├── fixtures/
│   ├── unit/
│   └── integration/
├── data/
│   └── synthetic/
├── results/
│   ├── datasets/
│   ├── combined/
│   │   └── pre_integration/
│   └── phase1_manifest.json
├── reports/
│   ├── milestones/
│   ├── datasets/
│   ├── audits/
│   ├── FIGURE_INDEX.tsv
│   ├── PROVENANCE.md
│   ├── provenance.json
│   ├── PHASE1_HANDOFF.md
│   └── INTEGRATION_PREPARATION.md
├── logs/
│   ├── slurm/
│   ├── workflow/
│   └── milestones/
├── benchmarks/
├── docs/
└── .github/
    └── workflows/
```

Modify the architecture only when there is a clear engineering reason.

Document significant changes.

Do not commit the real RDS.

Do not commit large generated Seurat objects.

---

# 8. FIGURE AND DIAGNOSTIC REQUIREMENTS

Figures and machine-readable summaries are mandatory.

Every scientifically meaningful transformation must produce sufficient diagnostics to evaluate its effect.

For every constituent dataset, generate figures for:

- raw/pre-filter QC
- post-filter QC
- pre-doublet state when applicable
- post-doublet state when applicable
- normalization diagnostics
- variable-feature diagnostics
- PCA diagnostics
- variance explained
- cumulative variance explained
- PC loadings
- PC selection
- UMAP for evaluated PC configurations when computationally reasonable
- every clustering resolution from 0.1 through 1.0
- cluster transitions
- cluster stability
- cluster sizes
- sample composition by cluster
- patient composition by cluster when patient identifiers exist
- marker heatmaps
- marker dot plots
- top-marker summaries
- recommended PC range
- recommended clustering resolution

For the combined pre-integration object, generate:

- PCA
- UMAP
- dataset composition
- sample composition
- patient composition when available
- candidate-batch composition
- biological-condition composition when available
- neighborhood composition summaries
- quantitative pre-integration batch/confounding diagnostics when feasible

All figures must:

- have informative titles
- identify the dataset
- identify the processing stage
- identify important parameters
- use deterministic filenames
- never silently overwrite results from different parameters
- be saved as PDF
- also be saved as PNG when computationally reasonable
- have machine-readable source summaries when feasible

Generate and maintain:

`reports/FIGURE_INDEX.tsv`

Required columns:

- figure_path
- dataset
- processing_stage
- analysis_method
- parameters
- input_object_checksum
- generating_script
- snakemake_rule
- git_commit

Phase 1 is incomplete if a major transformation lacks sufficient diagnostics.

---

# 9. SNAKEMAKE-FIRST DEVELOPMENT

Snakemake must be developed from the beginning.

Do not first build a collection of manual scripts and wrap them in Snakemake at the end.

Every milestone that introduces production analysis functionality must introduce or update:

- modular analysis code
- corresponding Snakemake rules
- configuration
- explicit inputs
- explicit outputs
- logs
- benchmarks for expensive rules
- tests
- synthetic-data execution where applicable

The same production workflow logic must be exercised by synthetic tests.

Do not create a separate toy workflow that bypasses production logic.

Requirements:

- modular rules
- configuration-driven paths
- wildcard-based dataset execution
- wildcard-based resolution execution when appropriate
- explicit dependencies
- explicit outputs
- one log per job/rule
- benchmark files for expensive rules
- rerun safety
- meaningful failure propagation
- no hard-coded analysis paths inside scripts
- local lightweight test execution
- SLURM production execution

Before real-data execution:

1. validate configuration
2. run applicable unit tests
3. run applicable synthetic tests
4. run Snakemake dry-run
5. inspect planned jobs
6. submit production work through SLURM

---

# 10. LOGGING AND ERROR DIAGNOSIS

Every workflow stage must generate logs.

Record when applicable:

- timestamp
- milestone
- dataset
- processing stage
- input
- output
- parameters
- random seed
- environment
- package versions
- Git commit
- warnings
- errors
- runtime
- peak memory
- SLURM JobID

Do not silently suppress warnings.

Do not silently catch exceptions and continue with incomplete outputs.

After failure:

1. inspect stdout
2. inspect stderr
3. inspect workflow log
4. inspect SLURM ExitCode
5. inspect SLURM State
6. inspect MaxRSS
7. inspect elapsed time
8. inspect partial outputs
9. identify root cause
10. explain root cause
11. propose correction
12. modify code only after evidence-based diagnosis
13. rerun the smallest appropriate test
14. then rerun the failed production stage

Do not repeatedly patch symptoms.

Do not repeatedly increase memory or walltime without evidence.

---

# 11. SYNTHETIC DATA AND CI/CD

Create a deterministic synthetic Seurat dataset.

It must contain:

- multiple constituent datasets
- multiple samples
- multiple patients
- at least two biological conditions
- known cell populations
- known marker genes
- simulated technical structure
- known biological structure

Use fixed random seeds.

Document generation.

Synthetic data tests software correctness.

Synthetic data does not establish biological validity.

GitHub Actions CI must eventually:

1. install or create the CI environment
2. validate configuration
3. generate synthetic data
4. run unit tests
5. run Snakemake dry-run
6. execute a reduced PC/resolution configuration
7. execute the synthetic workflow through production logic
8. verify expected outputs
9. verify critical files are non-empty
10. verify marker files exist
11. verify logs exist
12. verify provenance exists
13. fail when required outputs are missing

CI must not require:

- real patient data
- HPC access
- secrets
- absolute HPC paths

---

# 12. DATASET AUDIT REQUIREMENTS

The first real-data milestone must audit the RDS read-only through SLURM.

Determine:

- object class
- object version
- Seurat compatibility
- object size on disk
- object size in memory
- dimensions
- assays
- assay classes
- layers
- slots
- default assay
- raw-count availability
- normalized-data availability
- scaled-data availability
- variable features
- feature names
- duplicate features
- cell names
- duplicate cell names
- metadata columns
- metadata types
- metadata cardinality
- missing metadata
- `sample_id` cardinality
- `orig.ident` cardinality
- cross-tabulation of `sample_id` and `orig.ident`
- candidate patient identifiers
- candidate sample identifiers
- candidate dataset identifiers
- candidate batch variables
- biological conditions
- existing annotations
- QC metrics
- existing doublet information
- dimensional reductions
- graphs
- active identities
- existing clustering results
- cells per candidate dataset
- cells per sample
- cells per patient when available
- object provenance when inferable

Explicitly determine whether raw RNA counts are available.

Explicitly determine whether independent reprocessing from counts is scientifically and computationally possible.

Generate:

`reports/DATA_AUDIT.md`

and:

`reports/DATASET_STRUCTURE.tsv`

`DATA_AUDIT.md` must separate:

- observed facts
- inferences
- unknowns
- scientific risks
- computational risks
- legacy-analysis contamination risks
- recommendations
- researcher decisions required

---

# 13. INDEPENDENT DATASET PROCESSING

After researcher approval of the audit and Cursor audit gate, identify the validated constituent-dataset field.

Every constituent dataset must be processed independently.

Do not manually duplicate scripts per dataset.

Use common modular functions and Snakemake wildcard logic.

For every dataset:

1. create a deterministic immutable working input object
2. record cell count
3. record feature count
4. record sample composition
5. record patient composition when available
6. record assay structure
7. record metadata summary
8. calculate checksum
9. preserve provenance

---

# 14. QUALITY CONTROL

## 14.1 Pre-Filter QC

Generate:

- nFeature_RNA distributions
- nCount_RNA distributions
- percent.mt distributions
- percent.ribo distributions
- nCount_RNA versus nFeature_RNA
- sample-stratified distributions
- patient-stratified distributions when available
- cells per sample
- cells per patient when available
- detected genes per cell
- library-size summaries

Save machine-readable summaries.

---

## 14.2 QC Recommendation

QC thresholds must be configuration-driven.

Do not silently choose thresholds.

Generate per dataset:

`reports/datasets/<DATASET_ID>/QC_RECOMMENDATION.md`

Include:

- observed distributions
- outliers
- proposed thresholds
- justification
- expected removed cell count
- expected removed fraction
- sample-specific effects
- patient-specific effects
- possible biological risks
- researcher decisions required

Filtering must use approved/configured thresholds.

---

## 14.3 Post-Filter QC

Generate the same core figures after filtering.

Generate direct pre/post comparisons.

Record:

- cells before
- cells after
- fraction removed
- removal reasons
- removals per sample
- removals per patient when available

---

# 15. DOUBLET ASSESSMENT

Determine whether valid doublet results already exist.

Document provenance.

If doublet detection is required, implement it as a configurable workflow stage.

Doublet detection must operate per independent capture/sample/library when scientifically appropriate.

Never blindly run doublet detection on the entire merged object.

Never silently remove doublets.

Generate per dataset:

`reports/datasets/<DATASET_ID>/DOUBLET_REPORT.md`

---

# 16. NORMALIZATION AND VARIABLE FEATURES

Support configuration-defined:

- LogNormalize
- SCTransform

Recommend an approach based on observed dataset structure.

Do not run multiple normalization strategies solely to generate more outputs.

Record all parameters.

Save deterministic intermediate objects.

Generate normalization and variable-feature diagnostics.

---

# 17. PCA AND PC EVALUATION

Run PCA independently for every dataset.

Maximum PCs must be configurable.

Evaluate candidate PC ranges when available:

- 1:10
- 1:15
- 1:20
- 1:25
- 1:30
- 1:40
- 1:50

Generate:

- elbow plot
- variance explained
- cumulative variance explained
- top positive loadings
- top negative loadings
- PC heatmaps when reasonable
- PCA by sample
- PCA by patient when available
- PCA by candidate batch
- PCA by biological condition when available

PC recommendations must consider:

- variance structure
- elbow behavior
- biological signal
- technical structure
- neighborhood stability
- clustering stability
- marker coherence
- sensitivity to PC choice

An elbow plot alone must not select PCs.

Generate:

`reports/datasets/<DATASET_ID>/PC_RECOMMENDATION.md`

---

# 18. CLUSTERING RESOLUTION SWEEP

For every dataset evaluate:

- 0.1
- 0.2
- 0.3
- 0.4
- 0.5
- 0.6
- 0.7
- 0.8
- 0.9
- 1.0

The resolution grid must be configuration-driven.

Preserve every clustering solution using explicit metadata names.

Never overwrite prior solutions.

For every resolution generate:

- UMAP labeled by cluster
- UMAP colored by sample
- UMAP colored by patient when available
- cluster-size table
- sample composition per cluster
- patient composition per cluster when available
- cluster transitions
- clustering stability metrics when feasible
- silhouette or related diagnostics when appropriate
- warnings for tiny clusters
- warnings for sample-dominated clusters
- warnings for patient-dominated clusters
- warnings for unstable splitting

Generate:

- multi-resolution comparison figure
- clustree-style visualization when feasible
- machine-readable resolution comparison table

---

# 19. MARKER DISCOVERY

Perform exploratory marker discovery for every clustering resolution.

Record:

- assay
- layer/slot
- statistical test
- min.pct
- log-fold-change threshold
- only.pos
- latent variables when used
- random seed when relevant

Generate for every resolution:

- complete marker table
- filtered marker table
- top markers per cluster
- marker heatmap
- marker dot plot
- marker diagnostic summary

Use the appropriate expression assay.

Do not use integrated embeddings.

Do not use legacy integrated results.

Do not interpret cluster markers as condition-level inferential DE.

---

# 20. DATASET-SPECIFIC RECOMMENDATIONS

Generate per dataset:

`reports/datasets/<DATASET_ID>/ANALYSIS_RECOMMENDATION.md`

Recommend:

1. QC thresholds
2. normalization strategy
3. PC range
4. alternative PC range
5. clustering resolution
6. alternative clustering resolution
7. whether technical structure remains
8. whether sample dominance affects clusters
9. whether patient dominance affects clusters
10. marker coherence
11. small-cluster stability
12. limitations
13. unresolved questions

Resolution recommendation must consider:

- stability
- marker coherence
- cluster sizes
- sample representation
- patient representation
- over-fragmentation
- under-clustering
- sensitivity to PC choice
- biological interpretability

No single metric may automatically choose the resolution.

Generate a ranked comparison table.

---

# 21. COMBINED PRE-INTEGRATION BASELINE

After all independent datasets are completed and approved, construct a scientifically appropriate combined pre-integration object.

Do not simply reuse the previously processed merged object without justification.

Preserve:

- dataset provenance
- sample provenance
- patient provenance
- batch provenance
- object checksums
- processing configuration

Generate:

- PCA
- UMAP
- dataset composition
- sample composition
- patient composition when available
- candidate-batch composition
- biological-condition composition when available
- neighborhood composition summaries
- quantitative batch/confounding diagnostics when feasible

Save under:

`results/combined/pre_integration/`

Generate:

`reports/PRE_INTEGRATION_ASSESSMENT.md`

---

# 22. INTEGRATION PREPARATION

Phase 1 must evaluate whether integration should be attempted in Phase 2.

Determine:

- candidate technical variables
- candidate integration variables
- confounding with biological condition
- confounding with patient
- confounding with dataset
- risk of removing biological signal
- whether Harmony should be evaluated in Phase 2
- whether non-integrated analysis should be preserved as primary

Generate:

`reports/INTEGRATION_PREPARATION.md`

Do not execute integration.

---

# 23. PHASE 1 HANDOFF

Generate:

`reports/PHASE1_HANDOFF.md`

and:

`results/phase1_manifest.json`

For every dataset record:

- dataset identifier
- original provenance
- filtered object path/checksum
- doublet-handled object path/checksum when applicable
- normalized object path/checksum
- normalization method
- variable-feature parameters
- PCA object path/checksum
- recommended PC range
- alternative PC range
- clustering results 0.1-1.0
- recommended resolution
- alternative resolution
- marker paths for every resolution
- QC thresholds
- cell counts before/after transformations
- unresolved warnings
- unresolved confounding
- scientific limitations

For the combined pre-integration object record:

- construction method
- path
- checksum
- constituent datasets
- cell counts
- sample counts
- patient counts when available
- normalization status
- PCA configuration
- pre-integration diagnostics
- candidate integration variables
- confounding assessment
- integration recommendation

Also record:

- Git commit
- Git status
- configuration checksum
- input checksum
- Conda environment
- package versions
- Snakemake version
- workflow completion status
- failed rules
- skipped rules
- SLURM jobs
- provenance paths

Phase 2 must consume the Phase 1 manifest.

Phase 2 must not silently overwrite Phase 1 outputs.

After handoff:

STOP.

---

# 24. MILESTONE PLAN

## M0 — Infrastructure, Environment, and SLURM Safety

No real RDS access.

Required work:

- inspect repository
- inspect Git state
- activate `R_env`
- capture environment baseline
- inspect Snakemake availability
- dry-run solve if installation required
- safely install Snakemake only if permitted by dependency-safety rules
- verify critical R packages
- create/adapt repository architecture
- create `PROGRESS.md`
- create base configuration
- create configuration validation
- create Snakemake skeleton
- create SLURM profile/templates
- create structured logging skeleton
- create provenance skeleton
- create deterministic synthetic smoke fixture
- create initial unit/smoke tests
- execute lightweight tests
- execute Snakemake dry-run
- submit a trivial SLURM smoke job using `sbatch`
- monitor it
- inspect stdout/stderr/accounting
- generate `reports/milestones/M0_REPORT.md`
- update `PROGRESS.md`
- STOP

Do not access the real RDS.

---

## M1 — Real-Data Audit

First real-data milestone.

Required work:

- read current repository state
- inspect M0 outputs
- write modular audit code
- add Snakemake audit rule
- test audit logic on synthetic data
- run tests
- run dry-run
- submit real audit through `sbatch`
- use Tier 4 resources initially because the object is known to be very large
- monitor job
- inspect logs
- inspect accounting
- diagnose failures
- generate `DATA_AUDIT.md`
- generate `DATASET_STRUCTURE.tsv`
- recommend validated constituent-dataset field
- identify researcher decisions
- generate `reports/milestones/M1_REPORT.md`
- update `PROGRESS.md`
- STOP

After M1:

Cursor independent audit gate.

Do not proceed until researcher approval.

---

## M2 — Dataset Extraction and Pre-Filter QC

Required work:

- implement deterministic dataset extraction
- use validated dataset identifier
- preserve provenance
- add Snakemake rules
- add synthetic tests
- submit real jobs through SLURM
- generate dataset-specific immutable working inputs
- generate pre-filter QC figures
- generate machine-readable QC summaries
- update FIGURE_INDEX.tsv
- generate milestone report
- update PROGRESS.md
- STOP

---

## M3 — QC Filtering and Doublet Assessment

Required work:

- generate QC recommendations
- require configured/approved thresholds
- implement filtering
- generate post-filter QC
- generate pre/post comparisons
- assess existing doublet information
- implement configurable doublet stage if required
- add Snakemake rules
- add tests
- run real jobs through SLURM
- inspect logs/accounting
- update figures/index
- generate milestone report
- update PROGRESS.md
- STOP

After M3:

Cursor independent audit gate.

Do not proceed until researcher approval.

---

## M4 — Normalization and Variable Features

Required work:

- recommend normalization strategy
- implement configured normalization
- implement variable-feature selection
- generate diagnostics
- add Snakemake rules
- add tests
- run real jobs through SLURM
- inspect logs/accounting
- update figure index
- generate milestone report
- update PROGRESS.md
- STOP

---

## M5 — PCA and PC Evaluation

Required work:

- run PCA
- generate all PCA diagnostics
- evaluate candidate PC ranges
- evaluate sensitivity where feasible
- recommend PC ranges
- add Snakemake rules
- add tests
- run real jobs through SLURM
- inspect logs/accounting
- update figure index
- generate milestone report
- update PROGRESS.md
- STOP

---

## M6 — Clustering Resolution Sweep

Required work:

- construct neighbors
- run UMAP
- run resolution sweep 0.1-1.0
- preserve every clustering solution
- generate all required figures
- calculate cluster-quality diagnostics
- generate transition/stability outputs
- add Snakemake rules
- add tests
- run real jobs through SLURM
- inspect logs/accounting
- update figure index
- generate milestone report
- update PROGRESS.md
- STOP

---

## M7 — Marker Discovery and Dataset Recommendations

Required work:

- perform marker discovery for every resolution
- generate complete and filtered marker tables
- generate heatmaps
- generate dot plots
- generate marker diagnostics
- evaluate clustering solutions
- recommend PC/resolution per dataset
- generate ranked comparisons
- add Snakemake rules
- add tests
- run real jobs through SLURM
- inspect logs/accounting
- update figure index
- generate milestone report
- update PROGRESS.md
- STOP

After M7:

Cursor independent audit gate.

Do not proceed until researcher approval.

---

## M8 — Combined Pre-Integration Baseline

Required work:

- construct justified combined pre-integration object
- preserve provenance
- generate pre-integration PCA
- generate pre-integration UMAP
- generate composition figures
- generate batch/confounding diagnostics
- generate PRE_INTEGRATION_ASSESSMENT.md
- generate INTEGRATION_PREPARATION.md
- add Snakemake rules
- add tests
- run real jobs through SLURM
- inspect logs/accounting
- update figure index
- generate milestone report
- update PROGRESS.md
- STOP

---

## M9 — Workflow Hardening, CI/CD, Provenance, and Phase 1 Freeze

Required work:

- audit complete Snakemake DAG
- harden configuration validation
- harden rerun safety
- harden failure propagation
- complete unit tests
- complete synthetic integration tests
- complete GitHub Actions CI
- verify CI uses production workflow logic
- validate environment specifications
- finalize README
- finalize logging
- finalize benchmarks
- finalize provenance
- run end-to-end synthetic workflow
- validate required outputs
- inspect real-data workflow completeness
- generate PHASE1_HANDOFF.md
- generate phase1_manifest.json
- generate M9_REPORT.md
- update PROGRESS.md
- freeze Phase 1
- STOP

After M9:

Cursor final independent audit gate.

Do not begin Phase 2.

---

# 25. INDEPENDENT CURSOR AUDIT GATES

Cursor audits must use fresh sessions when feasible.

Cursor must initially operate read-only.

Cursor must not assume access to Antigravity conversations.

Audit gates:

- after M1
- after M3
- after M7
- after M9

Cursor must inspect:

- PROJECT.md compliance
- PROGRESS.md
- Git diff/history
- code
- configuration
- Snakemake DAG
- tests
- CI
- logs
- SLURM records
- figures
- reports
- provenance

Cursor must evaluate:

- scientific validity
- data leakage
- legacy-analysis contamination
- pseudoreplication risks
- incorrect statistical units
- QC bias
- doublet handling
- normalization misuse
- PC-selection quality
- clustering stability
- marker-analysis validity
- sample/patient dominance
- hard-coded paths
- environment drift
- missing dependencies
- non-determinism
- stale outputs
- incorrect Snakemake dependencies
- weak tests
- CI bypass of production logic
- swallowed errors
- missing logs
- missing provenance
- documentation mismatch

For every finding report:

1. Severity: Critical, High, Medium, or Low.
2. Exact file and line.
3. Evidence.
4. Scientific or computational consequence.
5. Recommended correction.
6. Regression test.

Cursor must not fix findings until the researcher reviews the audit.

---

# 26. DEVELOPMENT RULES

1. Inspect before editing.
2. Execute only the authorized milestone.
3. Never load the real RDS on the login node.
4. Use `sbatch` for real-data computation.
5. Use `srun` only after explicit interactive allocation.
6. Use `R_env`.
7. Install Snakemake safely and reproducibly.
8. Record every environment modification.
9. Never modify the original RDS.
10. Never fabricate metadata.
11. Never silently assume the dataset identifier.
12. Never use legacy integrated results to influence Phase 1.
13. Process constituent datasets independently.
14. Build Snakemake from the beginning.
15. Use production workflow logic in synthetic tests.
16. Generate figures at every scientifically meaningful stage.
17. Generate machine-readable summaries.
18. Evaluate multiple PC ranges.
19. Sweep resolutions 0.1-1.0.
20. Preserve every clustering solution.
21. Generate markers for every resolution.
22. Do not perform cell-type annotation.
23. Do not perform condition-level inferential DE.
24. Do not execute integration.
25. Do not use UMAP appearance as the sole selection criterion.
26. Do not use one metric as the sole clustering-resolution criterion.
27. Run small tests before expensive jobs.
28. Inspect logs after failures.
29. Diagnose root cause before fixes.
30. Do not repeatedly patch symptoms.
31. Do not repeatedly increase resources without evidence.
32. Do not weaken tests merely to make CI pass.
33. Do not delete failing tests without justification.
34. Do not hard-code analysis paths inside scripts.
35. Keep real data out of Git.
36. Keep large generated objects out of Git.
37. Update PROGRESS.md at every milestone.
38. Maintain FIGURE_INDEX.tsv.
39. Maintain provenance.
40. Ensure a fresh agent can reconstruct the project from the repository.
41. STOP at every milestone boundary.
42. Wait for explicit researcher approval before continuing.

---

# 27. DEFINITION OF PHASE 1 COMPLETION

Phase 1 is complete only when:

- M0-M9 are completed
- all required Cursor audit gates have been reviewed
- approved audit findings are resolved
- environment is reproducibly specified
- Snakemake production workflow is functional
- SLURM execution is functional
- synthetic tests are functional
- CI is functional
- real-data audit is complete
- constituent datasets are independently processed
- pre/post QC figures exist
- doublet assessment exists
- normalization diagnostics exist
- variable-feature diagnostics exist
- PCA diagnostics exist
- PC recommendations exist
- resolution sweeps 0.1-1.0 exist
- clustering diagnostics exist
- markers exist for every resolution
- dataset recommendations exist
- combined pre-integration baseline exists
- integration-readiness assessment exists
- FIGURE_INDEX.tsv is complete
- logs are complete
- benchmarks are complete
- provenance is complete
- PHASE1_HANDOFF.md exists
- phase1_manifest.json exists
- PROGRESS.md is current
- Phase 1 outputs are frozen

After Phase 1 completion:

STOP.

Do not begin Harmony or any Phase 2 integration analysis without explicit researcher authorization.

---
---

# PHASE 2 PROJECT ARCHITECTURE (added 2026-09-03 · amendment: CCC-oriented annotation)

Phase 1 above remains the authoritative Phase 1 specification and is frozen. This section
records the **actual final Phase 2 architecture** as built, including the approved
CCC-oriented annotation amendment.

## P2.1 Pipeline architecture

```text
Phase 1
  QC / doublet removal / SCTransform / per-sample PCA / per-sample clustering / markers
  -> frozen handoff: results/combined/pre_integration/combined_preintegration.rds
        |
        v
Phase 2
  Harmony (embedding-level, group.by.vars = sample_id, dims 1:30, package defaults)   [M11]
        |
        v
  pre/post integration assessment (mixing + biological preservation)                  [M12]
        |
        v
  neighbours -> UMAP -> Louvain resolution sweep 0.1-1.0 -> primary resolution 1.0     [M13]
        |
        v
  cluster markers (SCT assay, PrepSCTFindMarkers, NOT Harmony coordinates)             [M14]
        |
        v
  literature-grounded DETAILED annotation (Level 1/2/3 + confidence)                   [M15]
        |
        v
  computational refinement + descriptive composition                                  [M16]
        |
        v
  CCC-ORIENTED ANNOTATION (annotation_ccc)                                            [M15A]
        |
        v
  final validated Seurat object + manifest + handoff                                   [M17]
```

**Snakemake is not used in Phase 2.** Each milestone is a standalone argument-driven R
script under `scripts/R/phase2/` invoked by an explicit `sbatch` script under
`scripts/shell/phase2/`.

## P2.2 Hierarchical annotation strategy

The Phase 2 final object contains **both** annotation layers. Neither replaces the other.

```text
Detailed annotation  (postint_celltype_level1/2/3, postint_annotation_confidence)
        |
        v
Evidence-supported malignant populations
        |
        v
Collapse malignant populations
        |
        v
MPNST-Tumor                                     <-- annotation_ccc

Immune / stromal / endothelial populations
        |
        v
retain biologically meaningful identities        <-- annotation_ccc
```

1. **Detailed literature-supported biological annotation** — the primary biological
   annotation and the **source of truth**. Assigned per cluster from the M14 markers plus
   cited primary literature, recorded in `config/phase2/annotation_map_M15.tsv` and
   `results/phase2/annotation/ANNOTATION_EVIDENCE.tsv`.
2. **Simplified CCC-oriented annotation** (`annotation_ccc`) — a derived layer built for
   downstream tumour–microenvironment cell–cell communication analysis, in which the
   evidence-supported malignant states are collapsed into a single `MPNST-Tumor` population
   while immune, stromal and endothelial populations retain their biologically meaningful
   identities.

### P2.3 The defining constraint on `MPNST-Tumor`

> **`MPNST-Tumor` is NOT defined as "everything that is not immune."**

Only populations with sufficient positive evidence of malignant / MPNST identity may be
collapsed into `MPNST-Tumor`. The following logic is **prohibited**:

```r
if (!immune) annotation_ccc <- "MPNST-Tumor"   # FORBIDDEN
```

Such logic would incorrectly absorb fibroblasts, endothelial cells, pericytes,
non-malignant Schwann-lineage cells, ambiguous populations and low-quality clusters into
the malignant compartment.

The collapse rule is instead **evidence-based and mechanical**, applied to the detailed
Level 2 label, which already encodes the strength of the marker and literature evidence:

| Detailed Level 2 label | Collapsed to `MPNST-Tumor`? | Basis |
| --- | :---: | --- |
| `MPNST-like malignant (SCP-like)` | **YES** | Positive MPNST lineage-marker evidence (L1CAM, MPZ) plus cited literature |
| `MPNST-like malignant (NC-like)` | **YES** | Positive evidence (SHH, GAL3ST1, KLK6) plus cited literature |
| `Schwann-lineage tumour-like` | **YES** | Coherent Schwann/neural tumour-like programme with no alternative lineage identity |
| `Cycling tumour-like` | **YES** | Proliferating tumour population within the malignant compartment |
| `Candidate malignant` | **NO** | Label itself is provisional; Low confidence; evidence is patient-private neural expression only. Retained separately as `Candidate-Malignant-Unresolved` |
| Anything not in Level 1 `Malignant / tumour` | **NO** | Never eligible |

Ambiguous and low-quality populations are **never** forced into `MPNST-Tumor`; they retain
`Uncertain` or an explicit exclusion label. Every decision, including every *exclusion*, is
recorded with its rationale in `results/phase2/annotation/CCC_ANNOTATION_MAPPING.tsv`.

### P2.4 Tumour heterogeneity is preserved, not destroyed

Collapsing malignant states for communication analysis must not permanently hide tumour
heterogeneity. The final object therefore retains the detailed malignant state labels
alongside `annotation_ccc`, so both of the following are possible **without re-running
annotation**:

```text
Broad CCC analysis                 Later: tumour-state-specific CCC
MPNST-Tumor                        MPNST-like malignant (SCP-like)
     <->                           MPNST-like malignant (NC-like)
Macrophage / CD8 T / NK /          Schwann-lineage tumour-like
Fibroblast / Endothelial ...       Cycling tumour-like
                                        <-> immune / TME populations
```

### P2.5 Statistical rule carried into Phase 3

Phase 3 must **not** pool every cell across every patient and treat cells as biological
replicates. `sample_id` is simultaneously the dataset and the patient, and it is preserved
verbatim in the final object precisely so that communication analysis can be computed
**per sample** and only then compared across samples:

```text
Per sample:  MPNST-Tumor -> Macrophage
Per sample:  MPNST-Tumor -> CD8 T
Per sample:  Macrophage  -> MPNST-Tumor
        then: patient/sample-aware comparison
```

With four tumours, one sample per patient and **no biological-condition variable anywhere
in the data**, there is no replication structure for condition-level inference. Differential
abundance and pseudobulk differential expression are therefore not supportable on this
dataset and remain out of scope.

## P2.6 Annotation metadata fields in the final object

| Field | Content |
| --- | --- |
| `postint_celltype_level1` / `_level2` / `_level3` | Detailed annotation (M15) — source of truth |
| `postint_annotation_confidence` | High / Moderate / Low / Uncertain (M15) |
| `postint_celltype_level*_initial`, `postint_annotation_confidence_initial` | M15 values, preserved |
| `postint_celltype_level*_refined`, `postint_annotation_confidence_refined` | M16 review outputs |
| `postint_annotation_initial` / `postint_annotation_refined` | Level 2 before / after review |
| `postint_annotation_review_rule` / `_outcome` | Why the M16 review confirmed or downgraded |
| **`annotation_ccc`** | **CCC-oriented identity (this amendment)** |
| `annotation_ccc_compartment` | Broad CCC compartment (MPNST-Tumor / Immune / Fibroblast-Stromal / Endothelial / Other / Uncertain) |
| `annotation_ccc_is_tumor` | TRUE only for cells collapsed into `MPNST-Tumor` |
| `annotation_ccc_ccc_ready` | FALSE for populations flagged as unsuitable for communication inference |

## P2.7 Scientific prohibitions still in force

No condition-level pseudobulk DE · no cell-level condition DE treated as replication · no
differential abundance inference · no survival analysis · no trajectory inference · no RNA
velocity · **no CNV inference** · no CellChat / CellPhoneDB / LIANA / NicheNet / LochNESS ·
no perturbation, predictive, deep-learning or Transformer modelling · no Snakemake in
Phase 2. These are Phase 3 candidates and require separate authorization.

---
---

# PHASE 3 PROJECT ARCHITECTURE (added 2026-09-03)

Phases 1 and 2 above are frozen. This section records the Phase 3 architecture.

## P3.1 Biological objective

Determine how MPNST tumour cells interact with immune and other tumour-microenvironment
populations, identify recurrent and biologically plausible signalling circuits, establish
which findings are concordant across independent CCC frameworks, and connect candidate
tumour-derived signals to receiver-cell transcriptional states.

**Phase 3 is not a methods-benchmarking study.** Multiple tools are used because CCC
inference is method-dependent; an interaction supported by several independent frameworks
earns higher confidence. Method disagreement is preserved and reported.

## P3.2 Pipeline

```text
Phase 1  QC / preprocessing
Phase 2  Harmony -> clustering -> markers -> annotation -> CCC-oriented annotation
Phase 3
  sample-aware CCC input preparation                                     [M19]
        |
        v
  multi-method inference: LIANA | CellChat | CellPhoneDB                 [M20]
        |
        v
  concordance (canonical keys, per-method support flags, classes)        [M21]
        |
        v
  MPNST tumour <-> immune / stromal / endothelial prioritisation         [M22]
        |
        v
  receiver-response modelling (NicheNet: ligand -> receiver programme)   [M23]
        |
        v
  LochNESS receiver-state analysis (orthogonal neighbourhood phenotype)  [M24]
        |
        v
  integrated evidence matrix                                             [M25]
        |
        v
  robustness / sample-patient consistency                                [M26]
        |
        v
  integrated tumour-immune interaction map + freeze                      [M27]
```

## P3.3 Role of each method

| Method | Role | Why it is not interchangeable with the others |
| --- | --- | --- |
| **LIANA** | Primary multi-resource / multi-score **consensus LR** framework | Aggregates several scoring schemes over a consensus resource; output is a rank |
| **CellChat** | Independent **pathway- and network-oriented** framework | Adds pathway-level communication and outgoing/incoming signalling roles that LR lists do not express |
| **CellPhoneDB** | Independent **permutation-based LR** framework | Gives an empirical p-value against a label-shuffled null |
| **NicheNet** | **Receiver-response modelling**, not a fourth LR list | Tests whether a sender ligand explains a receiver's *transcriptional programme*, a different claim from receptor co-expression |
| **LochNESS** | **Orthogonal transcriptional-neighbourhood phenotype** | Detects local enrichment in receiver-state space. **It is not a CCC method** — see P3.5 |

**Raw scores are never averaged across tools.** A CellChat probability, a CellPhoneDB
p-value, a LIANA rank and a NicheNet ligand activity are different quantities on different
scales. Each contributes a boolean support flag; those flags form a transparent concordance
class. NicheNet is represented as orthogonal receiver-response support, not as a fourth
equivalent LR framework.

**Known partial non-independence:** LIANA's method set includes a CellPhoneDB-style score,
so LIANA and the standalone CellPhoneDB run are not fully independent. This is stated in the
concordance report rather than glossed over.

## P3.4 Sample-aware analysis is mandatory

`sample_id` is simultaneously the dataset and the patient, and there is **no
biological-condition variable**. CCC is therefore computed **per sample** and only then
assessed for recurrence:

```text
Per sample:  MPNST-Tumor -> Macrophage
Per sample:  MPNST-Tumor -> CD8-T
Per sample:  Macrophage  -> MPNST-Tumor
        then: patient/sample-aware recurrence, not cell-level pooling
```

Minimum-cell policy (stated, not silently chosen): populations flagged not CCC-ready in
Phase 2 are excluded outright; a population is EVALUABLE in a sample only with **≥ 10 cells**
(CellChat's `min.cells` default); a sample is ANALYSABLE with ≥ 3 evaluable populations; and
**a population below the minimum is recorded NOT EVALUABLE, never as absence of signalling.**

## P3.5 LochNESS — role and the rule that constrains it

LochNESS (Huang et al., *Nature* 623, 772–781, 2023; doi 10.1038/s41586-023-06548-w) detects
local enrichment or depletion of a population in transcriptional-neighbourhood space.

> **LochNESS is not a ligand–receptor inference method. Phase 3 never claims that high
> LochNESS means cell–cell communication.**

Its role is strictly downstream and orthogonal:

```text
CCC -> candidate sender->receiver signalling -> receiver-response analysis
    -> LochNESS receiver-state structure -> additional evidence about state heterogeneity
```

Implementation follows the official MMCA R formulation — `k = round(0.5·√N)`, L2-normalised
PCA, exact `FNN::get.knnx`, **disjoint query/reference sets so same-sample neighbours are
excluded**, global fraction computed over the reference set, and a label-permutation null
with the exclusion structure preserved. The perturb-seq Python implementation is retained as
an independent cross-check. `external/MMCA_ref/` holds a pinned read-only copy of the
official scripts (commit `af629c49`).

## P3.6 Expression basis

CCC uses the **`RNA` assay, joined, LogNormalize**, not SCT: the SCT assay carries four
SCTransform models, so SCT values are not on a common footing across samples.
**Harmony coordinates are never used as expression** — they inform neighbourhood structure
for LochNESS only.

## P3.7 Statistics

No cell-level pseudoreplication. Emphasis on sample recurrence, patient recurrence, effect
size, concordance and receiver response rather than p-values over 17,327 cells. Permutation
nulls are sample-aware.

## P3.8 Interpretation language

Conclusions use *predicted / inferred / supported / consistent with / concordant /
associated with*. Every Phase 3 report states:

> scRNA-seq ligand–receptor analysis infers communication **potential** from expression and
> associated transcriptional programmes. It does not establish physical cell adjacency or
> direct signalling.

## P3.9 Prohibitions

No survival analysis · no treatment-response modelling · no deep learning or Transformers ·
no trajectory inference · no RNA velocity · **no CNV inference** · **no spatial inference** ·
no pseudobulk condition DE unrelated to receiver-response questions · no large unrelated
pathway screens · no arbitrary composite scores without stated justification.

## P3.10 Environment

Phase 3 remains centred on `R_env` (R 4.4.3, Seurat 5.4.0, SeuratObject 5.3.0, harmony 1.2.4,
Matrix 1.7.4 — all unchanged from Phase 2). CellPhoneDB runs in a deliberately isolated
`cpdb_env` because its numpy/pandas pins would destabilise `R_env`.

# PHASE 4 PROJECT ARCHITECTURE (added 2026-09-03)

## P4.1 The unresolved problem Phase 4 exists to solve

Phase 2 defined the malignant compartment **conservatively and deliberately**. Only
populations with *positive* MPNST/Schwann/neural-crest lineage evidence were collapsed into
`MPNST-Tumor`: 3,420 cells, **17.35%** of 19,716. Phase 2 §P2.3 records the reason — the
prohibited inference was `if (!immune) annotation_ccc <- "MPNST-Tumor"`, which would have
absorbed fibroblasts, pericytes, endothelium and every ambiguous cluster into the tumour.

That decision was correct with the evidence then available, and it is not reversed here. But
it leaves a specific, testable residue:

| Population | Cells | Status entering Phase 4 |
| --- | ---: | --- |
| `MPNST-Tumor` | 3,420 | Called malignant on marker evidence alone |
| `Fibroblast` | 5,064 | Called non-malignant; **MPNST Mes-NC-like malignant cells can carry ECM/mesenchymal programmes** |
| `Candidate-Malignant-Unresolved` | 1,231 | Explicitly provisional; excluded from Phase 3 |
| `Uncertain` | 720 | C12 hypoxic-vs-perineurial, C23 low-complexity; excluded from Phase 3 |
| `Pericyte-VSMC` | 435 | Called non-malignant; a Phase 3 Notch partner |
| `Low-quality-excluded` | 438 | Technical artefact; excluded from Phase 3 |

Marker expression cannot settle this. **Copy-number evidence is orthogonal to it.** That is
the entire rationale for Phase 4.

The uncertainty is not academic: it propagates into every tumour-centric Phase 3 conclusion,
most sharply `MPNST-Tumor ↔ Fibroblast`, `MPNST-Tumor ↔ Pericyte-VSMC` and
`MPNST-Tumor ↔ Endothelial`.

## P4.2 Pipeline architecture

```text
Phase 1   QC / preprocessing
Phase 2   Harmony integration → clustering → hierarchical annotation
Phase 3   sample-aware CCC → multi-method concordance → receiver response → LochNESS

Phase 4   SCEVAN malignancy inference
                ↓
          CNV-based malignant refinement
                ↓
          tumour clone / tumour state resolution
                ↓
          refined CCC tumour annotation
                ↓
          targeted Phase 3 CCC sensitivity analysis
                ↓
          mechanistic tumour-state → TME model
```

Milestones **M28 → M35**.

## P4.3 Phase 4 does NOT invalidate Phase 2 or Phase 3

This must be stated plainly because the temptation to read it the other way is strong.

```text
Phase 2 annotation  =  historical biological annotation, derived from marker evidence
Phase 4 malignancy  =  additional ORTHOGONAL evidence, derived from inferred copy number
```

Phase 2 answered *"what lineage does this cell's transcriptome resemble?"*. Phase 4 answers
*"does this cell carry large-scale copy-number alterations?"*. These are different questions
with different failure modes, and neither is the ground truth for the other. Where they
disagree, **both are reported** (§P4.6).

Phase 2 and Phase 3 outputs are **frozen artefacts**. Phase 4 reads them and writes only into
`results/phase4/` and `reports/phase4/`. No Phase 2 or Phase 3 label is overwritten, moved or
recomputed. Phase 3 conclusions are subjected to a **sensitivity analysis**, not a retraction:
the question is whether the Phase 3 architecture survives a better tumour definition.

## P4.4 Why SCEVAN

De Falco, Caruso, Sudmant & Ceccarelli, *A variational algorithm to detect the clonal copy
number substructure of tumors from scRNA-seq data*, **Nature Communications 14, 1074 (2023)**,
doi:10.1038/s41467-023-36790-9 · https://github.com/AntonioDeFalco/SCEVAN

SCEVAN is appropriate here because it starts from raw counts, identifies confident normal
cells *automatically*, separates malignant from non-malignant cells, infers large-scale CNA
profiles, characterises clones and subclones, supports multi-sample clonal comparison, and is
substantially cheaper than inferCNV in the datasets the original study evaluated. It is also
native R, which matters for a Seurat-centred project.

**SCEVAN is not claimed to be universally superior to inferCNV.** The choice is practical and
scientific: adequate malignant classification + CNA architecture + subclonal capability +
lower computational burden. No second CNV method is added for benchmarking (§P4.9).

## P4.5 Anti-circularity — the rule that makes the analysis mean anything

Two circular arguments are available and both are prohibited.

1. Declaring the disputed fibroblast-like cells to be known normals, then reporting that
   SCEVAN "confirms" they are normal.
2. Declaring all `MPNST-Tumor` cells to be fixed tumour, then reporting that SCEVAN
   "independently confirms" them.

Therefore:

* The **primary** run passes `norm_cell = NULL`, so SCEVAN's own confident-normal detection
  operates with **no input from Phase 2 labels at all**. This is the run that carries
  inferential weight.
* The **sensitivity** run supplies only high-confidence *immune* cells as `norm_cell`, with
  `FIXED_NORMAL_CELLS = FALSE`.
* `FIXED_NORMAL_CELLS = TRUE` is **prohibited in this project**. Reading the SCEVAN source,
  it executes `cellType_pred[!cellType_pred %in% norm_cell_names] <- "malignant"` — every
  cell outside the reference set is forced to malignant. That is precisely the
  non-immune-equals-tumour inference Phase 2 banned, arriving through a function argument.
* `Fibroblast`, `Pericyte-VSMC`, `Candidate-Malignant-Unresolved` and `Uncertain` are **never**
  used as fixed normal references. They are the question.

## P4.6 Disagreement is a result, not a defect

| Phase 2 | SCEVAN | Reading |
| --- | --- | --- |
| `MPNST-Tumor` | malignant | strong concordance |
| `Fibroblast` | normal | supports a true fibroblast |
| `Fibroblast` | malignant | **candidate Mes-NC-like malignant population** |
| `MPNST-Tumor` | normal | investigate — may be copy-number quiet |

Neither method is forced to match the other. `malignancy_refined` integrates Phase 2
annotation, SCEVAN classification, CNA evidence, lineage markers, absence/presence of
canonical immune-stromal identity, patient consistency and clone structure. SCEVAN carries
substantial weight but **not absolute authority**, and every call carries an explicit
`malignancy_confidence` and a printed `decision_reason`. Rules are written to
`reports/phase4/MALIGNANCY_DECISION_RULES.md` before they are applied.

An `Ambiguous` call is a legitimate outcome. Ambiguous cells do **not** become tumour.

## P4.7 Sample-wise SCEVAN is primary

`sample_id = patient = dataset` in this project. SCEVAN therefore runs **separately** for
MPNST_1 … MPNST_4. Pooling four patients first would confuse patient-specific CNV architecture
with subclonal structure within a tumour. Cross-patient comparison happens afterwards, via
`multiSampleComparisonClonalCN()`, and is described as a **shared or recurrent CNA pattern** —
never as "the same clone" across different patients.

## P4.8 Clone ≠ state

```text
tumor_clone_phase4  =  CNV / clonal architecture   (from SCEVAN)
tumor_state_phase4  =  transcriptional phenotype   (from expression)
```

Their relationship is investigated and tabulated; one-to-one correspondence is never assumed.
The Phase 2 Harmony and PCA reductions are frozen — a malignant-only reduction, if
scientifically preferable, is created as a **new** reduction under a new name.

## P4.9 Prohibitions still in force

No spatial analysis · no trajectory inference · no RNA velocity · no survival analysis · no
treatment-response modelling · no deep learning or Transformers · no new pan-cancer analysis ·
no large condition-level DE · **no inferCNV** · **no CopyKAT** — the last two only if a
critical SCEVAN failure made independent CNV validation scientifically necessary, and never
for benchmarking. Phase 5 is not initiated.

## P4.10 Scientific limitations that must travel with every Phase 4 claim

**A.** SCEVAN infers copy number from gene-expression patterns. It is **not DNA sequencing**.
**B.** CNV inference is reliable for broad chromosomal, arm-level and large-segment events;
single-gene CNV calls from scRNA-derived inference are not asserted.
**C.** Some MPNST malignant cells may be relatively copy-number quiet, therefore a **SCEVAN
normal call is not proof of non-malignancy**.
**D.** Tumour state and CNV clone are not equivalent.
**E.** There are **four patients**. No population-level or epidemiological claim follows.
**F.** `sample_id = patient = dataset`, so biological and technical effects cannot be fully
separated — the same confound Phase 2 §P2.5 recorded.
**G.** SCEVAN's annotation covers chromosomes 1–22; sex-chromosome events are not assessable.

## P4.11 Metadata fields Phase 4 adds

Phase 4 **only adds**. Every field below is new; nothing in Phase 1–3 metadata is altered.

```text
scevan_call            malignant / non-malignant / filtered   (SCEVAN, per-sample)
scevan_confident_normal  SCEVAN's own confident-normal flag
scevan_clone           SCEVAN subclone label, patient-scoped
scevan_sample          the SCEVAN run that produced the call
cnv_summary            compact per-cell CNA burden summary

malignancy_phase2      malignant / non-malignant / ambiguous, from annotation_ccc
malignancy_scevan      as called by SCEVAN
malignancy_refined     Malignant / Non-malignant / Ambiguous  (integrated)
malignancy_confidence  High / Moderate / Low
malignancy_reason      printed rule that produced the call

annotation_ccc_phase3  verbatim copy of the frozen Phase 3 annotation_ccc
annotation_ccc_refined Phase 4 CCC layer built on malignancy_refined

tumor_state_phase4     transcriptional state, malignant cells only
tumor_clone_phase4     CNV clone, malignant cells only
```

## P4.12 Environment

Phase 4 runs in `R_env`. **SCEVAN 1.0.3 and yaGST 2017.8.25 were already present**, so no
installation was required and no dependency moved: R 4.4.3, Seurat 5.4.0, SeuratObject 5.3.0,
harmony 1.2.4, Matrix 1.7.4 are unchanged from Phase 2 and Phase 3. Phase 3's CCC stack
(liana 0.1.14, CellChat 2.2.0.9001, CellPhoneDB 5.0.1 in `cpdb_env`, nichenetr 2.2.1.1) is
reused **at identical versions** so that the Phase 3 → Phase 4 comparison is interpretable
rather than confounded by tool drift.

---

# PHASE 5 PROJECT ARCHITECTURE (added 2026-09-03)

## P5.1 The unresolved problem Phase 5 exists to solve

Phase 4 delivered a clean malignant compartment (6,434 cells) and one uncomfortable negative
result: **0 of 8 discrete malignant transcriptional states were recurrent** across ≥3 patients,
and **97.2% of malignant cells sat in patient-private states**. With n = 4 and
`sample_id` = patient = dataset, that is the expected outcome of asking a *clustering* question.

Phase 5 asks a different question, and it is the only one the data can still answer:

> **Do recurrent *continuous* malignant transcriptional programs exist across patients even when
> discrete malignant clusters do not?**

Discrete clusters and continuous programs are not the same object. A cell belongs to exactly one
cluster; a cell carries a *loading* on every program. Two patients can share an ECM program while
their ECM cells cluster separately, because clustering is dominated by patient-specific
covariance that program decomposition can factor out.

**Phase 5 must NOT attempt to rescue the Phase 4 states with another clustering sweep.** That
question was asked, answered, and the answer was negative.

## P5.2 Pipeline architecture

```text
Phase 4 refined malignant cells (6,434)
        ↓  M36  reconstruction, feasibility, historical-field verification
continuous program discovery (cNMF)              ← M37
        ↓
program stability · K selection · patient-balanced sensitivity
        ↓  M38  annotation, pathways, relationship to Phase 4 discrete states
cross-patient recurrence classification
        ↓  M39  malignant ECM-like cells vs true fibroblasts
malignant-ECM signature
        ↓  M40  robustness
        ↓  M41  freeze
Phase 5 malignant-program atlas
```

## P5.3 Primary input, and what is excluded from it

Program discovery uses **only** `malignancy_refined == "Malignant"` (6,434 cells).
`Ambiguous`, `Non-malignant` and `Excluded-low-quality` cells **never define a program**.
`Ambiguous` fibroblasts are **projected** onto already-derived programs afterwards, as a purely
descriptive sensitivity analysis, and are never reclassified.

Expression input is **RNA**, not Harmony coordinates, not UMAP, not PCA scores, and never the
SCEVAN CNA matrix.

## P5.4 The historical annotation field is verified, never assumed

M39 compares *historical* fibroblasts that Phase 4 called malignant against *historical*
fibroblasts that stayed non-malignant. The grouping variable must therefore be the frozen
pre-Phase-4 annotation. `annotation_ccc_refined` is **prohibited** as that grouping variable: it
already encodes the Phase 4 conclusion, so using it would make the comparison circular.

M36 identifies the correct field by testing which frozen annotation column reproduces the frozen
Fibroblast split 5,064 → 4,036 Malignant / 908 Non-malignant / 120 Ambiguous. The expectation is
`annotation_ccc_phase3`, and the expectation is **checked**, not trusted.

## P5.5 Patient imbalance is a first-class analysis object

Malignant-cell abundance is severely unequal across the four patients. Every program is therefore
reported with its cells per patient, per-patient activity, dominant-patient fraction and number of
represented patients, and the full-data primary analysis is accompanied by a **patient-balanced
sensitivity analysis** (equal-cell downsampling). The balanced run is a test of whether the
primary solution is merely the largest patient; it does not replace the primary analysis.

## P5.6 Recurrence is defined before the biology is looked at

Recurrence criteria — minimum cells per patient, minimum fraction of each patient's malignant
cells, an activity threshold, and the number of qualifying patients — are declared **before**
program labels are assigned, and are reported across a grid rather than at a single threshold.
Programs are classified as `recurrent`, `shared-limited`, `patient-private` or `uncertain`.
A program is **not** recurrent merely because a handful of cells in several patients carry
non-zero loadings.

## P5.7 Cell cycle is a candidate program, not a nuisance to be removed

Cell-cycle genes are retained in the primary factorization, because proliferation is a genuine
malignant program. If cell cycle dominates factorization to the point that other programs cannot
resolve, a **declared sensitivity analysis** excluding canonical cell-cycle genes is run and
**both** results are reported. ECM, HLA, interferon, Schwann-lineage, neural-crest and
angiogenesis genes are never removed merely because they dominate a factor — that is the signal.

## P5.8 Statistics

Cells are not biological replicates. **n = 4 patients**, and `sample_id` = patient = dataset.
Conclusions rest on patient-level direction, patient-level recurrence, effect size, cross-patient
consistency and robustness — never on a cell-level p-value computed over thousands of
non-independent cells. Where a group is absent or near-absent in a patient, that patient's
comparison is reported as **NOT EVALUABLE** rather than pooled away.

The malignant-ECM signature is a **transparent descriptive signature**. No high-capacity
classifier is trained, and no random cell-level train/test split is used to report
pseudo-independent accuracy.

## P5.9 Interpretation language

A transcriptional program is a continuous pattern of co-regulated expression. It is **not**
automatically a cell type, a lineage, a state or a clone. `program`, `state`, `clone` and
`cell identity` remain four separate concepts throughout Phase 5 and Phase 6.

## P5.10 Prohibitions still in force

Everything prohibited in Phases 1–4 remains prohibited. Phase 5 adds: no new CCC method and no
rerun of LIANA, CellChat, CellPhoneDB, NicheNet or LochNESS; no RNA velocity; no pseudotime
without a validated biological question; no survival analysis; no condition DE; no patient-level
predictive ML; no deep learning; no spatial, pan-cancer or other-sarcoma comparison; no second
CNV method. Phase 4 is not modified in any respect.

## P5.11 Environment

`R_env` is **not** modified. cNMF is a Python tool and therefore lives in an isolated
`p5_cnmf_env`; nothing else in the project depends on it. Gene-set resources (MSigDB Hallmark,
Reactome, GO BP) are downloaded once to `external/genesets/` with recorded checksums rather than
installed as R packages, so the frozen stack stays frozen.

---

# PHASE 6 PROJECT ARCHITECTURE (added 2026-09-03)

## P6.1 Biological objective

```text
Phase 5 malignant programs  +  Phase 4 SCEVAN CNA clones
        ↓
clone ↔ program coupling                          ← M43
within-clone program diversity                    ← M44
TF / regulatory activity                          ← M45
pathway activity                                  ← M46
broad CNA → transcriptional consequences          ← M47
integrated malignant architecture                 ← M48
robustness                                        ← M49
freeze                                            ← M50
```

Central questions:

> Is malignant transcriptional heterogeneity primarily associated with distinct CNA-defined
> clones, with substantial transcriptional-program diversity *within* clones, or with a mixture?

> Which transcription-factor and pathway activities distinguish the recurrent malignant programs?

## P6.2 MPNST_3 is excluded from every clone-based inference

Phase 4 established that MPNST_3 failed the SCEVAN immune sanity gate: its two runs are inverted
(agreement 0.0864), six of nine canonical immune populations are called ~100% malignant, and its
three "clones" are T/NK, plasma/pDC/B and plasma-dominated. Its clone structure therefore carries
**no tumour-clone inferential weight**.

Primary Phase 6 clone analyses use **MPNST_1, MPNST_2, MPNST_4 only**. MPNST_3 appears solely in
explicitly labelled QC panels showing *why* it is excluded, and its Phase 5 transcriptional
programs may still be described — but never its clone architecture.

## P6.3 Clone labels are patient-scoped and never homologous

`MPNST_1_clone1` and `MPNST_4_clone1` are unrelated labels. Clone-program analyses are performed
**within each patient independently** and are never pooled as if clone numbers matched across
patients.

## P6.4 Plasticity wording

Cross-sectional scRNA-seq cannot demonstrate state transitions. The permitted primary description
is **within-clone transcriptional-program diversity**. It may be interpreted as *consistent with
transcriptional plasticity* only after technical confounders — clone size, library size, detected
genes, cell cycle, malignancy confidence — have been assessed. Observed state switching, transition
rates, directional transitions and lineage trajectories are **not claimable** from these data.

## P6.5 The CNA → expression analysis is not orthogonal validation

SCEVAN infers copy number *from expression*. Any association between an inferred broad CNA event
and the expression of genes on the affected chromosome is therefore an **internal transcriptional
consistency analysis**, not independent confirmation of the copy-number call. This limitation
travels with every Phase 6 CNA statement.

CNA claims are permitted only at chromosome, chromosome-arm and large-segment level. Permitted:
"cells assigned to a clone carrying an inferred broad chr17 loss whose interval contains NF1".
Prohibited: "NF1-deleted cells", "NF1 deletion confirmed", "EGFR amplification confirmed".

## P6.6 Regulatory and pathway inference

TF activity uses **decoupleR** with **CollecTRI** (and DoRothEA where compatible); pathway activity
uses **PROGENy** through decoupleR plus Hallmark gene-set scoring. Input is gene expression —
never Harmony, UMAP or SCEVAN CNA values. pySCENIC is **not** introduced merely because it is
popular. Program, TF activity and pathway activity are carried as **three independent evidence
layers**; they are never collapsed into an invented weighted composite score or a single numeric
rank.

A TF–program relationship is not recurrent merely because a pooled cell-level association is
strong: the association is reported within each patient, with the number of patients showing a
concordant direction and the dominant-patient fraction.

## P6.7 Model selection is not decided in advance

Three architectures are evaluated and whichever the data support is reported:

* **Model A** — clone-constrained transcriptional states;
* **Model B** — extensive within-clone program diversity, consistent with plasticity;
* **Model C** — mixed architecture.

With n = 4, deep per-patient tumour portraits are the honest deliverable, not population-level
statistical claims.

## P6.8 Environment

`R_env`, unchanged. decoupleR 2.12.0 and OmnipathR 3.14.0 are **already installed**, so CollecTRI
and PROGENy are retrieved through the existing stack and cached to disk with recorded provenance
for reproducibility.
