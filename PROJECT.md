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
