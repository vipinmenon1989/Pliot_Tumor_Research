# Phase 3 Installation Log

*Generated 2026-09-03 · Milestone M18*

Pre-Phase-3 snapshot: `R_env_PRE_PHASE3.yml` (525 lines, taken **before** any change).
Post-Phase-3 snapshots: `R_env_PHASE3.yml` (671 lines), `cpdb_env_PHASE3.yml`.

## Baseline that must not change — verified after every stage

| Package | Before | After | Status |
| --- | --- | --- | --- |
| R | 4.4.3 | 4.4.3 | ✅ unchanged (`r-base` explicitly pinned `==4.4.3` in the conda calls) |
| Seurat | 5.4.0 | 5.4.0 | ✅ |
| SeuratObject | 5.3.0 | 5.3.0 | ✅ |
| harmony | 1.2.4 | 1.2.4 | ✅ |
| Matrix | 1.7.4 | 1.7.4 | ✅ |
| sctransform | 0.4.3 | 0.4.3 | ✅ |
| Python (R_env) | 3.11.13 | 3.11.13 | ✅ |

## Installed for Phase 3

| Package | Version | Source | Command | Reason | Result |
| --- | --- | --- | --- | --- | --- |
| **liana** | 0.1.14 | GitHub `saezlab/liana` | `remotes::install_github("saezlab/liana", upgrade="never")` | Primary consensus LR framework | ✅ |
| OmnipathR | 3.14.0 | Bioconductor | `BiocManager::install` | liana's LR resource backend | ✅ |
| decoupleR | 2.12.0 | Bioconductor | `BiocManager::install` | liana dependency | ✅ |
| basilisk | 1.18.0 | Bioconductor | `BiocManager::install` | liana optional Python bridge | ✅ |
| **CellChat** | 2.2.0.9001 | GitHub `jinworks/CellChat` | `remotes::install_github("jinworks/CellChat", upgrade="never")` | Independent pathway/network framework | ✅ (after 2 retries) |
| **nichenetr** | 2.2.1.1 | GitHub `saeyslab/nichenetr` | `remotes::install_github("saeyslab/nichenetr", upgrade="never")` | Receiver-response modelling | ✅ (after 2 retries) |
| **cellphonedb** | 5.0.1 | PyPI, **isolated env `cpdb_env`** | `conda create -n cpdb_env python=3.10; pip install cellphonedb` | Independent permutation LR framework | ✅ |
| anndata / scanpy (cpdb_env) | 0.11.4 / 1.11.5 | PyPI dependency | — | CellPhoneDB input handling | ✅ |
| NMF | 0.28 | CRAN | `install.packages` | CellChat dependency | ✅ |
| ggraph / tidygraph | 2.2.2 / 1.3.0 | conda-forge | `conda install -c conda-forge` | network figures | ✅ |
| ggpubr / rstatix / car / doBy / pbkrtest / Deriv | 1.0.0 / 1.1.0 / 3.1.5 / 4.7.2 / 0.5.5 / 4.2.0 | conda-forge | `conda install -c conda-forge --freeze-installed` | CellChat + nichenetr dependency chain | ✅ |
| systemfonts | **1.3.2** | conda-forge | `conda install -c conda-forge "r-systemfonts>=1.3.2"` | svglite 2.2.2 requires ≥ 1.3.0 | ✅ |
| svglite | 2.2.2 | conda-forge | — | CellChat dependency | ✅ |
| shadowtext / units / sf / ggiraph / gdtools | 0.1.6 / 1.0.1 / 1.0.21 / 0.9.6 / 0.5.0 | conda-forge | `conda install -c conda-forge` | nichenetr dependency chain (needs udunits2, GDAL/PROJ/GEOS system libraries) | ✅ |
| fdrtool / ROCR / caTools / Hmisc | 1.2.18 / 1.0.12 / 1.18.4 / 5.2.6 | conda-forge | — | nichenetr dependencies | ✅ |
| mlrMBO / ParamHelpers / smoof / mlr / lhs / emoa | (CRAN latest) | CRAN | `install.packages` | nichenetr dependencies | ✅ |
| caret / randomForest / e1071 / DiagrammeR | 7.0.1 / 4.7.1.2 / 1.7.17 / 1.0.12 | conda-forge | — | nichenetr dependencies | ✅ |
| BiocNeighbors, SingleCellExperiment, SummarizedExperiment, scater, scran, limma, edgeR, AnnotationDbi, org.Hs.eg.db | Bioconductor 3.20 set | Bioconductor | `BiocManager::install` | SCE interoperability, annotation | ✅ |

## Dependency change to record honestly

| Package | Before | After | Assessment |
| --- | --- | --- | --- |
| **igraph** | 2.2.1 | **2.1.4** | **Downgraded** by the conda solver while satisfying the CellChat/nichenetr chain. Phase 2 is frozen and all its outputs exist on disk as validated artefacts, so no Phase 2 result changes. The risk is confined to *re-running* Phase 2 M13 clustering, which is not planned. Recorded rather than hidden. |

No other existing package changed version.

## Failures, root causes and fixes

### Attempt 1 — JobID 19894761 (COMPLETED, 00:43:55)

liana, OmnipathR, decoupleR, basilisk, NMF and the MMCA vendoring all succeeded.
CellPhoneDB installed in `cpdb_env`. **CellChat and nichenetr FAILED.**

**Root cause:** source compilation of `systemfonts` and `Deriv` failed, cascading to
`svglite`, `doBy`, `pbkrtest`, `car`, `rstatix`, `ggpubr` — all hard dependencies of
CellChat, and of nichenetr via its plotting stack.

### Attempt 2 — JobID 19895124 (COMPLETED, 00:04:45)

**Fix applied:** install the failing compile chain as **prebuilt conda-forge R 4.4 binaries**
instead of compiling from source, with `--freeze-installed` so no existing package could be
altered. 14 packages installed successfully.

**Still failed.** CellChat: `namespace 'systemfonts' 1.2.3 is being loaded, but >= 1.3.0 is
required` — conda's frozen solve had picked systemfonts 1.2.3 while svglite 2.2.2 needs
≥ 1.3.0. nichenetr: `dependency 'shadowtext' is not available`, plus
`configuration failed for package 'units'` (missing udunits2) and `sf` (missing GDAL/PROJ/GEOS).

### Attempt 3 — JobID 19895459 (COMPLETED, 00:11:02)

**Fix applied:** `conda install -c conda-forge "r-base==4.4.3" "r-systemfonts>=1.3.2"
r-shadowtext r-units r-sf r-ggiraph r-gdtools r-fdrtool r-rocr r-catools r-hmisc`, with
`r-base` pinned so the interpreter could not move under the frozen Phase 2 stack. Then
retried both GitHub installs.

**Result: CellChat 2.2.0.9001 ✅ and nichenetr 2.2.1.1 ✅.**

### CellPhoneDB API fix — JobIDs 19895344 (failed) → 19895438 (succeeded)

`db_utils.download_database(target_dir, version=...)` raised
`TypeError: unexpected keyword argument 'version'`. **Root cause:** CellPhoneDB 5.0.1's
signature is `download_database(target_dir, cpdb_version)`. Fixed the keyword; the retry
completed in 00:02:14.

## Discipline observed

- **No foundational package was upgraded or downgraded on purpose.** `r-base` was explicitly
  pinned; Seurat, SeuratObject, harmony, Matrix and sctransform were verified byte-version
  identical after every stage.
- **No failure was addressed by increasing memory or walltime.** Every failure was diagnosed
  from the log and fixed at its root cause.
- **CellPhoneDB was isolated on purpose.** Its numpy/pandas pins would have moved `R_env`'s
  numpy 2.3.5, which `reticulate` and other Phase 2 components sit on. `R_env` was left
  untouched and the isolation is recorded in the method plan.
- **MMCA was vendored, not installed.** It is a repository of scripts with no `DESCRIPTION`,
  so `remotes::install_github()` would fail; `external/MMCA_ref/` holds a pinned read-only
  copy at commit `af629c498f421d6f0bfb325e9fd6b2dab0e22837`.

## Compute-node network verification

JobID 19894719 on `ihc-grid-1-1-1` confirmed HTTP 200 from CRAN, GitHub API, Bioconductor,
PyPI and Zenodo, so all installation ran inside SLURM rather than on the login node.
