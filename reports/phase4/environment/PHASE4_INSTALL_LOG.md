# Phase 4 — Environment and Installation Log

**Bottom line: the frozen Phase 2/3 R stack did not move.** Verified before and after every
operation below.

| Package | Phase 2/3 | After all Phase 4 setup | Moved? |
| --- | --- | --- | :--: |
| R | 4.4.3 (2025-02-28) | 4.4.3 (2025-02-28) | no |
| Seurat | 5.4.0 | 5.4.0 | no |
| SeuratObject | 5.3.0 | 5.3.0 | no |
| Matrix | 1.7.4 | 1.7.4 | no |
| harmony | 1.2.4 | 1.2.4 | no |
| sctransform | 0.4.3 | 0.4.3 | no |

Environment snapshots: `R_env_PRE_PHASE4.yml` (671 lines, captured before anything was
touched) · `p4_umap_env_PHASE4.yml`.

---

## 1. SCEVAN was already installed — no installation attempted

```text
SCEVAN   1.0.3   (Date: 2025-02-12; built for R 4.4.3 on 2025-11-14)
yaGST    2017.8.25
umap     0.2.10.0   (R package)
reticulate 1.44.1
```

`requireNamespace("SCEVAN")` succeeded on first inspection, so **no `remotes::install_github`
call was made and no R dependency was resolved, upgraded or downgraded** for the CNV work
itself. This is the reason the table above shows no movement.

Exported API confirmed against the installed version rather than assumed from documentation:

```text
annotateGenes  annoteBandOncoHeat  classifyTumorCells  getBreaksVegaMC
getConfidentNormalCells  multiSampleComparisonClonalCN  pipelineCNA
plotAllClonalCN  plotAllSubclonalCN  plotCNA_withAnnotCells  top30classification
```

`pipelineCNA()` signature as installed:

```r
pipelineCNA(count_mtx, sample = "", par_cores = 20, norm_cell = NULL,
            SUBCLONES = TRUE, beta_vega = 0.5, ClonalCN = TRUE, plotTree = TRUE,
            AdditionalGeneSets = NULL, SCEVANsignatures = TRUE, organism = "human",
            ngenes_chr = 5, perc_genes = 10, FIXED_NORMAL_CELLS = FALSE,
            output_dir = "./output", orig.ident = NULL)
```

---

## 2. One addition: Python `umap-learn`, in an isolated env

### Why it was needed

`SCEVAN:::subcloneAnalysisPipeline` calls `plotTSNE()`, which calls
`umap::umap(..., method = "umap-learn")` — the **Python** implementation, reached through
reticulate. Reading the source established that this is **not** an optional cosmetic step:

```r
plotOncoHeatSubclones(...)
plotTSNE(count_mtx, res_class$CNAmat, ...)                      # <- line 93, fails without umap-learn
classDf[names(res_subclones$clustersSub), "subclone"] <- ...    # <- line 97, writes the CLONE LABELS
```

`plotTSNE()` runs *before* the line that writes subclone assignments into `classDf`, so its
failure loses the clone labels entirely. Subclonal structure is a Phase 4 deliverable (§18,
§32, §33), so the dependency had to be satisfied.

### Why an isolated environment

Installing into `R_env` risks the conda solver moving R packages, and Phase 2/3 froze that
stack. A Python-only environment cannot touch it. This follows the Phase 3 precedent of
isolating CellPhoneDB in `cpdb_env`.

```bash
conda create -y -n p4_umap_env -c conda-forge python=3.11 "umap-learn>=0.5" \
      numpy scipy scikit-learn numba pynndescent
```

```text
python       3.11
umap-learn   0.5.12
numpy        2.4.6
scikit-learn 1.9.0
numba        0.67.0
```

Reached from R by exporting `RETICULATE_PYTHON` in the Phase 4 SLURM scripts:

```bash
export RETICULATE_PYTHON=/local/.../anaconda3/envs/p4_umap_env/bin/python
```

Verified by executing the exact call SCEVAN makes (`umap(..., method = "umap-learn",
n_components = 2, n_neighbors = 15, min_dist = 0.1, metric = "euclidean", seed = 1)`) and
re-printing every `R_env` version afterwards. **`R_env` was not modified.**

**Dependency changes to `R_env` in Phase 4: none.**

---

## 3. Two SCEVAN 1.0.3 defects found and worked around without patching the package

### 3.1 `output_dir` is not honoured by the read-back and plotting helpers

`pipelineCNA()` accepts `output_dir` and forwards it to most helpers, but these hardcode
`path = "./output"` and their callers do not forward it:

```r
getScevanCNV(sample, path = "./output", ...)
getScevanCNVfinal(sample, path = "./output", ...)
plotAllClonalCN   plotAllSubclonalCN   plotConsensusCNA
plotCNA_withAnnotCells   analyzeSegm2
# and: plotCNclonal() calls getScevanCNV(paste0(sample, name))  -- no output_dir
```

**Symptom** (SLURM 19896576): classification completed and reported `"found 363 tumor cells"`,
then the pipeline died with `cannot open the connection`.

**Fix**: each SCEVAN run sets the *process working directory* to its own folder and lets
`output_dir` take SCEVAN's default `"./output"`, so the parameterised and hardcoded paths refer
to the same place. The package is not edited and no error is suppressed.

### 3.2 `plotTSNE` needs Python umap-learn

Root-caused in §2 above (SLURM 19896623: `Python module umap was not found`).

---

## 4. A shell defect of ours, for completeness

SLURM 19896654 failed in 11 s with
`qt-main_activate.sh: line 5: QT_XCB_GL_INTEGRATION: unbound variable`. Cause: `set -u`
combined with conda's `qt-main` activation hook, which dereferences an unset variable —
triggered by `conda deactivate`. Conda's own guidance is to disable `nounset` around
activation, so the Phase 4 shell scripts use `set -o pipefail` without `-u`, with the reason
written in the script.

---

## 5. What was NOT done

* No inferCNV, no CopyKAT — no second CNV method was added for benchmarking (§66).
* No Phase 2 or Phase 3 package was downgraded to make a Phase 4 tool work.
* No SCEVAN source file was edited.
* No resource request was raised in response to any failure above. M28's 96 G was set from the
  Phase 3 measured load (27.4 GiB) and used 30.1 GiB at peak; it was never increased.

---

## 6. A third SCEVAN 1.0.3 defect, found during M29 — the clone phylogeny plot

Every M29 run emitted this, caught by SCEVAN's own internal `tryCatch`:

```text
Error in `stat_tree()`: Problem while converting geom to grob.
Caused by error in `is.waive()`: could not find function "is.waive"
  SCEVAN:::subcloneAnalysisPipeline(...)
    SCEVAN:::plotCloneTree(sample, res_subclones, output_dir = output_dir)
      ggplot2::ggsave(...) -> ggtree stat_tree
```

`ggtree` calls `ggplot2:::is.waive()`, which **ggplot2 4.x no longer exports**. The environment
here has ggplot2 4.0.3 because Phase 3 required it.

**Impact: cosmetic only.** SCEVAN wraps `plotCloneTree` in its own `tryCatch`, so the run
continues and `pipelineCNA` returns normally. Every scientific output is intact — subclone
assignments in `classDf`, `*_Clonal_CN.seg`, `*_subcloneN_CN.seg`, `*_CNAmtx.RData`,
`*_CNAmtxSubclones.RData`, `*PlotOncoHeat.RData`, the CNA heatmaps, the consensus profile and
the onco-heatmaps. Only the clone **phylogeny** figure is missing.

**Not fixed, deliberately.** The fix would be downgrading ggplot2 below 4.0, which would
destabilise the Phase 3 figure suite and violate the rule against downgrading a frozen
dependency to make a tool work. A clone dendrogram is not a Phase 4 deliverable, and the
clone-versus-state relationship is reported instead as an explicit cross-tabulation
(`CLONE_VS_TUMOR_STATE.tsv`, figure `31_04_clone_vs_tumor_state.pdf`), which is more
interpretable than a tree anyway. Recorded as Phase 4 limitation K.

---

## 7. Reproduction

```bash
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh

# one-off: the isolated interpreter SCEVAN's subclone stage needs
sbatch scripts/shell/phase4/install_umap_learn.sh

# M28 -> M29 (per patient) -> M30 -> M32 -> M31/M33/M34 -> M35
J28=$(sbatch --parsable scripts/shell/phase4/run_m28_feasibility.sh)
J29=$(sbatch --parsable --dependency=afterok:$J28 scripts/shell/phase4/run_m29_scevan.sh)
J30=$(sbatch --parsable --dependency=afterok:$J29 scripts/shell/phase4/run_m30.sh)
JP=$(sbatch  --parsable --dependency=afterok:$J30 scripts/shell/phase4/run_m32_prepare.sh)
JL=$(sbatch  --parsable --dependency=afterok:$JP  scripts/shell/phase4/run_m32_liana.sh)
JC=$(sbatch  --parsable --dependency=afterok:$JP  scripts/shell/phase4/run_m32_cellchat.sh)
JD=$(sbatch  --parsable --dependency=afterok:$JP  scripts/shell/phase4/run_m32_cellphonedb.sh)
JK=$(sbatch  --parsable --dependency=afterok:$JL:$JC:$JD scripts/shell/phase4/run_m32_concordance.sh)
JT=$(sbatch  --parsable --dependency=afterok:$JK  scripts/shell/phase4/run_m31_m33_m34.sh)
sbatch --dependency=afterok:$JT scripts/shell/phase4/run_m35.sh
```

Every script sets `RETICULATE_PYTHON` to the isolated interpreter and matches `par_cores` to
`--cpus-per-task`. Seeds are fixed at 42 throughout (SCEVAN, clustering, UMAP, module scores,
subsampling).
