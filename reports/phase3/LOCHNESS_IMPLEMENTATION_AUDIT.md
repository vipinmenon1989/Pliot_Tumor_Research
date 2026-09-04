# LochNESS Implementation Audit

**Phase 3 · Milestone M18 · MPNST project**
*Generated: 2026-09-03*

Three sources were compared, as required before any LochNESS code is written for MPNST:

| Source | Location |
| --- | --- |
| **Publication** | Huang et al., *Single-cell, whole-embryo phenotyping of mammalian developmental disorders*, **Nature 623, 772–781 (2023)**, doi **10.1038/s41586-023-06548-w** |
| **Official R implementation** | `https://github.com/shendurelab/MMCA` — vendored read-only at `external/MMCA_ref/` (commit `af629c498f421d6f0bfb325e9fd6b2dab0e22837`, 2024-06-12, Chengxiang Qiu) |
| **Our Python implementation** | `PertTF-Virtual-Challeng-Weilab/perturbseq-pipeline/src/perturbseq_pipeline/lochness.py` (read-only reference; ported from pertTF `composition_change_analysis`) |

Files vendored from MMCA: `demo_lochness.R`, `Section_5_step_1_run_lochness_calculation.R`,
`Section_5_step_2_run_lochness_calculation_permutation.R`.

> **MMCA is a repository of analysis scripts, not an installable R package.** It has no
> `DESCRIPTION`, no `R/` directory and no namespace. `remotes::install_github("shendurelab/MMCA")`
> would therefore fail, and was **not** attempted. The scripts were vendored as a pinned
> read-only reference instead, and the MPNST implementation is written separately.

---

## 1. Mathematical definition

The published definition, confirmed identically in both implementations:

```
lochNESS(cell i, group g) = local_fraction_i(g) / global_fraction(g) − 1
```

- `0` → group `g` appears in cell *i*'s neighbourhood exactly as often as chance predicts
- `> 0` → locally over-represented
- `< 0` → locally under-represented

The official R code expresses it as (`demo_lochness.R`, `Section_5_step_1`):

```r
kadj = round(0.5 * sqrt(ncol(seurat_object)))
...
mutant_mask      = ifelse(knn$label == "WT", 0, 1)
mutant_neighbors = rowSums(mutant_mask)
pca_mutant_score = (mutant_neighbors / kadj) /
                   (1 - mt_counts$Freq[mt_counts$Var1=='WT'] / sum(mt_counts$Freq)) - 1
```

The denominator `1 − (n_WT / n_total)` is the global fraction of **non-WT (mutant)** cells,
computed **over the reference set only** (the cells that were eligible to be neighbours), not
over the whole dataset. This is a subtle but important detail.

The Python implementation states the same formula in its module docstring:

```
lochNESS(cell, g) = local_fraction(g) / overall_fraction(g) - 1
```

**The two agree on the definition. They differ in six implementation choices, below.**

---

## 2. Neighbourhood construction

| Aspect | MMCA R (official) | perturbseq Python |
| --- | --- | --- |
| Representation | **L2-normalised PCA** — `Seurat::L2Dim(reduction='pca')`, using `pca.l2` embeddings | `X_pca_harmony` if present, else `X_pca`; **no L2 normalisation** |
| Neighbour engine | `FNN::get.knnx(..., algorithm="kd_tree")` — **exact** kNN | `scanpy.pp.neighbors` — approximate (UMAP/annoy backend) |
| Graph reuse | none; kNN recomputed per query group | reuses a dedicated `lochness_nn` graph when present |
| Self inclusion | structurally impossible (query and reference sets are disjoint) | self excluded by scanpy convention |

## 3. Value of *k*

| MMCA R | perturbseq Python |
| --- | --- |
| **`kadj = round(0.5 * sqrt(N))`** where `N` = cells in the analysis subset. Data-adaptive. | **`k = min(cfg.lochness.n_neighbors, N-1)`** — a configured constant (pertTF default ≈ 300). Not data-adaptive. |

This is a real divergence. On a 19,716-cell object the official rule gives
`round(0.5·√19716) = 70`; the Python default would use ~300.

## 4. Reference population and same-sample exclusion — **the most important difference**

**MMCA R excludes same-embryo neighbours by construction.** The loop is over embryos
(`RT_group`); for each embryo *j* the query cells are that embryo's cells and the reference
matrix contains **only cells from other embryos**:

```r
for (j in unique(seurat_object$RT_group)) {
  mt_cells = colnames(seurat_object)[seurat_object$RT_group == j]   # query
  wt_cells = colnames(seurat_object)[seurat_object$RT_group != j]   # reference: OTHER embryos only
  knn = get.knnx(t(wt.mtx), t(mt.mtx), k = kadj, algorithm = "kd_tree")
}
```

Because the reference matrix never contains the query embryo's own cells, no neighbour can
come from the same embryo. Output files are named `..._lochness_excludeself_pool.rds`,
confirming the intent.

**perturbseq Python implements no such exclusion.** Grepping `lochness.py` for
`same_sample|same_embryo|exclude|batch_key|sample_key` returns only unrelated matches (a
docstring about self-exclusion in the scanpy graph, and a note about targets excluded for
low cell counts). Neighbours are drawn from the whole dataset.

**Consequence for MPNST.** Without exclusion, a per-sample score largely re-detects
within-sample structure — and in this dataset `sample_id` is simultaneously the dataset, the
patient and the batch (Phase 2 M12 measured 78% same-sample nearest neighbours even *after*
Harmony). Same-sample exclusion is therefore not a stylistic detail here; omitting it would
make the score a batch detector.

## 5. Scope of scores produced

| MMCA R | perturbseq Python |
| --- | --- |
| One score per **query cell** for the single contrast "is my neighbourhood enriched for mutant cells". Query cells are visited embryo by embryo. | Two modes: a full `cell × group` matrix, or (for very large data) `lochness_self` — each cell's score for its *own* label. |
| `demo_lochness.R` also offers an **alternative implementation**: compute `lochNESS_i` once per wildtype reference sample, then examine the correlation matrix across those scores to assess how much the result depends on which reference sample was used. | no equivalent |

That alternative is effectively a leave-one-reference-sample-in sensitivity analysis and is
directly adaptable to MPNST.

## 6. Normalisation

Neither implementation z-scores or rank-transforms the output. The only normalisation is the
division by the global fraction, which is intrinsic to the definition. The R version
additionally L2-normalises the PCA embedding *before* computing distances (§2), which the
Python version does not.

Minor: the R version divides the local count by the requested `kadj`; the Python version
divides by the *actual* neighbour count per row (its docstring notes the scanpy graph stores
`n_neighbors − 1` entries). At k = 70 that is a ~1.4% systematic difference in the local
fraction; at k = 300, ~0.3%.

## 7. Permutation / null strategy

`Section_5_step_2_run_lochness_calculation_permutation.R`:

```r
for (shuffle_ind in c(1:nn)) {           # nn = 10
  set.seed(shuffle_ind + 2022)
  seurat_object$sMutant = as.vector(seurat_object$Mutant[sample(length(seurat_object$Mutant))])
  ...                                     # identical computation using sMutant
}
```

- **The group label is permuted, not the cells or the embedding.** 10 shuffles.
- **The neighbourhood structure and the same-embryo exclusion loop are preserved** — the
  `RT_group` partition is untouched, so the geometry and exclusion logic are identical
  between observed and null.
- The shuffle is *global across cells* in the analysis subset, so the null destroys any
  association between label and sample. That is the intended null for a mutant-vs-WT design,
  but it means the null does **not** preserve per-sample label composition.
- **perturbseq Python implements no permutation null at all** (`rng` appears only for
  subsampling; `random_state` only for the neighbour graph).

## 8. Assumptions of the original method

1. There is a meaningful **binary-ish contrast** (mutant vs wildtype) with a designated
   **reference population** that is biologically "normal".
2. Multiple **independent replicates of the reference** exist, so that a query sample can be
   scored against other samples.
3. The **shared embedding is comparable across samples**, i.e. batch effects have already
   been handled, otherwise the score measures batch.
4. Cells are embedded in a **single continuous manifold** in which "local neighbourhood" is
   biologically meaningful (in MMCA, developmental trajectories).
5. The group of interest is a **minority** of the reference-set composition (the denominator
   `1 − n_WT/n_total` is small), which is what makes the ratio informative.

## 9. What can legitimately transfer to MPNST

| Transferable | Why |
| --- | --- |
| The **formula** `local/global − 1` | Definition is design-agnostic. |
| **`k = round(0.5·√N)`** | Data-adaptive and published; adopt rather than the Python constant. |
| **L2-normalised PCA-space kNN with exact `FNN::get.knnx`** | Published choice; exact kNN is cheap at ~3k–5k cells per receiver lineage. |
| **Same-sample exclusion by disjoint query/reference sets** | Essential here, because sample = patient = batch (Phase 2 M12). |
| **Global fraction computed over the reference set** | Needed for the score to be centred at 0 under the null. |
| **Label-permutation null with preserved neighbourhood/exclusion structure** | Directly adoptable; the Python version has no null. |
| **The `lochNESS_i` per-reference-sample variant** | Becomes a leave-one-sample-out robustness check, which Phase 3 M26 needs anyway. |

## 10. What cannot transfer

| Not transferable | Why |
| --- | --- |
| **"Mutant vs WT" group semantics** | MPNST has no mutant/wildtype axis and **no biological-condition variable at all** (Phase 2 M10 audited all 93 metadata columns). A target/reference contrast must be defined afresh — see `reports/phase3/LOCHNESS_MPNST_DESIGN.md`. |
| **The interpretation "genotype drives local enrichment"** | There is no genotype contrast. Any MPNST enrichment is between samples/contexts, which are confounded with patient and batch. |
| **A globally-shuffled label null, used naively** | With only four samples and label perfectly confounded with sample, a global shuffle would destroy the sample structure entirely. Phase 3 uses a **sample-aware** permutation instead (M24), and the deviation is documented. |
| **Trajectory-conditioned analysis** | MMCA computes LochNESS separately per developmental trajectory. MPNST has no trajectory (and trajectory inference is a Phase 3 prohibition). The analogous stratification here is **per receiver lineage**. |
| **The Python `lochness_self` large-mode** | Unnecessary at this scale (largest receiver lineage ≈ 5k cells) and it omits the exclusion logic that matters most here. |

## 11. Implementation decision

**Phase 3 will implement LochNESS natively in R**, following the MMCA formulation
(`0.5·√N`, `L2Dim` PCA, exact `FNN::get.knnx`, disjoint query/reference sets, global fraction
over the reference set, label-permutation null with preserved exclusion structure), because:

1. The MPNST workflow is Seurat/R-centred, so an R implementation avoids an
   AnnData conversion round-trip and a second environment.
2. The official implementation *is* R, so following it reduces the chance of silent drift.
3. Scale is small — the largest receiver lineage is ~5,064 cells, so exact kNN at
   `k = round(0.5·√5064) = 36` is trivial. The Python version's scalability advantages
   (Numba kernel, million-cell mode) solve a problem MPNST does not have.

The Python implementation is retained as an **independent cross-check**: M24 runs both on the
same deterministic subset and reports the agreement in
`results/phase3/lochness/LOCHNESS_IMPLEMENTATION_COMPARISON.tsv` and
`reports/phase3/LOCHNESS_IMPLEMENTATION_COMPARISON.md`.

## 12. Scientific guard carried forward

LochNESS is **not** a ligand–receptor inference method. Phase 3 will never write or imply
"high LochNESS = cell–cell communication". Its role is strictly as an **orthogonal
transcriptional-neighbourhood phenotype** used to characterise receiver-state structure
*after* CCC has nominated candidate signalling, per the Phase 3 specification §6.
