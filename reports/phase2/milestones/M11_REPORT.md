# Milestone M11 — Default Harmony Integration

**Phase 2 · MPNST single-cell integration**
*Generated: 2026-09-02*
*Status: **COMPLETE — STOPPED, awaiting authorization for M12***

> **SCIENTIFIC STATUS:** Harmony execution has completed, but integration quality has
> **not** yet been scientifically accepted. Pre/post-Harmony evaluation is reserved for
> M12. Nothing in this report asserts that the integration is biologically correct.

---

## 1. Input

| Property | Value |
| --- | --- |
| Phase 1 handoff object | `results/combined/pre_integration/combined_preintegration.rds` |
| Size | 6,041,513,977 bytes (6.04 GB) |
| md5 (verified at runtime) | `88a442688f912d882f6c6da01820e329` ✅ matches `prov_analysis.json` |
| sha256 (verified at runtime) | `c3fdce8b61602725977f989d8bbf10030ffe48ef6140372a13d2c59903a8dc66` ✅ matches `phase1_manifest.json` |
| Cells | 19,716 |
| Features (default assay) | 29,113 (`SCT`) |
| Features (`RNA` assay) | 31,764 |
| Default assay | `SCT` |
| PCA reduction consumed | `pca` — 19,716 cells × 50 dims, assay `SCT` |
| Phase 1 UMAP present | `umap_preintegration` |
| Phase 1 graphs present | `SCT_nn`, `SCT_snn` |
| Metadata columns | 93 |

The object was opened **read-only**. Both Phase 1 checksums were recomputed on the
compute node before loading and asserted equal; the job would have aborted on any
mismatch.

---

## 2. Harmony

### 2.1 Grouping variable

**`sample_id`** — 4 levels, sole entry in `group.by.vars`.

| `sample_id` | Cells | Share |
| --- | ---: | ---: |
| MPNST_1 | 7,615 | 38.62% |
| MPNST_2 | 2,284 | 11.58% |
| MPNST_3 | 2,940 | 14.91% |
| MPNST_4 | 6,877 | 34.88% |

Verified at runtime: 4 levels, no `NA`, per-level counts identical to the M10 record,
and cell-for-cell identical to `orig.ident` (`TRUE`).

### 2.2 Biological justification (inherited from M10)

From `reports/phase2/HARMONY_VARIABLE_DECISION.md`, approved at the M10 gate:

- `sample_id` is the only metadata field in the object that carries technical batch
  structure at all. `orig.ident` is identical to it; every other candidate is prohibited
  (legacy annotation, legacy clusters, legacy Harmony/CCA/RPCA/MNN output), degenerate
  (zero-variance QC flags), continuous (not a Harmony grouping variable), or non-existent
  (no clinical or technical covariate of any kind is recorded).
- **Dataset = sample = patient = presumed technical batch is a single 4-level variable.**
  Correcting technical batch is mathematically indistinguishable from erasing
  between-patient tumour biology. The researcher reviewed and approved this trade-off
  at the M10 gate; the risk analysis in `HARMONY_VARIABLE_DECISION.md` §6 stands
  unchanged and is carried into M12.
- Consequence carried forward: **no between-tumour, between-condition or
  differential-abundance claim can be supported from the integrated embedding**, because
  the between-tumour axis is the axis that was deliberately removed.

### 2.3 Dimensions

`dims.use = 1:30` of reduction `pca` — identical to the frozen Phase 1 baseline
(neighbour graph and `umap_preintegration` were both built on PCs 1–30, cumulative
variance 91.1%). Using the same 30 dimensions makes the M12 comparison a comparison of
*integration*, not of dimensionality. The reduction has 50 dims available; availability
of dims 1–30 was asserted before Harmony ran.

### 2.4 Version

`harmony` **1.2.4** · `Seurat` 5.4.0 · `SeuratObject` 5.3.0 · `Matrix` 1.7.4 · R 4.4.3
(`/local/projects-t3/lilab/vmenon/anaconda3/envs/R_env/bin/R`).
Nothing was installed, upgraded or downgraded.

### 2.5 Exact call

```r
set.seed(42)

harmony::RunHarmony(
  object         = obj,
  group.by.vars  = "sample_id",
  reduction.use  = "pca",
  dims.use       = 1:30,
  reduction.save = "postint_harmony",
  verbose        = TRUE
)
```

### 2.6 Explicit vs. default parameters

Machine-readable: `results/phase2/harmony/harmony_parameters.tsv` and
`harmony_parameters.json`.

**Explicitly supplied** (all technically required, or fixed by the M10 decision):

| Parameter | Value | Why explicit |
| --- | --- | --- |
| `group.by.vars` | `"sample_id"` | The M10 decision. |
| `reduction.use` | `"pca"` | Required — names the input embedding. |
| `dims.use` | `1:30` | Required — `NULL` would silently use all 50 PCs. |
| `reduction.save` | `"postint_harmony"` | Required by the Phase 1 → Phase 2 namespace contract. |
| `verbose` | `TRUE` | Diagnostics only; no effect on the result. |
| seed | `42` | Inherited from Phase 1 `config/config.yaml`. |

**Left at package defaults** (no scientific parameter was tuned):

| Parameter | Default value used | Source |
| --- | --- | --- |
| `project.dim` | `TRUE` | `RunHarmony.Seurat` default |
| `theta` | `2` | `NULL` → `rep(2, length(vars_use))` |
| `sigma` | `0.1` | harmony default |
| `lambda` | `1` | harmony default |
| `nclust` | **100** | `NULL` → `min(round(19716/30), 100)` |
| `max_iter` | `10` | harmony default |
| `early_stop` | `TRUE` | harmony default |
| `ncores` | `1` | harmony default |
| `plot_convergence` | `FALSE` | harmony default |
| `alpha` | `0.2` | `harmony_options()` default |
| `tau` | `0` | `harmony_options()` default |
| `block.size` | `0.05` | `harmony_options()` default |
| `max.iter.cluster` | `20` | `harmony_options()` default |
| `epsilon.cluster` | `1e-3` | `harmony_options()` default |
| `epsilon.harmony` | `1e-2` | `harmony_options()` default |

### 2.7 Convergence

Harmony reported **"converged after 9 iterations"** of a maximum 10 — `early_stop`
(default `TRUE`) triggered, so the iteration cap was not the binding constraint.

`results/phase2/harmony/harmony_convergence.tsv`:

| Iteration | Harmony objective | k-means rounds |
| ---: | ---: | ---: |
| 0 | 753.089 | — |
| 1 | 673.748 | 9 |
| 2 | 503.977 | 8 |
| 3 | 439.831 | 6 |
| 4 | 400.310 | 6 |
| 5 | 379.449 | 5 |
| 6 | 367.923 | 5 |
| 7 | 361.789 | 5 |
| 8 | **353.656** | 5 |
| 9 | 354.898 | 5 |

The objective falls monotonically from 753.09 to 353.66 across the first eight
iterations, then ticks up marginally at iteration 9 (+0.35%), which is the relative
change that satisfied `epsilon.harmony = 0.01` and stopped the run. This is ordinary
Harmony convergence behaviour. Per-clustering-step objectives are in
`harmony_kmeans_objective.tsv`.

**Convergence provenance.** `RunHarmony.Seurat` discards Harmony's internal state, so
the history above comes from a second, matrix-level call made with identical input,
identical seed and identical defaults. The two embeddings were compared element-wise:

> **max |primary − diagnostic| = 0.000e+00 → MATCH (bit-identical)**

The convergence history therefore describes the primary run exactly, and the Harmony
result is confirmed deterministic under `set.seed(42)`.

---

## 3. Execution

| Field | Value |
| --- | --- |
| SLURM JobID | **19886486** |
| Job name | `p2_m11_harmony` |
| Account / Partition / Node | `ihc` / `ihc` / `ihc-grid-1-1-1` |
| Requested CPUs | 8 |
| Requested memory | 64G |
| Requested walltime | 04:00:00 |
| **State** | **COMPLETED** |
| **ExitCode** | **0:0** |
| Elapsed | **00:08:24** |
| MaxRSS | **8,788,416 K (8.38 GiB)** |
| TotalCPU | 00:07:58 |
| Memory efficiency | 13.1% of the 64G request |
| Walltime efficiency | 3.5% of the 04:00:00 request |
| Submission | `sbatch scripts/shell/phase2/run_m11_harmony.sh` |

Internal timing (`harmony_parameters.json → timing_seconds`):

| Stage | Seconds |
| --- | ---: |
| Input checksums (md5 + sha256 of 6.04 GB) | 49 |
| `readRDS` of the Phase 1 object | 34.0 |
| **Harmony (primary)** | **10.1** |
| Harmony (diagnostic re-run) | 9.5 |
| `saveRDS` of the Phase 2 object | 256.7 |
| Output checksums | 61 |
| Round-trip re-read + validation | 32 |
| Total R runtime | 485.0 (8.08 min) |

Harmony itself is 10 seconds of the 8-minute job; the wall time is dominated by
serialising and checksumming a 5.6 GB object.

### 3.1 Warnings and errors

- **R warnings collected: 0.** `options(warn = 1)` was set, so any warning would have
  surfaced immediately and been recorded; none occurred.
- **R errors: none.** ExitCode 0:0.
- `logs/phase2/slurm/m11_harmony_19886486.err` contains **only Harmony's own progress
  output** (`Transposing data matrix`, `Initializing state using k-means centroids
  initialization`, the per-iteration progress bars, and `Harmony converged after 9
  iterations`) for both the primary and the diagnostic run. Harmony writes these to
  stderr as `message()` calls; none is an error or a warning.

---

## 4. Output

| Property | Value |
| --- | --- |
| Object | `results/phase2/harmony/phase2_harmony_integrated.rds` |
| Size | 6,048,637,630 bytes (5.63 GiB) |
| md5 | `cf63e84313a91de25fe2e41660f78f06` |
| sha256 | `6085976f788633b68c07c1ddff49e20e6f6f6727579cb7900854febd2ae6e69d` |
| Cells | **19,716** (in = out; no cell was added, dropped or reordered) |
| Reductions | `pca`, `umap_preintegration`, **`postint_harmony`** |
| Graphs | `SCT_nn`, `SCT_snn` |
| Assays | `RNA`, `SCT` (default `SCT`) |
| Metadata columns | 93 (unchanged; Harmony adds none) |
| Harmony reduction | `postint_harmony` — 19,716 cells × **30** dims |
| Harmony key | `postintharmony_` (Seurat strips the underscore from the reduction name when forming a key; this is Seurat's `Key()` sanitisation, not a naming error) |
| Harmony reduction assay | `SCT` |
| Feature loadings projected | `TRUE` (`project.dim = TRUE`, the package default) |

### 4.1 Supporting outputs

```
results/phase2/harmony/harmony_parameters.json                    (full parameter + provenance record)
results/phase2/harmony/harmony_parameters.tsv                     (explicit vs default table)
results/phase2/harmony/harmony_convergence.tsv                    (objective per Harmony iteration)
results/phase2/harmony/harmony_kmeans_objective.tsv               (objective per clustering step)
results/phase2/harmony/harmony_embedding_dimension_summary.tsv    (per-dimension summary, 30 rows)
results/phase2/harmony/harmony_grouping_composition.tsv           (cells per sample_id)
results/phase2/harmony/prov_m11_harmony.json                      (Phase 1-compatible provenance record)
```

### 4.2 Embedding summary (`harmony_embedding_dimension_summary.tsv`)

- 19,716 cells × 30 Harmony dimensions.
- **Non-finite values: 0** across all 30 dimensions.
- Row names of the embedding match the object's cell names, in order.
- Per-dimension spread, alongside the corresponding pre-Harmony PC for reference:

| Dim | Harmony sd | Pre-Harmony PC sd | Harmony min | Harmony max |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 31.343 | 31.446 | −49.32 | 142.77 |
| 2 | 16.702 | 17.958 | −30.93 | 156.21 |
| 3 | 18.205 | 16.119 | −46.96 | 112.31 |
| 4 | 15.956 | 14.676 | −55.35 | 80.93 |
| 5 | 11.507 | 12.559 | −38.56 | 88.42 |
| … | … | … | … | … |
| 30 | 4.179 | 4.710 | −64.82 | 29.74 |

These numbers are recorded as evidence that the embedding is well-formed. **They are
not an integration assessment** — no conclusion about mixing or biological preservation
is drawn from them here.

### 4.3 Validation results

Two independent validation passes were run. Every check returned `TRUE`; any `FALSE`
would have aborted the job with a non-zero exit status.

**(a) Phase 1 preservation proof** — in-memory digests taken before Harmony and
re-taken after:

| Check | Result |
| --- | --- |
| `pca` embedding unchanged (md5 `ba1a61faacc9ae02060b5960b25fc6d7`) | ✅ TRUE |
| `umap_preintegration` embedding unchanged | ✅ TRUE |
| Full metadata frame unchanged (md5 `030b8d676fe22987b8a1e9a75965d294`) | ✅ TRUE |
| Metadata column names unchanged | ✅ TRUE |
| Cell names unchanged | ✅ TRUE |
| Phase 1 graphs `SCT_nn`, `SCT_snn` present | ✅ TRUE |
| Phase 1 reductions present | ✅ TRUE |
| Assays unchanged | ✅ TRUE |
| Default assay unchanged (`SCT`) | ✅ TRUE |

**(b) Round-trip validation** — the saved RDS was re-read from disk after the in-memory
object was freed, because existence of a file is not proof of success:

| Check | Result |
| --- | --- |
| Cell count = 19,716 | ✅ TRUE |
| Cell names unchanged | ✅ TRUE |
| Metadata unchanged (digest matches pre-Harmony) | ✅ TRUE |
| Metadata columns unchanged | ✅ TRUE |
| Phase 1 `pca` present **and** byte-identical | ✅ TRUE |
| Phase 1 `umap_preintegration` present **and** byte-identical | ✅ TRUE |
| Phase 1 graphs present | ✅ TRUE |
| Assays and default assay unchanged | ✅ TRUE |
| `postint_harmony` reduction present | ✅ TRUE |
| Harmony dims = 30 | ✅ TRUE |
| Harmony cells = 19,716 | ✅ TRUE |
| Serialised Harmony embedding identical to the in-memory one | ✅ TRUE |
| Harmony embedding contains no non-finite values | ✅ TRUE |
| Harmony embedding row names match object cell names in order | ✅ TRUE |

**(c) Pre-run guards** — all passed before Harmony was allowed to run: output path is new
and outside every protected Phase 1 location; both input checksums match; 19,716 cells
with no duplicate barcodes; grouping variable has the exact M10-validated structure;
`pca` exists with ≥30 dims and row names matching the cells in order; the reduction name
`postint_harmony` did not already exist.

---

## 5. M10 Finding 2 — CLOSED

M10 flagged that `percent.mt` is identically zero across all 7,615 MPNST_1 cells while
Phase 1 M8 passed `vars.to.regress = "percent.mt"` to SCTransform unconditionally,
regressing a zero-variance covariate for that model. M11 ran the integrity precheck:

| Check | Result |
| --- | --- |
| Non-finite values in `pca`, all 50 dims | **0** |
| Non-finite values in `pca` dims 1:30 (the Harmony input) | **0** |
| Non-finite values in `SCT@scale.data` (5,192 × 19,716) | **0** |
| — of which in MPNST_1 / MPNST_2 / MPNST_3 / MPNST_4 | 0 / 0 / 0 / 0 |
| Non-finite values in the Harmony embedding | **0** |

**The zero-variance regression did not corrupt the residuals or the PCA.** This finding
is closed and requires no correction.

Incidental observation for M14: `SCT@scale.data` holds **5,192** genes, not the 3,000
`variable.features.n` requested in Phase 1 M8. With four SCT models, Seurat retains
residuals for the union of the per-model variable-feature sets. This is expected
behaviour and is noted so that M14 marker plotting does not assume a 3,000-gene
scale.data matrix.

---

## 6. Scientific status

> **Harmony execution has completed, but integration quality has not yet been
> scientifically accepted. Pre/post-Harmony evaluation is reserved for M12.**

What M11 does establish:
- Harmony ran to convergence, deterministically, on the approved input with package
  defaults, and produced a well-formed 30-dimensional embedding for all 19,716 cells.
- The Phase 1 non-integrated baseline is provably intact and remains available for the
  M12 comparison.

What M11 explicitly does **not** establish:
- That technical structure was reduced.
- That biological structure was preserved.
- That over-correction did not occur, in particular for the patient-private malignant
  compartment and for the sample-restricted populations flagged in
  `HARMONY_VARIABLE_DECISION.md` §6.1.
- That default Harmony parameters are adequate for this dataset.

A falling objective function is an optimiser converging on its own criterion. It is
**not** evidence of correct integration, and it is not treated as such here.

---

## 7. Problems, limitations and unresolved questions

1. **The perfect batch/biology confounding is unchanged** and is now baked into the
   embedding. Between-tumour differences have been deliberately suppressed. This
   constrains every downstream interpretation and must be restated in M12–M17.
2. **Over-correction has not been tested.** The three sample-restricted populations
   identified from legacy labels at M10 (B cells ~95% MPNST_3, Malignant SCP-like ~97%
   MPNST_1, cycling tumour cells ~82% MPNST_4) are the specific populations at risk.
   M12 must check them explicitly.
3. **Four SCT models remain** in the object (M10 Finding 1, accepted by the researcher).
   Harmony operates at the embedding level and does not touch expression values, so this
   is untouched by M11. `PrepSCTFindMarkers()` remains mandatory before M14 marker
   discovery, and `SCT@scale.data` spans 5,192 genes rather than 3,000.
4. **Group sizes are unbalanced 3.3×** (MPNST_1 7,615 vs MPNST_2 2,284) and default
   `theta = 2` is applied uniformly. Whether this biased the correction toward the larger
   samples is an open question for the M12 diagnostics.
5. **`nclust` defaulted to 100** — the package cap (`min(round(N/30), 100)` = 100 for
   19,716 cells). Harmony's soft-cluster resolution is therefore at its ceiling for this
   dataset size. Recorded, not changed.
6. **Phase 1 documentation defects D1–D9** (M10 §7) remain uncorrected pending the
   researcher's instruction.
7. **No figures were produced in M11**, so `reports/FIGURE_INDEX.tsv` is unchanged. All
   Harmony figures belong to M12.

None of these blocks M12; items 1, 2 and 4 are precisely what M12 is designed to test.

---

## 8. Recommendation

**The project is technically ready for M12.**

Every M11 completion-gate criterion is met: the approved input object was used and
checksum-verified, the approved grouping variable and dimensions were used, Harmony ran
with package defaults, the Phase 1 state is provably preserved byte-for-byte, the new
`postint_harmony` reduction exists with the expected shape and no non-finite values, the
saved object round-trips, checksums are recorded, and SLURM execution was inspected with
no warnings or errors.

Recommended M12 scope (for authorization, not started):

1. Matched pre/post figures on identical axes and colour scales: `pca` and
   `umap_preintegration` versus a Harmony UMAP, coloured by dataset/sample, by
   `percent.mt`/`percent.ribo`/`nCount_RNA`, and by the Phase 1 `preint_recommended_cluster`
   (provenance only). *(A UMAP on `postint_harmony` is required for M12's paired figures;
   it was deliberately deferred out of M11 per the milestone's diagnostic-only scope.)*
2. Quantitative mixing, computed natively (no `lisi`/`kBET` installed): same-sample kNN
   fraction and neighbourhood entropy in the Harmony space, **directly comparable to the
   frozen Phase 1 values of 0.9613 and 0.0966**, plus an inverse-Simpson (iLISI-equivalent)
   index on `sample_id`.
3. Over-correction diagnostics: per-cell retention of pre-Harmony k nearest neighbours
   after integration; correlation of pre- vs post-Harmony cell–cell distances; silhouette
   width of `sample_id` in both spaces.
4. Targeted checks on the three at-risk sample-restricted populations from §7.2.
5. `reports/phase2/HARMONY_ASSESSMENT.md` separating observed evidence, interpretation,
   over-correction risk, under-correction risk, residual confounding, and a
   recommendation — including, if default Harmony proves inadequate, a *proposed*
   sensitivity analysis for researcher approval rather than any automatic tuning.

---

## 9. STOP

> **STATUS: STOPPED.**
> M11 is complete. No clustering, no neighbours, no UMAP, no markers, no annotation and
> no pre/post scientific assessment were performed. Phase 2 will not proceed to M12
> without explicit researcher authorization.
