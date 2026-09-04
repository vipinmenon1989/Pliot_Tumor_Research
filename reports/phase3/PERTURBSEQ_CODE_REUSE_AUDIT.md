# Perturb-seq Pipeline Code Reuse Audit

**Phase 3 · Milestone M18 · MPNST project**
*Generated: 2026-09-03*

Reference codebase (**READ-ONLY**):
`/local/projects-t3/lilab/vmenon/PertTF-Virtual-Challeng-Weilab/perturbseq-pipeline/src/perturbseq_pipeline`

> **Nothing in the perturb-seq pipeline was modified, and Phase 3 creates no runtime
> dependency on it.** Where logic is adopted it is re-implemented inside
> `scripts/phase3/` with attribution, so `Pilot_MPNST` remains independently reproducible
> (Phase 3 specification §37).

**Guiding principle applied throughout:** scientific suitability outranks software reuse. A
component was rejected wherever the perturb-seq design assumption (a designed perturbation
screen with non-targeting controls, tens of thousands of cells per condition, and a
guide-level replication structure) does not hold for a four-patient tumour atlas.

Decisions use: `REUSE` · `ADAPT` · `REFERENCE ONLY` · `NOT APPLICABLE`.

---

| source_file | function/class | original_purpose | potential_phase3_use | scientific_compatibility | technical_compatibility | reuse_decision | modifications_required | notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `lochness.py` | module docstring + `lochness_score`, `_compute_self_lochness` | lochNESS per cell per perturbation: `local_fraction(g)/overall_fraction(g) − 1` | The core LochNESS formula for receiver-state analysis | **High** — the formula is design-agnostic and matches Huang et al. exactly | Medium — Python/AnnData; MPNST is Seurat/R | **ADAPT** | Re-implemented in R. Must add **same-sample exclusion** (absent here, present in official MMCA) and switch *k* to `round(0.5·√N)` | The single most valuable component. See `LOCHNESS_IMPLEMENTATION_AUDIT.md` §4–5 for the two corrections required |
| `lochness.py` | `_adjacency` | binarise a kNN graph and return the true per-row neighbour count | Same idea when converting an R kNN result to a fraction | High | Low — trivial logic, not worth a language port | **REFERENCE ONLY** | — | The "divide by *actual* neighbour count, not requested k" observation was adopted as a documented choice; MMCA divides by requested `kadj` and Phase 3 follows MMCA for fidelity |
| `lochness.py` | `_count_same_label_neighbors_numba` | Numba parallel kernel for same-label neighbour counting at million-cell scale | Speed | High | **Unnecessary** — largest MPNST receiver lineage is 5,064 cells; exact kNN at k≈36 is instant | **NOT APPLICABLE** | — | Solves a scale problem MPNST does not have; adopting it would add a Numba dependency for no benefit |
| `lochness.py` | `_build_neighbor_graph` | dedicated large-k (~300) scanpy neighbour graph, prefers `X_pca_harmony` | Neighbourhood construction | **Medium** — the "use the batch-corrected embedding" instinct is right, but k=300 is not the published rule and there is no L2 normalisation | Medium | **REFERENCE ONLY** | — | Phase 3 follows MMCA: `Seurat::L2Dim` on `postint_harmony`, exact `FNN::get.knnx`, `k = round(0.5·√N)` |
| `lochness.py` | `compute_lochness` (orchestration), `LochnessResults`, `attach_scores` | end-to-end orchestration, min-cell target filtering, result attachment | Output-shape and min-cell-filtering patterns | High | Low — tied to `Config`/AnnData | **REFERENCE ONLY** | — | The "exclude targets with too few cells and record *why*" pattern was adopted conceptually into the M19 minimum-cell policy |
| `compute.py` | `detect_slurm_cpus`, `detect_available_cpus`, `resolve_worker_count` | resolve worker count from SLURM allocation rather than machine size | Sizing `future`/`BiocParallel` workers inside SLURM jobs | High | Medium — the *idea* ports; `SLURM_CPUS_PER_TASK` is read directly in R | **ADAPT** | One-line R equivalent | Prevents oversubscription; Phase 3 shell scripts pass `SLURM_CPUS_PER_TASK` into R explicitly |
| `compute.py` | `limit_blas_threads` | context manager pinning BLAS threads during parallel work | Avoiding BLAS × worker thread explosion in CCC jobs | High | Medium | **ADAPT** | `RhpcBLASctl::blas_set_num_threads` equivalent | Same failure mode exists in R; harmony 1.2.4 already does this internally |
| `compute.py` | `derive_seed` | derive a reproducible per-item seed from a base seed | Per-sample and per-permutation seeding | **High** | High — trivially portable | **ADAPT** | R equivalent | Adopted for the LochNESS permutation loop so every shuffle is independently reproducible |
| `compute.py` | `run_parallel`, `_parallel_worker_shim`, `resolve_stage_backend`, `ComputeDecision` | generic parallel map with backend selection (CPU/GPU) | Parallelising per-sample CCC | Medium | Low — Python multiprocessing; no GPU in the Phase 3 envelope | **REFERENCE ONLY** | — | Phase 3 parallelises at the **SLURM job** level instead (independent LIANA / CellChat / CellPhoneDB jobs), which is simpler and auditable |
| `distance.py` | `distance_test_permutation` | permutation test for perturbation distance, preserving group sizes | Structure of the LochNESS permutation null | **Medium** — the permutation *machinery* is sound, but the perturb-seq null shuffles guide labels against NTCs, which has no MPNST analogue | High | **ADAPT** | Null redesigned to be **sample-aware** (see `LOCHNESS_MPNST_DESIGN.md`); only the loop/aggregation scaffolding is borrowed | MMCA's own null shuffles the label globally; with four confounded samples that is inappropriate here, so the deviation is documented |
| `distance.py` | `compute_energy_distance`, `compute_mmd`, `energy_distance_from_cdist` | energy distance / MMD between perturbation and control cell clouds | Quantifying receiver-state shift between signalling contexts | **Low** — these compare a perturbation against a *designed control*. MPNST has no control population, and any "context" contrast is confounded with patient | High | **NOT APPLICABLE** | — | Rejected on scientific, not technical, grounds. Using them would manufacture a control/treatment framing the data cannot support |
| `distance.py` | `_sample_cell_indices` | deterministic stratified cell subsampling with a fixed seed | Balanced subsampling for the R-vs-Python LochNESS cross-check | High | High | **ADAPT** | R equivalent | Used in M24 to build the deterministic comparison subset |
| `distance.py` | `compute_pcoa_coordinates`, `compute_distance_space` | PCoA over a perturbation-distance matrix | — | Low | High | **NOT APPLICABLE** | — | No perturbation space exists in MPNST |
| `enrichment.py` | `omnibus_test`, `_cmh_test` | Cochran–Mantel–Haenszel test stratified by lane/batch | Testing whether a receiver state is enriched across samples **while stratifying by sample** | **Medium–High** — CMH stratification is exactly the right instinct for sample-aware inference | High — `stats::mantelhaen.test` is in base R | **ADAPT** | Use base R `mantelhaen.test`; strata = `sample_id` | Adopted only where a stratified count comparison is genuinely warranted; **not** used to manufacture cell-level p-values (Phase 3 §51) |
| `enrichment.py` | `_odds_ratio`, `_reference_mask` | odds ratio of a label in a cluster vs reference | Receiver-state composition summaries | Medium | High | **ADAPT** | R equivalent | Effect size preferred over p-values, per Phase 3 §51 |
| `enrichment.py` | `test_cluster_enrichment`, `_test_cluster_enrichment_standard/_large`, `_build_large_count_tables` | full cluster-enrichment engine with large-data paths | — | Low — built around targeting-vs-NTC contrasts | Low | **REFERENCE ONLY** | — | Conceptual only; the MPNST question is not "is this guide enriched in this cluster" |
| `enrichment.py` | `_guide_concordance` | agreement between independent guides targeting the same gene | **Conceptual template for multi-method CCC concordance** — independent evidence streams for the same underlying claim | **High as a pattern** | Low as code | **REFERENCE ONLY** | — | The M21 concordance design (per-method support flags, then a class, never averaged raw scores) follows this pattern |
| `gene_sets.py` | `benjamini_hochberg` | BH FDR correction | Multiplicity control where tests are performed | High | High — `stats::p.adjust(method="BH")` exists | **NOT APPLICABLE** | — | Base R already provides it; porting would be pointless |
| `gene_sets.py` | `run_ora_enrichment`, `run_program_enrichment`, `parse_gmt` | over-representation analysis of gene programs | Receiver-response pathway characterisation | **Low** — Phase 3 §69 prohibits "large unrelated pathway screens"; receiver response is handled by NicheNet's ligand→target framework | High | **NOT APPLICABLE** | — | Deliberately declined to avoid scope creep |
| `gene_sets.py` | `clean_term_name`, `format_display_label`, `normalize_species_name`, `adapt_gene_set_to_species` | label hygiene and human/mouse gene-symbol mapping | Gene-symbol normalisation when reconciling LR resources across methods | **High** — symbol/complex normalisation is a real M21 problem | Medium | **ADAPT** | R equivalent, human-only | Adopted as the basis for the canonical `sender|receiver|ligand|receptor` key |
| `ps_score.py` | `_prepare_for_pertps`, `_expression_layer`, `_expression_cut`, `_classify_quadrants` | perturbation-score computation via the external `pertps` package | — | Low — no perturbation scores in MPNST | Medium | **NOT APPLICABLE** | — | |
| `ps_score.py` | `_pct_control_expressing_standard/_large` | fraction of control cells expressing a gene, with a large-data path | **Expression-fraction thresholds for CCC inputs** — every LR method needs "is this ligand expressed in enough of the sender cells" | **High** | Medium | **ADAPT** | R equivalent computed on the SCT `data` layer per sample × population | Adopted for the M19 input-QC tables |
| `modules.py` | `select_genes`, `_stratified_marker_sample`, `_build_effect_matrix_*` | gene-module discovery over a perturbation × gene effect matrix | — | Low | Low | **NOT APPLICABLE** | — | NMF-style module discovery on a perturbation effect matrix has no MPNST analogue |
| `modules.py` | `_as_sparse_float64`, `_sparse_expm1`, `_sparse_square`, `_dense_layer` | safe sparse-matrix helpers that avoid accidental densification | Guarding against dense `cells × genes` blow-ups in CCC input prep | **High as a discipline** | Low as code — R has `Matrix` | **REFERENCE ONLY** | — | The discipline (never densify; `expm1` on `@x` in place) is adopted; Phase 3 §52 forbids dense `cells × genes` / `cells × cells` matrices |
| `io.py` | `load_data`, `_load_h5ad`, `_load_mtx`, `apply_layer_choices` | data loading and layer selection | — | Medium | Low — MPNST input is a Seurat RDS | **NOT APPLICABLE** | — | |
| `io.py` | `attach_sample_metadata`, `_lanes_from_obs` | attach and validate per-sample metadata; derive lane IDs | Validating `sample_id` integrity before per-sample CCC | High | Low | **REFERENCE ONLY** | — | Phase 2 already validated this metadata exhaustively (M17) |
| `io.py` | `_check_matching_vars` | assert identical feature spaces before concatenating objects | Asserting gene-space consistency across per-sample CCC inputs | High | High | **ADAPT** | R equivalent | Adopted as an M19 assertion |
| `qc.py` | (module) | perturb-seq QC: guide calling rates, MOI, doublets | — | Low | Low | **NOT APPLICABLE** | — | Phase 1 already handled QC |
| `guides.py` | (module) | guide-to-target assignment, NTC classification | — | **None** | — | **NOT APPLICABLE** | — | No guides in MPNST |
| `plots.py` | UMAP/heatmap/volcano helpers (178 KB) | plotting | Figure conventions | Medium | Low — matplotlib; Phase 3 uses ggplot2 | **REFERENCE ONLY** | — | Only the convention of emitting PDF + PNG with deterministic filenames and indexing every figure was carried across, which Phase 2 already established |
| `report.py` | markdown report assembly | reporting | Report structure | Medium | Low | **REFERENCE ONLY** | — | Phase 2's milestone-report convention is already established and is retained |
| `config.py` | `Config` dataclass tree (58 KB) | typed hierarchical configuration | Phase 3 configuration | Medium | Low — Python dataclasses | **REFERENCE ONLY** | — | Phase 3 follows the established repository convention: TSV/YAML config plus explicit CLI arguments |
| `cli.py` | command-line entry points | orchestration | — | Medium | Low | **REFERENCE ONLY** | — | Phase 3 uses one R script per milestone plus an `sbatch` wrapper, matching Phase 2 |
| `meta.py` | `build_perturbation_meta` | build a perturbation metadata table | — | Low | Low | **NOT APPLICABLE** | — | |
| `cluster.py` | (module) | Leiden clustering and neighbour graphs | — | Low | Low | **NOT APPLICABLE** | — | Phase 2 M13 already produced the clustering |
| `data_access.py` | accessor helpers | — | — | Low | Low | **NOT APPLICABLE** | — | |

---

## Summary of decisions

| Decision | Count | Components |
| --- | ---: | --- |
| **ADAPT** | 9 | LochNESS core formula · SLURM-aware worker count · BLAS thread pinning · `derive_seed` · permutation scaffolding · deterministic stratified subsampling · CMH stratified test · odds ratio · gene-symbol normalisation · expression-fraction thresholds |
| **REFERENCE ONLY** | 11 | adjacency helper · neighbour-graph builder · LochNESS orchestration · `run_parallel` · cluster-enrichment engine · guide-concordance pattern · sparse-matrix discipline · sample-metadata validation · plotting · reporting · config |
| **NOT APPLICABLE** | 13 | Numba kernel · energy distance / MMD · PCoA / distance space · BH (base R has it) · ORA / program enrichment · pertps scoring · module discovery · IO · QC · guides · meta · clustering · data access |

## Two components rejected on scientific rather than technical grounds

1. **Energy distance / MMD (`distance.py`).** Technically portable and well written, but they
   measure the displacement of a perturbed cell cloud from a *designed control*. MPNST has no
   control population, and any "signalling context" contrast is perfectly confounded with
   patient identity. Using them would impose a treatment/control framing the data cannot
   support.
2. **ORA / gene-program enrichment (`gene_sets.py`).** Would be easy to run, but Phase 3 §69
   prohibits large unrelated pathway screens, and receiver response is properly addressed by
   NicheNet's ligand→target regulatory framework rather than by generic set enrichment.

## Attribution

Where logic was adapted, the Phase 3 source file carries a header comment naming the
originating module and function, e.g.:

```r
# LochNESS core formula adapted from:
#   perturbseq-pipeline/src/perturbseq_pipeline/lochness.py (read-only reference)
# corrected against the official implementation:
#   shendurelab/MMCA @ af629c49 — Section_5_step_1_run_lochness_calculation.R
# Corrections applied: same-sample exclusion added; k = round(0.5*sqrt(N)); L2-normalised PCA.
```
