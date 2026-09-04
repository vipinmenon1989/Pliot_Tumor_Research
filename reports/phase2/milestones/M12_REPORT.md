# Milestone M12 — Pre/Post Harmony Evaluation

**Phase 2 · MPNST single-cell integration**
*Generated: 2026-09-02*
*Status: **COMPLETE — STOPPED, awaiting authorization for M13***

> **Harmony recommendation: `ACCEPT DEFAULT HARMONY WITH CAVEATS`**
> Full scientific argument: [`reports/phase2/HARMONY_ASSESSMENT.md`](../HARMONY_ASSESSMENT.md)

---

## 1. Scripts created and modified

### Created

| Path | Purpose |
| --- | --- |
| `scripts/R/phase2/evaluate_harmony.R` | M12 evaluation: matched pre/post UMAPs, exact-kNN mixing metrics, silhouette widths, within-sample structure preservation, canonical-program biological references, sample-restricted population check, 12 figures, 12 machine-readable outputs. Includes `--validation-mode` for synthetic smoke testing. |
| `scripts/shell/phase2/run_m12_evaluate.sh` | SLURM launcher (8 CPUs, 96G, 04:00:00, `ihc` / `ihc-grid-1-1-1`). |
| `reports/phase2/HARMONY_ASSESSMENT.md` | The M12 scientific assessment (sections A–J). |
| `reports/phase2/milestones/M12_REPORT.md` | This file. |

### Modified

| Path | Change |
| --- | --- |
| `reports/FIGURE_INDEX.tsv` | Added the columns `phase`, `milestone` and `slurm_job_id`; back-filled all 260 Phase 1 rows as `phase1` / `NA_backfilled_phase1`; appended 24 M12 rows. 284 rows total. `snakemake_rule` retained for backwards compatibility; Phase 2 rows carry `NA_no_snakemake_in_phase2`. |
| `PROGRESS.md`, `CHANGELOG.md` | M12 status and record. |

### Deliberately not modified

`processed_mpnst.rds`, `results/datasets/**`, `results/combined/**` (Phase 1 handoff mtime
still 2026-07-19 23:25), `results/phase1_manifest.json`,
`results/phase2/harmony/phase2_harmony_integrated.rds` (the M11 object was opened
read-only; its checksums were re-verified and its embedding digest re-checked after all
computation), `config/config.yaml`, `workflow/Snakefile`, every Phase 1 report.

---

## 2. Commands

```bash
# validated on a synthetic structural mirror first (login node, 800 cells, no real data)
Rscript scripts/R/phase2/evaluate_harmony.R --validation-mode --input <synthetic>.rds ...

# production
sbatch scripts/shell/phase2/run_m12_evaluate.sh
#   -> Rscript scripts/R/phase2/evaluate_harmony.R \
#        --input results/phase2/harmony/phase2_harmony_integrated.rds \
#        --out-dir results/phase2/harmony/evaluation \
#        --fig-dir reports/phase2/figures/m12 \
#        --group-by sample_id --dims 30 --k-primary 15 --k-robust 50 \
#        --random-seed 42 --expected-cells 19716
```

No Snakemake. Harmony was neither re-run nor re-tuned.

---

## 3. SLURM execution

| Field | Value |
| --- | --- |
| JobID | **19886628** |
| Job name | `p2_m12_evaluate` |
| Account / Partition / Node | `ihc` / `ihc` / `ihc-grid-1-1-1` |
| Requested CPUs | 8 |
| Requested memory | 96G |
| Requested walltime | 04:00:00 |
| **State** | **COMPLETED** |
| **ExitCode** | **0:0** |
| Elapsed | **00:04:09** |
| MaxRSS | **8,338 M (8.14 GiB)** |
| TotalCPU | 00:03:35 |
| Memory efficiency | 8.5% of request |
| Walltime efficiency | 1.7% of request |

Requested resources were sized from M11 accounting (8.38 GiB) plus the two 19,716 × 19,716
distance matrices needed for the silhouette widths (~1.55 GB each, computed sequentially).
The request was conservative; actual usage came in under the M11 peak because the Seurat
object is released from memory before the metric stage.

Internal timing: input checksums 48 s · `readRDS` 33.6 s · Harmony UMAP 16.0 s · kNN
(k=15 and k=50, both spaces) 16 s · within-sample structure 3 s · distance matrices and
four silhouettes 55 s · 24 figures 22 s · provenance 17 s. Total R runtime 216.4 s.

### Warnings and errors

- **R warnings collected: 0** (`options(warn = 1)` in force). **R errors: none.**
- `logs/phase2/slurm/m12_evaluate_19886628.err` contains a single informational Seurat
  message (`RunUMAP` now defaults to R-native uwot with the cosine metric). Not a warning
  about the data; it documents the UMAP backend, which is identical to the one Phase 1
  used for `umap_preintegration`, preserving the fairness of the comparison.
- Post-job inspection confirmed: all 24 figures present and non-trivial in size, all 12
  tables present, no `NA`/non-finite values in any metric column, all four sample levels
  represented in every breakdown, no metric outside its mathematically valid range.

---

## 4. Input and output paths

**Input:** `results/phase2/harmony/phase2_harmony_integrated.rds`
md5 `cf63e84313a91de25fe2e41660f78f06` · sha256 `6085976f788633b68c07c1ddff49e20e6f6f6727579cb7900854febd2ae6e69d`
— both re-verified at runtime against the M11 record; the job would have aborted on a mismatch.

**Outputs:** `results/phase2/harmony/evaluation/` (12 files) and
`reports/phase2/figures/m12/` (24 files). No new Seurat object was created; the Harmony
UMAP is persisted as a 2-column embedding only (`umap_harmony_m12_embedding.{tsv,rds}`),
avoiding a third multi-gigabyte object.

---

## 5. Analysis configuration

| Item | Value |
| --- | --- |
| Cells evaluated | **19,716 (all)** |
| Subsampling | **None anywhere in M12** |
| PRE space | `pca` dims 1:30 (frozen Phase 1 M8) |
| POST space | `postint_harmony` dims 1:30 (M11) |
| Pre-Harmony UMAP | `umap_preintegration` — Phase 1, **reused unmodified** |
| Post-Harmony UMAP | `umap_harmony_m12` — identical `RunUMAP` call, only the input reduction changed (dims 1:30, seed 42, `n.neighbors` 30, `min.dist` 0.3, metric cosine — all Seurat defaults) |
| Neighbourhood sizes | k = 15 (primary, matches Phase 1) and k = 50 (robustness) |
| kNN | `RANN::nn2`, exact |
| Silhouette | `cluster::silhouette` on full euclidean pairwise distances |
| Seed | 42 |
| Harmony re-tuned | **No** |

**Fairness.** Identical cells, dimensions, metadata, neighbour parameters, UMAP
parameters and seed on both sides. Nothing was optimised for appearance. Digests of the
Phase 1 PCA, the Phase 1 UMAP and the M11 Harmony embedding were re-checked after all
computation and were unchanged.

**Unavailable packages.** `lisi`, `kBET`, `clustree` are not installed and were **not**
installed. Substitutes: inverse Simpson computed natively on exact kNN (the same estimator
`lisi` uses) in place of iLISI; a dominance ratio against global composition in place of
kBET. Documented in `harmony_evaluation_parameters.tsv`.

---

## 6. Scientific findings

### 6.1 PRE-Harmony metrics reproduce Phase 1 exactly

Per-sample same-sample neighbour fractions match Phase 1 M8 to 15 significant figures
(e.g. MPNST_1 `0.988321295688334` in both). The two headline numbers differ only in
presentation: Phase 1's `0.9613` is the *unweighted* mean of four per-sample means
(cell-weighted = `0.9728`, the M12 value), and Phase 1's entropy `0.0966` is in **log2**
(cell-weighted 0.0707 bits = 0.0490 nats, the M12 value). This validates the M12 pipeline
against the frozen baseline.

### 6.2 Technical mixing — improved substantially

Because the samples are unequal in size, "fully mixed" means neighbourhoods matching
global composition: same-sample fraction 0.3065, entropy 1.2683 nats, inverse Simpson
3.2627 — **not** 0, log K and 4.

| Metric (k = 15, all cells) | PRE | POST | Fully-mixed ref. | Gap closed |
| --- | ---: | ---: | ---: | ---: |
| Same-sample neighbour fraction | 0.9728 | **0.7784** | 0.3065 | 29.2% |
| Neighbourhood entropy (nats) | 0.0490 | **0.3565** | 1.2683 | 25.2% |
| Inverse Simpson (iLISI-equiv.) | 1.0549 | **1.4408** | 3.2627 | 17.5% |
| Dominance ratio | 3.30 | **2.53** | 1.00 | — |

At k = 50 the picture is the same and slightly stronger (same-sample 0.9443 → 0.7053,
37.5% of the gap closed). Every sample improved; MPNST_2 most (0.910 → 0.504), MPNST_1
least (0.988 → 0.849).

### 6.3 The decisive result — mixing is compartment-specific

| Compartment (canonical program) | n | Same-sample fraction PRE → POST | Δ |
| --- | ---: | --- | ---: |
| Panleukocyte | 1,396 | 0.933 → **0.527** | −0.406 |
| Myeloid | 2,441 | 0.958 → **0.611** | −0.347 |
| T/NK | 1,093 | 0.953 → **0.689** | −0.264 |
| Endothelial | 966 | 0.962 → **0.709** | −0.253 |
| Mural | 932 | 0.951 → **0.724** | −0.227 |
| Fibroblast | 4,075 | 0.990 → 0.839 | −0.151 |
| **Schwann/neural-crest** | 2,087 | **0.990 → 0.900** | **−0.090** |
| B/plasma | 836 | 0.986 → 0.918 | −0.067 |

Harmony mixed the shared immune and vascular compartments hard while leaving the
Schwann/neural-crest compartment — the presumptive malignant lineage, 86.5% confined to
MPNST_1 — nearly as patient-private as it began. **This is the desired behaviour**, and it
explains why the *global* mixing figure looks modest: the cell-weighted average is
dominated by compartments in which patient-private structure is expected to survive.

### 6.4 Biological preservation — retained

Two non-legacy references were used. **Legacy `orig.anno` was deliberately excluded** as
legacy integration-derived annotation (M12 §12, and the M10 default).

*Reference 1 — the 54 independently derived Phase 1 per-sample clusters:*

| Measure | Value |
| --- | ---: |
| Within-sample kNN retention (k = 15) | **0.8206** |
| Phase 1 cluster coherence PRE → POST | 0.8644 → **0.8427** (−2.5% relative) |
| Global kNN retention | 0.6598 |

Global retention falling to 0.66 *is* the correction. The discriminating figure is that
82% of each cell's own-sample neighbours survived and agreement with integration-free
Phase 1 structure fell only 2.5%.

*Reference 2 — ten canonical broad-lineage program scores* (documented gene sets, used as
a sanity check, **not** an annotation; 26.3% of cells left `Unassigned`):

Overall program silhouette **rose** 0.0627 → 0.0729. Six of ten programs became more
cohesive, the largest gains being T/NK (+0.209), endothelial (+0.218),
Schwann/neural-crest (+0.170) and myeloid (+0.104) — cells of the same broad lineage from
different patients were brought together.

### 6.5 Over-correction — present but bounded

- **B/plasma cohesion fell 44%** (silhouette 0.437 → 0.243, n = 836, 93.5% MPNST_3). Its
  same-sample fraction barely moved, so the cells were not dispersed across samples; they
  became less separable from adjacent lymphoid/myeloid territory. Classic failure mode for
  a batch-private population with no counterpart to align to.
- **Fibroblast cohesion fell 44%** (0.174 → 0.097, n = 4,075, 70.6% MPNST_4), with the
  lowest within-sample retention of any program (0.795).
- **MPNST_4** shows the largest within-sample damage (retention 0.781, coherence −3.6%),
  consistent with the Phase 1 flag on its `percent.mt`-correlated cluster C09.
- Mast cells (n = 40) drop most in relative terms but are too few to weigh.

Against widespread over-correction: overall program cohesion rose, within-sample retention
is 0.82, the Phase 1 clusters remain discrete islands, and the Schwann/neural-crest
compartment was neither dispersed nor collapsed — it became *more* cohesive.

### 6.6 Under-correction — real in absolute terms, largely appropriate

Only 29.2% of the achievable mixing gap was closed and every sample's dominance ratio
remains above 2.2, highest for the two small samples (MPNST_2 4.35, MPNST_3 4.95) —
uniform default `theta = 2` under a 3.3× size imbalance corrects minority samples
proportionally less. But §6.3 shows the residual sits where it belongs. The one genuine
concern is the weaker relative correction of MPNST_2 and MPNST_3.

### 6.7 A metric that failed

Global silhouette of `sample_id` was **−0.0104 pre-Harmony** despite 97% same-sample
nearest neighbours, because each sample spans the full range of malignant, immune and
stromal states, so within-sample and between-sample mean distances nearly cancel. Batch
structure here is **local**; only the kNN metrics detect it. Reported as a methodological
caveat, contributing nothing to the conclusion. The same machinery *is* informative for
the compact program groups.

### 6.8 Confounding

`sample_id` is simultaneously dataset, patient and the only batch proxy; no clinical or
technical covariate exists. Integration quality therefore cannot be cleanly separated from
biological preservation, and no metric can prove the removed variance was technical. The
compartment-resolved result in §6.3 is the strongest available evidence precisely because
it is *differential* — but it is a consistency argument, not proof. **The integrated
embedding cannot support any between-tumour, between-patient, between-condition or
differential-abundance claim.**

---

## 7. Harmony recommendation

# `ACCEPT DEFAULT HARMONY WITH CAVEATS`

Both required criteria are met: technical mixing improved substantially and in the right
compartments, and biological structure was preserved. Two caveats must travel into M13–M17:

- **C1.** Any **B/plasma** cluster arising downstream may be distorted (44% cohesion loss);
  cross-check it against the preserved non-integrated baseline before annotating.
- **C2.** Any **fibroblast/stromal** cluster, and MPNST_4-derived structure generally,
  warrants the same cross-check.

**A sensitivity analysis is not recommended and was not run.** Raising `theta` would push
hardest on the Schwann/neural-crest and B/plasma populations that are already most
fragile; lowering it would undo the immune/vascular alignment that justifies integrating
at all. If the researcher wants one anyway, the two experiments worth running — neither
executed, both needing explicit approval — are a `theta` ∈ {1, 2, 4} sweep re-scored
through this same metric suite, and a size-imbalance-aware `theta`/`lambda` weighting
targeting the under-correction of MPNST_2 and MPNST_3.

---

## 8. Figures (24 files: 12 figures × PDF + PNG) in `reports/phase2/figures/m12/`

| File | Content |
| --- | --- |
| `01_pre_post_umap_by_sample_id` | **Primary paired figure.** `sample_id` is simultaneously the dataset, the patient and the Harmony grouping variable, so M12 spec figures A, B, C and E are one comparison — stated on the figure rather than duplicated four times. |
| `02_pre_post_umap_faceted_by_sample_id` | Per-sample occupancy / density (spec figure F). |
| `03_pre_post_umap_by_qc_covariates` | `nCount_RNA`, `nFeature_RNA`, `percent.mt`, `percent.ribo`. |
| `04_pre_post_umap_by_phase1_independent_clusters` | The 54 integration-free Phase 1 clusters; legend suppressed for readability. |
| `05_pre_post_neighborhood_entropy` | With the global-composition reference line at 1.2683. |
| `06_pre_post_same_sample_neighbor_fraction` | |
| `07_pre_post_inverse_simpson` | iLISI-equivalent. |
| `08_pre_post_technical_silhouette` | Technical and biological panels side by side. |
| `09_knn_retention_global_vs_within_sample` | The key over-correction diagnostic. |
| `10_biological_preservation_summary` | |
| `11_pre_post_umap_by_canonical_program` | Sanity-check programs, not annotation. |
| `12_local_sample_dominance_vs_global` | Observed vs expected from global abundance. |

**No figure for biological condition** (spec figure D) — that field does not exist in the
data and none was fabricated.

---

## 9. Tables in `results/phase2/harmony/evaluation/`

`pre_post_mixing_summary.tsv` (every metric × k × space with mean/sd/quartiles, desirable
direction and a written definition) · `neighborhood_mixing_metrics_by_sample.tsv` ·
`technical_silhouette_summary.tsv` · `biological_preservation_summary.tsv` ·
`sample_restricted_population_check.tsv` · `biological_program_gene_sets.tsv` ·
`program_by_sample_composition.tsv` · `harmony_evaluation_parameters.tsv` ·
`m12_headline_metrics.json` · `prov_m12_evaluation.json` ·
`umap_harmony_m12_embedding.{tsv,rds}` · `figure_index_m12.tsv`

---

## 10. Unresolved questions

1. **Residual `sample_id` structure cannot be apportioned** between uncorrected batch
   effect and genuine between-patient tumour biology. Permanent, by study design.
2. **C1/C2** — whether the B/plasma and fibroblast cohesion losses materially affect the
   M13 clustering and M15 annotation. Testable only once clusters exist.
3. **Small-sample under-correction** — MPNST_2 and MPNST_3 retain the highest dominance
   ratios. Does the researcher want the size-imbalance sensitivity analysis?
4. **Four SCT models** (M10 Finding 1) remain; `PrepSCTFindMarkers()` stays mandatory for
   M14, and `SCT@scale.data` spans 5,192 genes rather than 3,000.
5. **Phase 1 documentation defects D1–D9** (M10 §7) remain uncorrected pending instruction.
6. **Only Harmony was evaluated.** CCA/RPCA/MNN, proposed in Phase 1
   `INTEGRATION_PREPARATION.md`, were not benchmarked and are out of Phase 2 scope as
   currently authorised.

---

## 11. Recommended M13 configuration (proposed, not started)

| Item | Proposal | Rationale |
| --- | --- | --- |
| Reduction | `postint_harmony` | The M12-accepted embedding. |
| Dimensions | `1:30` | Consistent with M11/M12 and the frozen baseline. |
| UMAP | Reuse the exact `umap_harmony_m12` parameters and seed | Already computed and validated in M12; recomputation is deterministic. Store as `postint_umap_harmony` in the M13 object per the Phase 1 namespace contract. |
| Neighbours | `FindNeighbors(reduction = "postint_harmony", dims = 1:30)`, graphs `postint_harmony_nn` / `postint_harmony_snn` | Phase 1 graphs `SCT_nn`/`SCT_snn` preserved untouched. |
| `k.param` | 20 | Matches Phase 1 `config/config.yaml: clustering.k_param`. |
| Algorithm | Louvain (algorithm 1) | Matches Phase 1. |
| Resolutions | 0.1 … 1.0 in steps of 0.1, **all preserved** | Per the Phase 2 specification. |
| Metadata names | `postint_harmony_clusters_res_0.1` … `_res_1.0` | Namespace contract. |
| Per-resolution outputs | UMAP, cluster sizes, sample/dataset composition, transition/stability diagnostics, tiny-cluster and sample-dominated-cluster warnings, ranked comparison table | Per the Phase 2 specification. |
| **M12-specific additions** | Flag every cluster that is (a) >60% one sample, (b) enriched for the **B/plasma** program (C1), or (c) enriched for the **fibroblast** program (C2), and report each against the preserved non-integrated baseline | Carries the M12 caveats forward rather than leaving them in a report. |
| Resources | 8 CPUs, 96G, 06:00:00 | M12 used 8.14 GiB; a 10-resolution sweep plus 10 UMAPs and a new object save is the main additional cost. |

M13 must not select the final resolution autonomously — it recommends a primary and an
alternative and stops for researcher approval.

---

## 12. STOP

> **STATUS: STOPPED.**
> No clustering sweep, no resolution selection, no `FindMarkers`/`FindAllMarkers`, no
> annotation, no differential expression or abundance, no pathway/trajectory/CNV analysis
> was performed. Harmony was not re-run or re-tuned. M13 requires explicit researcher
> authorization.
