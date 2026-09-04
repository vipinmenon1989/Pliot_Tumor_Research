# Phase 5 Feasibility — Malignant Transcriptional Programs

*Milestone M36 · 2026-09-03 · SLURM 19899335 (COMPLETED 00:02:05, 4 CPUs, ReqMem 96 G, MaxRSS 29.28 GiB)*

> Phases 1–4 are complete and frozen. The Phase 4 object is opened **read-only**; every Phase 5
> artefact is written under `results/phase5/`. Its md5 was re-verified after M36 and is unchanged.

## 1. Input integrity

`results/phase4/phase4_final_object.rds` — 6,057,010,311 bytes — was verified by **both** checksums
before it was loaded:

```text
md5     e85ba8486e456917e2483f2773bdbaf3   (matches)
sha256  a658447744e5ee8621fc5a3127492ceb11bce37bac4d6f3f0f7255bdcb647ae4   (matches)
```

19,716 cells · RNA (31,764 features) + SCT · 4 reductions · 183 metadata columns.
**13 frozen Phase 4 counts re-asserted and all passed**: refined 6,434 / 9,078 / 3,766 / 438;
Fibroblast 5,064 → 4,036 / 908 / 120; clones 7 / 4 / 3 / 8.

## 2. The historical annotation field was VERIFIED, not assumed

M39 compares historical fibroblasts that Phase 4 called malignant against historical fibroblasts
that stayed non-malignant. The grouping variable must therefore be the pre-Phase-4 annotation, and
`annotation_ccc_refined` is **barred** from that role because it already encodes the Phase 4
conclusion — using it would make the comparison circular.

Every `annotation*` column was tested against the frozen fibroblast split:

| column | levels | n Fibroblast | → Malignant | → Non-malignant | → Ambiguous | reproduces frozen split | eligible |
| --- | ---: | ---: | ---: | ---: | ---: | :--: | :--: |
| `annotation_ccc` | 17 | 5,064 | 4,036 | 908 | 120 | **yes** | yes |
| `annotation_ccc_compartment` | 6 | 0 | 0 | 0 | 0 | no | yes |
| `annotation_ccc_is_tumor` | 2 | 0 | 0 | 0 | 0 | no | yes |
| `annotation_ccc_ccc_ready` | 2 | 0 | 0 | 0 | 0 | no | yes |
| **`annotation_ccc_phase3`** | 17 | **5,064** | **4,036** | **908** | **120** | **yes** | **yes** |
| `annotation_ccc_refined` | 17 | 908 | 0 | 908 | 0 | no | **NO — encodes the conclusion** |

**Selected: `annotation_ccc_phase3`** — the expected field, and it qualified on the evidence rather
than by assumption. `annotation_ccc` is identical to it and would have served equally.
Audit: `results/phase5/validation/M36_HISTORICAL_ANNOTATION_FIELD_AUDIT.tsv`.

## 3. Required fields

All 16 required Phase 4 fields are present and populated. Two points a reader should know:

* `tumor_clone_phase4` is non-NA for **6,256** of the 6,434 malignant cells — Phase 4 assigned a
  clone only to cells its SCEVAN run placed, so 178 malignant cells carry no clone and are
  NOT EVALUABLE for every Phase 6 clone analysis.
* `cnv_burden` and the other CNA metrics are non-NA for 18,449 cells — the 1,267 SCEVAN-filtered
  cells were recorded as *not assessed*, never as non-malignant.

## 4. Patient imbalance — severe, and the reason a balanced analysis is mandatory

| patient | malignant cells | share | median nCount | median nFeature | High | Moderate | clones | clone-reliable |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--: |
| **MPNST_4** | **3,685** | **57.3%** | 5,250 | 2,249 | 936 | 2,749 | 8 | yes |
| MPNST_1 | 1,886 | 29.3% | 16,034 | 5,361 | 1,861 | 25 | 7 | yes |
| MPNST_2 | 651 | 10.1% | 18,901 | 3,978 | 464 | 187 | 4 | yes |
| MPNST_3 | 212 | 3.3% | 9,327 | 3,299 | 0 | 212 | 3 | **no** |

**Imbalance ratio 17.4×.** MPNST_4 alone holds 57.3% of the malignant compartment, and it is also
the shallowest library (median 5,250 counts against MPNST_2's 18,901). Any factorization run on
these cells will be exposed to both effects at once, which is why the patient-balanced sensitivity
analysis is a requirement of the design rather than an optional extra.

Confidence over the malignant compartment: 3,261 High, 3,173 Moderate, 0 Low.

## 5. Phase 4 discrete states, for the M38 comparison

| state | MPNST_1 | MPNST_2 | MPNST_3 | MPNST_4 |
| --- | ---: | ---: | ---: | ---: |
| Mesenchymal_ECM-1 | **1,851** | 0 | 0 | 0 |
| Mesenchymal_ECM-2 | 1 | 3 | 7 | **1,587** |
| Mesenchymal_ECM-3 | 0 | 0 | 0 | **1,291** |
| Mesenchymal_ECM-4 | 0 | **645** | 1 | 2 |
| Mesenchymal_ECM-5 | 0 | 0 | 0 | **505** |
| Cycling | 34 | 3 | 12 | 241 |
| Schwann_like | 0 | 0 | 120 | 59 |
| Interferon | 0 | 0 | 72 | 0 |

Each of the five `Mesenchymal_ECM` states is effectively one patient's: ECM-1 is MPNST_1, ECM-4 is
MPNST_2, and ECM-2/-3/-5 are all MPNST_4. That is the concrete form of Phase 4's negative result,
and it is the structure M38 tests against the continuous programs.

## 6. Expression basis and gene availability

RNA counts across four per-sample layers (`counts.MPNST_1.1` … `counts.MPNST_4.4`), integer-verified,
joined for Phase 5. Cell order is consistent between the object and its metadata.

30,762 of 31,764 genes are detected in at least one malignant cell:

| detected in ≥ | genes |
| --- | ---: |
| 0.1% of malignant cells | 25,329 |
| **0.5%** | **19,676** |
| 1% | 17,413 |
| 5% | 12,377 |

95 canonical Tirosh S/G2M cell-cycle genes are present, retained in the primary run and removed
only in the declared sensitivity run.

## 7. Working artefacts exported, so the 6 GB object is opened once

| artefact | contents |
| --- | --- |
| `phase5_cell_metadata.rds` | all 19,716 cells × 183 Phase 1–4 metadata columns (4.3 MB) |
| `phase5_working_counts.rds` | RNA counts, 31,764 genes × **7,462** cells (malignant ∪ historical fibroblast) |
| `phase5_working_lognorm.rds` | the frozen RNA `data` layer for the same cells |
| `phase5_malignant_gene_detection.rds` | per-gene detection over the malignant compartment |
| `phase5_cellcycle_genes.txt` | the 95 cell-cycle genes actually present |

The working set is malignant (6,434) ∪ historical Fibroblast (5,064), overlapping in 4,036 — so
M39's non-malignant and Ambiguous fibroblasts are available without reopening the Phase 4 object.

## 8. Feasibility verdict

**Phase 5 is feasible as specified**, with three constraints that the design must carry rather than
discover later:

1. **Patient imbalance is severe (57.3% / 29.3% / 10.1% / 3.3%)** and is confounded with sequencing
   depth. The patient-balanced run is therefore mandatory, and the smallest patient (MPNST_3, 212
   cells) is too small to serve as the downsampling cap — the cap is MPNST_2's 651.
2. **MPNST_3 contributes only 212 malignant cells and 0 High-confidence ones.** It can enter program
   discovery but cannot carry weight in any per-patient conclusion, and its SCEVAN clone structure
   is already excluded from Phase 6 by the frozen Phase 4 reliability flag.
3. **178 malignant cells have no clone label.** They are NOT EVALUABLE for Phase 6 clone analyses
   and are reported as such rather than dropped silently.

Facts: `results/phase5/validation/m36_feasibility_facts.json`.
