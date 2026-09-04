# Phase 4 — Malignancy Refinement Report

**Milestone M30 · SLURM 19897148 · COMPLETED · 00:01:00 · MaxRSS 1.38 GiB**

> SCEVAN infers copy number from gene-expression patterns. **This is not DNA sequencing.**
> Inference is reliable for broad chromosomal, arm-level and large-segment events, not for
> single genes. A SCEVAN non-malignant call is **not proof** of non-malignancy: some MPNST
> malignant cells may be copy-number quiet.

---

## 1. The headline answer

### Was the Phase 2 estimate of 17.35% too conservative?

**Yes — but not in the direction of a simple upward correction.**

| Stage | Malignant cells | Fraction of 19,716 |
| --- | ---: | ---: |
| **Phase 2, conservative** (marker evidence only) | 3,420 | **17.35%** |
| **SCEVAN primary** (copy-number evidence only) | 8,794 of 18,449 assessed | **47.7% of assessed** · 44.6% of all |
| **Phase 4 refined** (integrated) | **6,434** | **32.63%** |
| — of which **High** confidence | 3,261 | 16.54% |
| — of which Moderate confidence | 3,173 | 16.09% |

Alongside: **9,078 Non-malignant**, **3,766 Ambiguous (19.10%)**, 438 Excluded-low-quality.

The refined malignant compartment is **1.88× the Phase 2 estimate**. But that headline hides a
two-directional correction, and both directions matter:

* Phase 4 **adds** 4,036 `Fibroblast`, 836 `Candidate-Malignant-Unresolved`, 92 `Uncertain` and
  65 `Pericyte-VSMC` cells that Phase 2 had classified as non-malignant or unresolved.
* Phase 4 **withdraws confidence from** 2,015 of the 3,420 cells Phase 2 called `MPNST-Tumor`,
  moving them to `Ambiguous` because SCEVAN calls them copy-number normal.

So the correct statement is not "Phase 2 undercounted". It is: **Phase 2 looked in the wrong
place.** The malignant compartment is roughly twice as large as Phase 2 thought, and it is
substantially composed of cells Phase 2 read as stroma — while a majority of the cells Phase 2
did call malignant are not corroborated by copy-number evidence.

---

## 2. Per-patient refined malignancy

| Patient | cells | Phase 2 malignant | SCEVAN malignant | **refined Malignant** | Non-malignant | Ambiguous | Excluded |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| MPNST_1 | 7,615 | 2,428 | 1,952 | **1,886** | 3,303 | 2,415 | 11 |
| MPNST_2 | 2,284 | 97 | 694 | **651** | 1,500 | 127 | 6 |
| MPNST_3 | 2,940 | 231 | 1,414 | **212** | 1,792 | 936 | 0 |
| MPNST_4 | 6,877 | 664 | 4,734 | **3,685** | 2,483 | 288 | 421 |

MPNST_3's refined count (212) is *lower* than its Phase 2 count (231) because its SCEVAN
partition failed the sanity gate and malignant promotions from it were disabled. That is the
amendment working as designed, not a biological result.

---

## 3. Threshold sensitivity — the knife-edges, disclosed

`pop_frac_low` was fixed at **0.25 a priori**. Two strata sit close to it, so the headline was
recomputed across a grid. **The 0.25 row is the reported result; the grid is disclosure, not
tuning.**

| `pop_frac_low` | Malignant | fraction | Fibroblast malignant | `MPNST-Tumor` retained | Candidate resolved |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.15 | 8,295 | 42.07% | **4,036** | 3,266 | 836 |
| 0.20 | 8,295 | 42.07% | **4,036** | 3,266 | 836 |
| **0.25 (a priori)** | **6,434** | **32.63%** | **4,036** | **1,405** | **836** |
| 0.30 | 6,434 | 32.63% | **4,036** | 1,405 | 836 |
| 0.40 | 6,434 | 32.63% | **4,036** | 1,405 | 836 |

Two things follow, and they point in opposite directions:

* **The fibroblast finding is completely threshold-independent** — 4,036 malignant fibroblasts
  at every threshold from 0.15 to 0.40. It does not rest on the choice of cut-off.
* **The `MPNST-Tumor` retention is entirely threshold-dependent** — 3,266 versus 1,405 across
  the 0.20/0.25 boundary, because MPNST_1's `MPNST-Tumor` pop_frac is **0.222**, just under
  the a priori 0.25. Had the threshold been 0.20, those 1,861 cells would have been
  `Malignant/Moderate` rather than `Ambiguous/Low`. This is the single most fragile number in
  Phase 4 and it should be read as such.

---

## 4. Which rules fired

Full coverage: **R99 = 0 cells**, so no cell reached the fail-safe and there is no rule gap.

| Rule | Outcome | n | Reading |
| --- | --- | ---: | --- |
| R13 | Malignant / High | 1,227 | **disputed stromal promoted on strong CNA evidence** |
| R14 | Malignant / Moderate | 2,874 | disputed stromal promoted, one criterion short |
| R17 | Malignant / High | 907 | `Candidate-Malignant-Unresolved` / `Uncertain` resolved as malignant |
| R2 | Malignant / High | 1,127 | Phase 2 and SCEVAN agree, both reference strategies |
| R4 | Malignant / Moderate | 178 | Phase 2 malignant, SCEVAN normal, population predominantly malignant |
| R3 | Malignant / Moderate | 100 | agree on call, sensitive to reference choice |
| R18 | Malignant / Moderate | 21 | provisional promoted on weaker support |
| R6 | Non-malignant / High | 4,711 | concordant immune |
| R16 | Non-malignant / Moderate | 811 | stromal called normal inside a mixed stratum |
| R9 | Non-malignant / High | 741 | concordant endothelium |
| R19 | Non-malignant / Moderate | 476 | provisional resolved as non-malignant |
| R7p | Non-malignant / Moderate | 1,125 | **plasma-cell immunoglobulin-locus artefact** (Amendment A2) |
| R1a | Non-malignant / Moderate | 742 | SCEVAN-filtered but canonical immune |
| R12 | Non-malignant / High | 445 | **true fibroblasts / pericytes** |
| R7, R10 | Non-malignant / Moderate | 27 | minority malignant call inside a canonical population |
| **R5** | **Ambiguous / Low** | **1,933** | **Phase 2 malignant, SCEVAN normal — genuine disagreement, preserved** |
| R8 | Ambiguous / Low | 851 | majority-malignant immune population (almost all MPNST_3) |
| R1b | Ambiguous / Low | 445 | no CNV evidence and no decisive lineage evidence |
| R20 | Ambiguous / Low | 326 | provisional, mixed stratum, stays unresolved |
| R11 | Ambiguous / Low | 168 | majority-malignant endothelium (MPNST_4) |
| R13x, R17x | Ambiguous / Low | 37 | promotions disabled — unreliable sample |
| R15 | Ambiguous / Low | 6 | isolated malignant calls in a normal stromal population |
| R0 | Excluded-low-quality | 438 | mitochondrial-dominated technical artefact |

`R2u` fired for **0** cells: every Phase 2-malignant cell with SCEVAN agreement turned out to
be in a sanity-passing sample.

---

## 5. Disagreement, preserved rather than resolved (§25)

| Situation | Cells | What Phase 4 did |
| --- | ---: | --- |
| Phase 2 `MPNST-Tumor` + SCEVAN malignant | 1,227 | **Malignant / High** — strong concordance |
| Phase 2 `Fibroblast` + SCEVAN normal | 908 | **Non-malignant** — supports a true fibroblast |
| Phase 2 `Fibroblast` + SCEVAN malignant | 4,036 | **Malignant** — candidate Mes-NC-like population |
| Phase 2 `MPNST-Tumor` + SCEVAN normal | 2,015 | **Ambiguous** — flagged for investigation, not overridden |

Neither method was forced to match the other, and no disagreement was averaged away.

---

## 6. Guards

```text
GUARD PASSED: population_class assigned to every cell
GUARD PASSED: no Ambiguous or Non-malignant cell is labelled MPNST-Tumor
GUARD PASSED: annotation_ccc_phase3 is a verbatim copy of the frozen Phase 3 layer
R99 fail-safe: 0 cells  (complete rule coverage)
```

`Ambiguous` cells are **never** promoted to tumour, in this milestone or in
`annotation_ccc_refined`.

---

## 7. Refined CCC annotation

`annotation_ccc_refined` alongside the verbatim `annotation_ccc_phase3`:

| Population | MPNST_1 | MPNST_2 | MPNST_3 | MPNST_4 | total |
| --- | ---: | ---: | ---: | ---: | ---: |
| **MPNST-Tumor** | 1,886 | 651 | 212 | 3,685 | **6,434** |
| Ambiguous-unresolved | 2,415 | 127 | 936 | 288 | 3,766 |
| Macrophage | 1,133 | 361 | 704 | 867 | 3,065 |
| Plasma-cell | 416 | 72 | 484 | 737 | 1,709 |
| Fibroblast | 303 | 503 | 101 | 1 | **908** |
| Endothelial | 335 | 112 | 151 | 151 | 749 |
| Dendritic | 100 | 168 | 144 | 200 | 612 |
| Uncertain | 379 | 97 | 0 | 0 | 476 |
| Low-quality-excluded | 11 | 6 | 0 | 421 | 438 |
| CD4-T | 64 | 16 | 38 | 161 | 279 |
| CD8-T | 146 | 26 | 8 | 67 | 247 |
| Pericyte-VSMC | 250 | 26 | 68 | 4 | 348 |
| B-cell | 54 | 15 | 4 | 71 | 144 |
| T-cell-other | 43 | 18 | 11 | 83 | 155 |
| Monocyte | 37 | 27 | 72 | 54 | 190 |
| Plasmacytoid-DC | 30 | 48 | 0 | 44 | 122 |
| NK | 13 | 11 | 7 | 43 | 74 |

`Ambiguous-unresolved`, `Uncertain`, `Low-quality-excluded` are **excluded** from the refined
CCC analysis as not CCC-ready. Excluded is **NOT EVALUABLE**, which is not the same as "no
signalling".

Note the consequence for downstream comparability: `Fibroblast` now has **1 cell in MPNST_4**
and `Pericyte-VSMC` **4**, so those populations become NOT EVALUABLE in that patient, and
`Plasmacytoid-DC` and `Uncertain` vanish from MPNST_3/MPNST_4. The M32 sensitivity analysis
therefore distinguishes "lost" from "no longer testable".

---

## 8. Figures

`results/phase4/figures/` — `01_phase2_vs_scevan_malignancy` · `02_scevan_malignancy_umap` ·
`03_refined_malignancy_umap` · `04_malignancy_by_phase2_cluster` ·
`05_malignancy_by_annotation` · `06_malignancy_by_patient` ·
`07_ambiguous_population_malignancy` · `08_cnv_profiles_by_annotation` (PDF + PNG each).

## 9. Tables

`results/phase4/malignancy/PHASE4_MALIGNANCY_CALLS.tsv` (19,716 × 37, one row per cell with
its rule and printed reason) · `SCEVAN_VS_PHASE2_ANNOTATION.tsv` ·
`AMBIGUOUS_POPULATION_MALIGNANCY_AUDIT.tsv` · `MALIGNANCY_SUMMARY_BY_CLUSTER.tsv` ·
`MALIGNANCY_SUMMARY_BY_SAMPLE.tsv` · `SCEVAN_SAMPLE_RELIABILITY.tsv` ·
`SCEVAN_POPULATION_MALIGNANT_FRACTION.tsv` · `MALIGNANCY_THRESHOLD_SENSITIVITY.tsv` ·
`LINEAGE_MARKER_EVIDENCE_BY_POPULATION.tsv` · `SCEVAN_CELL_CLASSIFICATION.tsv`.
