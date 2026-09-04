# M29 — Per-Patient SCEVAN Execution

**Status: COMPLETE** · 2026-09-03 · all four patients, both reference strategies, 8/8 runs SUCCESS

## Design

SCEVAN was run **independently for each patient** (§14), because `sample_id` = patient =
dataset in this project and pooling first would confuse patient-specific CNV architecture with
intra-tumour subclones. Each patient got **two** runs:

* **primary** — `norm_cell = NULL`. SCEVAN's own confident-normal detection, with **no input
  from Phase 2 labels**. This is the run `malignancy_scevan` reports.
* **sensitivity** — the 7,448 high-confidence immune cells as `norm_cell`, with
  `FIXED_NORMAL_CELLS = FALSE` (never `TRUE`; see `SCEVAN_FEASIBILITY.md` §5).

Identical parameters across all eight runs: `SUBCLONES = TRUE`, `ClonalCN = TRUE`,
`plotTree = TRUE`, `beta_vega = 0.5`, `organism = "human"`, `ngenes_chr = 5`,
`perc_genes = 10`, `par_cores = 8`, seed 42.

## Results

| Patient | cells | assessed | **malignant** | non-malignant | filtered | subclones | confident normals | primary↔sensitivity agreement |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| MPNST_1 | 7,615 | 7,043 | **1,952** | 5,091 | 572 | 7 | 30 | **0.9951** |
| MPNST_2 | 2,284 | 2,205 | **694** | 1,511 | 79 | 4 | 30 | **0.9663** |
| MPNST_3 | 2,940 | 2,749 | 1,414 | 1,335 | 191 | 3 | 25 | **0.0864** |
| MPNST_4 | 6,877 | 6,452 | **4,734** | 1,718 | 425 | 8 | 30 | **0.9974** |
| **total** | **19,716** | **18,449** | **8,794** | 9,655 | 1,267 | 22 | 115 | — |

`filtered` cells are those SCEVAN dropped in preprocessing (<200 genes, or <5 genes on some
chromosome). They are recorded as **not-assessed**, never as non-malignant.

The reference strategy makes almost no difference in three of four patients (0.966–0.997),
which is the point of running both: the primary result is not an artefact of how the baseline
was chosen. MPNST_3 is the exception and is treated below.

## SCEVAN calls against the frozen Phase 2 annotation

Malignant fraction **of assessed cells**, primary run:

| Phase 2 `annotation_ccc` | MPNST_1 | MPNST_2 | MPNST_3 | MPNST_4 |
| --- | ---: | ---: | ---: | ---: |
| **Fibroblast** | **0.628** | **0.559** | 0.082 | **1.000** |
| **Candidate-Malignant-Unresolved** | **0.733** | 0.167 | 0.789 | 0.429 |
| MPNST-Tumor | 0.225 | 0.089 | 0.429 | **0.912** |
| Pericyte-VSMC | 0.016 | 0.071 | 0.081 | **0.942** |
| Endothelial | 0.000 | 0.009 | 0.046 | **0.527** |
| Uncertain | 0.008 | 0.058 | 0.438 | **0.988** |
| CD4-T / CD8-T / NK / T-other | 0.000 | 0.000 | ~0.99 | ≤0.031 |
| Macrophage / Monocyte / Dendritic | 0.000 | 0.000 | ≤0.055 | 0.000 |
| B-cell | 0.000 | 0.000 | 1.000 | 0.036 |
| Plasma-cell | 0.197 | 0.635 | 0.993 | 0.982 |
| Low-quality-excluded | 0.000 | 0.000 | — | 0.950 |

Three things are visible here and all three are reported, not just the convenient one.

**1. The fibroblast signal is real and replicated.** 62.8% and 55.9% malignant in MPNST_1 and
MPNST_2 — the two patients that pass every version of the sanity gate — and 100% in MPNST_4.
This is the direct answer to the Phase 4 question in §22 and §43.

**2. SCEVAN does not corroborate most cells Phase 2 called malignant.** `MPNST-Tumor` is only
22.5% and 8.9% malignant in MPNST_1 and MPNST_2. Two readings are available and the data here
cannot separate them: some MPNST malignant cells are **copy-number quiet** (§65C, which is
exactly why a SCEVAN normal call is not proof of non-malignancy), or SCEVAN under-calls in
samples where the malignant compartment is small. This disagreement is preserved by rules
R4/R5 rather than resolved by fiat, and it is the reason the refined malignant fraction is not
simply "Phase 2 plus fibroblasts".

**3. Two samples fail a basic sanity check, for different reasons.**

## MPNST_3 — a genuine failure, reported as one

Both MPNST_3 runs assign large blocks of unambiguous immune cells to the malignant class, in
**opposite directions**:

```text
              primary        sensitivity
B-cell        89/89 tumour   0/89
CD4-T        327/328 tumour  4/328
CD8-T        135/136 tumour  0/136
NK            32/32 tumour   0/32
Plasmacytoid-DC 165/165      0/165
Macrophage     4/680         678/680 tumour
Dendritic      8/144         142/144 tumour
Monocyte       1/72          72/72 tumour
```

Agreement 0.0864. SCEVAN found only **25 confident normal cells** there, and MPNST_3 has the
lowest depth of the four (median 4,265 counts, **1,594 genes** per cell). With the baseline
estimated from 25 cells, `classifyCluster()` partitions the top-level CNA clustering along a
**lymphoid/myeloid** axis rather than a malignant/normal one, and the two runs label those
clusters oppositely depending on which reference they started from.

**Six of nine** testable immune populations individually exceed 25% malignant. This run's
malignant/normal partition is structurally wrong, and **malignant promotions from MPNST_3 are
disabled** (`MALIGNANCY_DECISION_RULES.md` Amendments A1/A2, rules R13x and R17x). The cost is
stated plainly: any Mes-NC-like malignant fibroblast population in MPNST_3 cannot be detected
by this analysis.

## MPNST_4 — one artefact-prone population, not a failed run

MPNST_4's high overall immune rate (0.293) is driven almost entirely by **Plasma-cell:
608 of 619 assessed, 98.2%**. Every other immune population is clean (B-cell 0.036, CD4-T
0.007, CD8-T 0.000, Dendritic 0.000, Macrophage 0.000, Monocyte 0.000, NK 0.024, pDC 0.000,
T-cell-other 0.031); **zero of nine** populations fail; and primary↔sensitivity agreement is
**0.9974**, the highest of the four.

Plasma cells express immunoglobulin loci (IGH 14q32, IGK 2p11, IGL 22q11) at extraordinary
levels, which expression-derived CNV inference reads as large segmental gain. That is a
documented artefact class for this family of methods, and the same signature appears at lower
magnitude in the passing samples too (MPNST_1 0.197, MPNST_2 0.635) — evidence that it is a
property of the population, not of the sample.

Amendment A2 therefore excludes plasma cells from the sanity denominator, adds a **breadth**
criterion so a global inversion cannot dilute itself past the gate, and adds rule **R7p** so
those plasma cells are called `Non-malignant` with the artefact named rather than dropped into
`Ambiguous`. **This amendment was written after observing that the flat gate excluded MPNST_4,
and the rules document says so explicitly**, along with the two independent checks a sceptical
reader should apply.

## Subclone structure

22 subclones across the four patients — MPNST_1 7, MPNST_2 4, MPNST_3 3, MPNST_4 8. Clone
labels are **patient-scoped** by construction (`<PATIENT>_cloneN`): a clone in one patient is
not the same clone as one in another, and cross-patient similarity is described as a **shared
or recurrent CNA pattern** (§50).

## SLURM accounting

| JobID | Patient | State | Elapsed | AllocCPUS | ReqMem | MaxRSS |
| --- | --- | --- | --- | ---: | --- | --- |
| 19896712_1 | MPNST_1 | COMPLETED | 00:34:20 | 8 | 64 G | **60.87 GiB** |
| 19896712_2 | MPNST_2 | COMPLETED | 00:08:13 | 8 | 64 G | 6.24 GiB |
| 19896712_3 | MPNST_3 | COMPLETED | 00:11:20 | 8 | 64 G | 18.10 GiB |
| 19896712_4 | MPNST_4 | COMPLETED | 00:42:07 | 8 | 64 G | 50.49 GiB |

**No failures, no retries, no resource increases.** Worth recording honestly: MPNST_1 peaked at
**95% of its 64 G request**. The estimate came from the M28 dense-matrix measurement (1.71 GiB)
multiplied for SCEVAN's smoothed/relative/segmented copies plus `parallelDist`, and it was
closer than intended. It did not need raising because it did not fail, but a larger sample
would need a larger request, and that is now measured rather than guessed.

All four tasks ran concurrently: 4 × 8 = 32 CPUs (the full node allocation) and 4 × 64 G =
256 G of the 450 G envelope.

## Warnings

Every run emitted a caught error from `plotCloneTree` — ggtree calls `ggplot2:::is.waive()`,
removed in ggplot2 4.x. SCEVAN wraps it in its own `tryCatch`, so all eight runs returned
normally and **every scientific output is intact**: subclone assignments, `*_Clonal_CN.seg`,
`*_subcloneN_CN.seg`, `*_CNAmtx.RData`, `*_CNAmtxSubclones.RData`, `*PlotOncoHeat.RData`, CNA
heatmaps, consensus profiles and onco-heatmaps. Only the clone **phylogeny figure** is missing.
**ggplot2 was deliberately not downgraded** — that would destabilise the Phase 3 figure suite,
and a dendrogram is not a Phase 4 deliverable. Recorded as limitation K.

## Scientific decisions

1. Per-patient runs, never pooled (§14).
2. Primary run non-circular by construction (`norm_cell = NULL`).
3. `FIXED_NORMAL_CELLS = TRUE` prohibited project-wide.
4. Both reference strategies run and their agreement reported per patient — the reason MPNST_3's
   failure was detectable at all.
5. Disagreement with Phase 2 preserved, in both directions.
6. Two sanity amendments, both dated and both labelled as post-hoc.

## Next milestone

**M30 — malignancy integration**: apply the a priori rules plus Amendments A1/A2, produce
`PHASE4_MALIGNANCY_CALLS.tsv`, the Phase 2 × SCEVAN cross-tabulation, the ambiguous-population
audit, and the eight §27 figures.
