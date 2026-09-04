# M39 — Malignant ECM-like MPNST Cells versus True Fibroblasts

**Date** 2026-09-03 · **Status** COMPLETE · **SLURM** 19899354

**Scientific question.** What transcriptional features distinguish CNA-associated malignant ECM-like
MPNST cells from genuine non-malignant fibroblasts, despite shared extracellular-matrix expression?

## Groups, from the field M36 verified

Grouping uses **`annotation_ccc_phase3`** — the frozen pre-Phase-4 annotation.
`annotation_ccc_refined` was barred because it already encodes the Phase 4 conclusion.

```text
A  historical Fibroblast + malignancy_refined == Malignant       4,036
B  historical Fibroblast + malignancy_refined == Non-malignant     908
C  historical Fibroblast + malignancy_refined == Ambiguous         120   (projection only)
```

All three counts re-asserted against the frozen Phase 4 values.

## Evaluability was checked BEFORE any test was run

| patient | A malignant ECM-like | B true fibroblast | evaluable | reason |
| --- | ---: | ---: | :--: | --- |
| MPNST_1 | 512 | 303 | **yes** | evaluable |
| MPNST_2 | 637 | 503 | **yes** | evaluable |
| MPNST_3 | 0 | 101 | no | no malignant ECM-like cells — Phase 4 disabled its promotions |
| MPNST_4 | 2,887 | **1** | no | too few non-malignant fibroblasts to compare against |

Declared minimum: 30 cells on **both** sides. **A four-patient paired test was not constructed** —
the groups do not coexist in 2 of the 4 patients, and inventing one would have been the kind of
false precision the phase brief prohibits.

## Effect sizes, per patient

Pseudobulk log2 fold change (A vs B) was computed per evaluable patient over 21,864 genes, and
cell-level AUC per gene per patient with `presto`. **The cross-patient Spearman correlation of the
pseudobulk log2FC is only 0.091** — most of the A-vs-B difference is patient-specific. That number
is reported prominently because it bounds the whole comparison: a small consistent core exists
inside a largely patient-specific difference.

## The transparent signature

Declared criteria, applied in **every** evaluable patient: |pseudobulk log2FC| ≥ 1, cell-level AUC
≥ 0.65 (or ≤ 0.35), detected in ≥ 10% of the higher group. **No classifier was trained and no
random cell-level train/test split was used**, so no pseudo-independent accuracy is reported.

**61 genes up in malignant ECM-like cells** — CA12, SCG2, IGFBP3, TMEM176A, TMEM176B, TUBB2B, CAPG,
COL14A1, ACKR3, IGFBP2, GJA1, CTSH, ANGPT1, CPE, FOXP1, TNFSF9, COL11A1, ENPP2, EFNA5, S100A10,
CCDC3, BASP1, PTGIS, SPOCK1, BOC …

**31 genes up in true fibroblasts** — CDH19, APOD, SCN7A, ABCA6, MCTP1, ABCA10, VIT, PAMR1, FOXS1,
SPARCL1, ABCA9, GPC3, ABCA8, BST2, GSN, A2M, BAMBI, SLC12A2, SYNE2, LAMA2, SRPX, GAS7, CEBPD, CFH,
NID1 …

The fibroblast side is biologically coherent in a way worth stating: **CDH19, SCN7A, APOD and the
ABCA6/8/9/10 cluster are markers of nerve-associated / endoneurial fibroblasts and non-myelinating
Schwann-lineage stroma** — exactly the resident population a peripheral-nerve-sheath tumour would be
expected to retain. The malignant side is dominated by IGFBP2/IGFBP3, COL11A1 and COL14A1, GJA1 and
SPOCK1.

## Programs and CNA burden across the groups

The 5,064 historical fibroblasts — including the 908 Non-malignant and 120 Ambiguous cells that
were deliberately excluded from program discovery — were **projected** onto the fixed cNMF spectra
by non-negative least squares with the spectra held constant, so a projected cell cannot influence
a program.

Median program usage difference (A − B) in the evaluable patients shows the separation is
patient-specific in *which* program carries it: MPNST_1 separates on P4 (+0.337) and P5 (+0.246)
while MPNST_2 separates on P3 (+0.495). The direction of separation is consistent; the program is
not.

Phase 4 CNA burden by group, reported as **context, not as independent validation of the split that
defined the groups**: MPNST_1 median `cnv_burden` 0.369 (A) vs 0.184 (B); MPNST_2 0.311 vs 0.186;
MPNST_4 0.262 vs 0.230 (n = 1, not evaluable).

## Ambiguous fibroblasts

The 120 `Ambiguous` historical fibroblasts were projected and tabulated in
`AMBIGUOUS_FIBROBLAST_PROGRAM_PROJECTION.tsv`. **They were not reclassified** — an assertion in the
script confirms every one of them still carries `malignancy_refined == "Ambiguous"`.

**Outputs.** `MALIGNANT_ECM_SIGNATURE.tsv` (92 rows) · `MALIGNANT_ECM_VS_TRUE_FIBROBLAST.tsv`
(21,864 rows) · `AMBIGUOUS_FIBROBLAST_PROGRAM_PROJECTION.tsv` ·
`M39_EVALUABILITY_BY_PATIENT.tsv` · `M39_PROGRAM_USAGE_A_VS_B.tsv` ·
`M39_CNV_METRICS_BY_GROUP.tsv` · `M39_FIBROBLAST_PROGRAM_PROJECTION_ALL.tsv`

**Warnings.** The comparison rests on **2 of 4 patients**, and the cross-patient log2FC correlation
is 0.091. The signature is the intersection that survives both patients, not a general result.

**Next.** M40 — robustness.
