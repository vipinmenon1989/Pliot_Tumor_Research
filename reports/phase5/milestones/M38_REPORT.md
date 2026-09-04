# M38 — Program Annotation, Pathways and the Comparison with the Phase 4 States

**Date** 2026-09-03 · **Status** COMPLETE · **SLURM** 19899354, re-run as 19899356 after adding
explicit technical-content detection to the labels

**Scientific question.** What are the eight programs, which are recurrent across patients, and are
the five patient-private `Mesenchymal_ECM` states manifestations of one shared continuous program?

## The eight programs at K = 8

| id | label | top genes (first 8) | top gene set | technical? |
| --- | --- | --- | --- | :--: |
| P1 | **Translation_ribosomal** | RPL10 RPL19 RPS2 RPLP1 EEF1A1 RPS12 RPS4Y1 RPS4X | REACTOME_EUKARYOTIC_TRANSLATION_INITIATION | **yes** (54% of top 50 ribosomal) |
| P2 | Unassigned_NEURONAL_SYSTEM | TSHZ2 AUTS2 TENM2 AFF3 FRMD4A CHN1 NLGN1 CNTN4 | REACTOME_NEURONAL_SYSTEM | no |
| P3 | **Mesenchymal_ECM_mixed** | SERPINF1 CST3 FBLN1 C1R SFRP2 EFEMP1 CLU NNMT | HALLMARK_TNFA_SIGNALING_VIA_NFKB | no |
| P4 | Hypoxia_Angio | GBE1 IGF1R ADAMTS9-AS2 OSBPL3 VEGFA RNF144A | HALLMARK_UV_RESPONSE_DN | no |
| P5 | **Translation_ribosomal** | RPS7P1 RPS7P10 RPS3AP6 RPL3P4 RPL6P27 | REACTOME_TRANSLATION | **yes** (36% ribosomal, 66% pseudogene/lncRNA) |
| P6 | Schwann_like_mixed | S100B AP1S2 GPM6B SERPINE2 CRYAB S100A6 MIA | HALLMARK_INTERFERON_GAMMA_RESPONSE | no |
| P7 | **Cycling** | NUSAP1 PRC1 TOP2A UBE2C ASPM HMGB2 TPX2 MKI67 | REACTOME_CELL_CYCLE_MITOTIC | no |
| P8 | **Myeloid_ambient_like** | C1QC C1QB C1QA CD14 RNASE1 FOLR2 MS4A4A AIF1 | REACTOME_NEUTROPHIL_DEGRANULATION | **yes** (18 myeloid markers in top 50) |

**Three of eight programs are technical-dominated**, and they are named for what they are. The
first pass of the labelling rule called P1 `Mesenchymal_ECM` on a 0.20 marker-family overlap while
its own top 15 genes were ribosomal proteins and its top gene set was translation initiation. That
label would have misled every downstream reader, so an explicit technical-content check was added
(fraction of the top 50 genes that are ribosomal proteins, pseudogenes/lncRNAs, or canonical
myeloid markers) and M38 was re-run. **Nothing about the factorization, K selection or the
recurrence rule changed — only how a program is named once its own genes are inspected.**

## Recurrence — the headline result

Criteria declared before any program was labelled, and deliberately the same shape as Phase 4's
state-recurrence rule so the two are comparable: a cell is program-active at ≥ 0.20 relative usage;
a patient carries a program with ≥ 10 active cells **and** ≥ 5% of its malignant cells active;
recurrent needs ≥ 3 of 4 patients.

| program | active cells | carrying patients | n | dominant patient | dominant fraction | status |
| --- | ---: | --- | ---: | --- | ---: | --- |
| P1 Translation_ribosomal | 3,720 | MPNST_2, MPNST_4 | 2 | MPNST_4 | 0.974 | **shared-limited** |
| P2 Neuronal | 1,286 | MPNST_1 | 1 | MPNST_1 | 0.997 | patient-private |
| P3 Mesenchymal_ECM | 670 | MPNST_2 | 1 | MPNST_2 | 0.972 | patient-private |
| P4 Hypoxia_Angio | 973 | MPNST_1 | 1 | MPNST_1 | 0.996 | patient-private |
| P5 Translation_ribosomal | 828 | MPNST_1 | 1 | MPNST_1 | 1.000 | patient-private |
| P6 Schwann_like | 184 | MPNST_3 | 1 | MPNST_3 | 0.674 | patient-private |
| P7 Cycling | 229 | MPNST_4 | 1 | MPNST_4 | 0.873 | patient-private |
| P8 Myeloid_ambient | 89 | MPNST_3 | 1 | MPNST_3 | 0.989 | patient-private |

**0 recurrent · 1 shared-limited · 7 patient-private.** And the single shared-limited program is
one of the technical ones — **not a single biologically interpretable program is carried by even
two patients at the declared threshold.**

Dominant-program assignment is equally stark: P1 → 3,490 MPNST_4 cells; P2 → 1,002 MPNST_1;
P3 → 651 MPNST_2 (every malignant cell in that patient); P4 → 581 MPNST_1; P5 → 301 MPNST_1;
P6 → 122 MPNST_3 + 55 MPNST_4; P7 → 139 MPNST_4 + 3 MPNST_3 + 2 MPNST_1; P8 → 83 MPNST_3.

**Threshold sensitivity, reported across a grid rather than at one cut.** At the relaxed 0.10
activity threshold with the same 5% carriage rule, **two programs reach three patients: P1
(technical) and P3 (Mesenchymal_ECM)**. P3's per-patient active fractions at 0.10 are MPNST_2 1.000,
MPNST_4 0.108, MPNST_3 0.085, MPNST_1 0.007. So the ECM/mesenchymal program is the one
biologically interpretable program that comes close to recurrence — carried by 3 of 4 patients at a
10% relative-usage threshold, by 1 at the pre-declared 20%. That distinction is the honest result
and is reported as such rather than resolved in either direction.

## Patient-balanced support

Every one of the eight programs is recovered in the balanced run (dominant-patient fraction cut
from 0.573 to 0.301): best-match spectra cosine 0.50–0.98, top-100 gene Jaccard 0.01–0.87, all
eight above the declared 0.60 cosine **or** 0.25 Jaccard support threshold. **Patient imbalance
does not create these programs — but it does confine them**, which is a different statement and the
one the data support.

## Are the five Mesenchymal_ECM states one shared program?

| Phase 4 state | patient | dominant program | median relative usage |
| --- | --- | --- | ---: |
| Mesenchymal_ECM-1 | MPNST_1 (1,851) | **P2** Neuronal | 0.460 |
| Mesenchymal_ECM-2 | MPNST_4 (1,587) | **P1** Translation_ribosomal | 0.867 |
| Mesenchymal_ECM-3 | MPNST_4 (1,291) | **P1** | 0.894 |
| Mesenchymal_ECM-4 | MPNST_2 (645) | **P3** Mesenchymal_ECM | 0.887 |
| Mesenchymal_ECM-5 | MPNST_4 (505) | **P1** | 0.877 |

**Answer: no — and the structure is informative in both directions.** *Within* MPNST_4, three
separately-clustered ECM states (3,383 cells) collapse onto one continuous program, so Phase 4's
discrete partition was finer there than the continuous structure warrants. *Across* patients, the
ECM states map to three different programs, so the five states are not one shared biology.
**Possible result C of the three the brief anticipated: ECM programs remain substantially
patient-private.**

**Outputs.** `MALIGNANT_PROGRAMS.tsv` · `MALIGNANT_PROGRAM_TOP_GENES.tsv` (400 rows) ·
`MALIGNANT_PROGRAM_CELL_SCORES.tsv` (6,434 × 16) · `MALIGNANT_PROGRAM_PATIENT_DISTRIBUTION.tsv` ·
`MALIGNANT_PROGRAM_PATHWAYS.tsv` · `PROGRAM_VS_PHASE4_STATE.tsv` ·
`PROGRAM_PATIENT_BALANCE_SENSITIVITY.tsv` · `M38_PROGRAM_TECHNICAL_CONTENT.tsv` ·
`M38_PROGRAM_MARKER_FAMILY_OVERLAP.tsv`

**Resources.** 19899354 COMPLETED 00:02:37, 4 CPUs, 64 G, MaxRSS 2.19 GiB; 19899356 COMPLETED
00:02:13.

**Next.** M39 — malignant ECM-like cells versus true fibroblasts.
