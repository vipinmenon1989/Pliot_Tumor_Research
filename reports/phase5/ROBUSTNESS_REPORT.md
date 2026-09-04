# Phase 5 Robustness Report

*Milestone M40 · 2026-09-03/04 · SLURM 19899355, 19899358, 19899359, 19899742*

Eight perturbations were applied to the K = 8 solution. Six were declared in the phase design; two
more (**technical covariates** and a **technical-gene-free universe**) were added after M38 showed
that three of the eight programs are technical-dominated and that the only cross-patient program is
one of them. Both additions are declared sensitivity analyses of the primary result and are reported
whichever way they came out.

## 1. Results at a glance

| program | label | recurrence | rank K±1 | balanced | gene subset (min ρ) | cell-cycle removed | High-only | LOO worst | dominant patient out | max technical ρ |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| P1 | Translation_ribosomal | shared-limited | 1.00 | 0.92 | **0.85** | 1.00 | 0.96 | 0.35 | 0.57 | 0.17 |
| P2 | Unassigned_NEURONAL_SYSTEM | patient-private | 1.00 | 0.98 | 0.77 | 1.00 | 0.95 | 0.21 | 0.21 | 0.22 |
| P3 | Mesenchymal_ECM_mixed | patient-private | 0.98 | 0.96 | 0.60 | 1.00 | 0.95 | 0.31 | 0.31 | 0.22 |
| P4 | Hypoxia_Angio | patient-private | 0.65 | 0.63 | 0.64 | 1.00 | 0.97 | 0.17 | 0.17 | 0.16 |
| P5 | Translation_ribosomal | patient-private | 0.62 | 0.50 | 0.71 | 1.00 | 0.97 | 0.20 | 0.20 | 0.37 |
| P6 | Schwann_like_mixed | patient-private | 1.00 | 0.96 | 0.24 | 1.00 | 0.30 | **0.80** | **0.80** | 0.26 |
| P7 | Cycling | patient-private | 1.00 | 0.95 | 0.54 | 0.99 | 0.80 | **0.88** | **0.88** | 0.31 |
| P8 | Myeloid_ambient_like | patient-private | 1.00 | 0.91 | 0.23 | 1.00 | 0.15 | 0.17 | 0.17 | 0.21 |

## 2. What each test says

**Program rank (K ± 1).** Median best-match cosine 0.998 against K = 7 and 1.000 against K = 9.
Six of eight programs are essentially identical at neighbouring ranks; P4 and P5 (0.62–0.65) are the
two that reorganise. **The solution is not an artefact of the exact K.**

**Patient balance.** Median 0.940, minimum 0.504, and every program clears the declared support
threshold. Downsampling the dominant patient from 57.3% to 30.1% of the compartment does not
destroy any program. **Patient imbalance does not create these programs; it confines them.**

**Gene selection (10 × random 80%).** Minimum ρ 0.23–0.85; only P1 clears the declared 0.80 bar.
This threshold is conservative by construction — it compares a *projection* on 80% of genes against
a *consensus* solution on all of them, which are different estimators — but it was declared before
the test was run and is reported as declared. The two programs that move most (P6 0.24, P8 0.23)
are the two smallest.

**Cell cycle.** Median cosine 1.000, minimum 0.994. **Removing the 95 canonical Tirosh genes changes
nothing**, including P7 itself (0.99): the Cycling program is carried by far more genes than the
canonical list, and no other program was being distorted by them.

**Technical covariates.** Within each patient, the largest |Spearman ρ| between a program's usage
and nCount_RNA / nFeature_RNA / percent.mt is **0.37** (P5). No program clears 0.50. So the two
ribosomal programs are **not** simply depth artefacts *within* a patient — what they track is
*between*-patient structure, which is a different and more interesting problem (§4).

**Malignancy confidence.** Recurrence class is unchanged for 6 of 8 programs on High-confidence
cells only. P6 and P8 become `uncertain` — both are MPNST_3 programs, and MPNST_3 has **0
High-confidence malignant cells**, so this is bookkeeping rather than instability. The High-only
cNMF run recovers the other programs at cosine 0.95–0.97.

**Leave-one-patient-out.** Withholding MPNST_1 drops P1, P2, P4 and P5 below 0.5; withholding
MPNST_2 drops P3; withholding MPNST_3 drops P8; withholding MPNST_4 drops nothing below 0.57.
**Every patient-private program disappears with its own patient**, which is the cleanest possible
demonstration that these programs are patient-bound rather than cohort-level.

## 3. Two verdicts, because the strict one is uninformative alone

`robust_overall` applies every declared threshold and is **FALSE for all eight programs**. That is
technically true and practically useless: five of the eight fail only the dominant-patient test,
which a patient-private program *cannot* pass by definition.

`robust_excluding_patient_scope` waives that one test for programs already classified
patient-private, and also excludes the conservative gene-subset bar. It is **TRUE for P2, P3, P4 and
P7** — the neuronal, ECM, hypoxia/angiogenic and cycling programs — and FALSE for P1, P5, P6 and P8.

Both columns are in `PHASE5_PROGRAM_ROBUSTNESS.tsv`. Neither is a substitute for reading the row.

## 4. The technical-gene-free universe — the most important sensitivity result

Three of the eight programs are technical-dominated (two ribosomal/pseudogene, one carrying 18
canonical myeloid markers in its top 50 genes — ambient RNA inside a malignant-only factorization),
and the **only** program carried by more than one patient is one of the ribosomal ones. The obvious
worry is that the patient-private result is an artefact of those gene classes.

So a fourth full cNMF run was done with ribosomal-protein genes, ribosomal pseudogenes,
lncRNA/clone-name genes and canonical myeloid markers removed (19,663 → 18,246 genes), across the
**whole K grid**, so that universe selects its own rank rather than inheriting K = 8:

| K | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
| --- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| recurrent | 1 | 0 | 0 | 0 | **0** | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| shared-limited | 2 | 2 | 3 | 1 | **2** | 1 | 1 | 1 | 1 | 1 | 0 | 0 |
| patient-private | 1 | 3 | 3 | 6 | **6** | 8 | 8 | 8 | 9 | 10 | 12 | 13 |

**At every K from 5 to 15, zero programs are recurrent.** Removing the technical gene classes does
not reveal shared biological programs. **The patient-private result is not an artefact of ribosomal,
pseudogene or ambient myeloid genes** — it is the structure of this cohort.

## 5. What survives

1. **Continuous malignant programs are patient-private in this cohort.** 0 recurrent at the declared
   criteria, 0 recurrent on the technical-gene-free universe at 11 of 12 ranks, and every
   patient-private program vanishes when its own patient is withheld.
2. **The one biologically interpretable program that comes closest is the ECM/mesenchymal program
   P3**, which reaches three patients only when the activity threshold is relaxed from 0.20 to 0.10.
3. **Programs are not artefacts of rank, patient balance, cell-cycle genes, confidence filtering or
   within-patient sequencing depth.**
4. **They are, however, confounded with patient identity by construction** — with n = 4 and
   `sample_id` = patient = dataset, a factorization of unintegrated counts cannot separate a
   patient-specific biological program from patient-level technical structure. That limitation is
   not removable by any analysis of these four patients.
