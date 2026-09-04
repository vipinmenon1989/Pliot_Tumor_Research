# Phase 4 — Robustness Report

**Milestone M34** · SLURM 19897639

> **n = 4 patients.** No population-level or epidemiological claim follows from anything here.
> `sample_id` = patient = dataset, so biological and technical effects cannot be fully
> separated — the same confound recorded in Phase 2 and Phase 3.

---

## Q1 — Does each malignant state occur in multiple patients?

**No. This is the central negative result of Phase 4's state analysis.**

```text
recurrent (>=3 patients, each contributing >=5% and >=10 cells) :  0 of 8
shared between 2 patients                                       :  1
patient-specific or patient-dominated                           :  7
fraction of malignant cells in patient-private states           :  0.972
```

Only `Schwann_like` (179 cells; MPNST_3 120, MPNST_4 59; dominant share 0.67) is genuinely
shared, and it partly rests on the sample whose CNV partition failed the sanity gate. Every
other state has a dominant patient holding 83–100% of its cells.

Recurrence is judged by patient **share**, not patient presence: `Mesenchymal_ECM-2` is present
in all four patients but 1,587 of its 1,598 cells are MPNST_4, with 1, 3 and 7 cells from the
others. Calling that recurrent would violate §36, and an earlier version of this analysis did
exactly that before being corrected.

## Q2 — Are the SCEVAN CNV patterns patient-specific?

**Predominantly yes, with a small shared core.** Pairwise Jaccard of broad (≥10 Mb) clonal
event sets:

| pair | Jaccard |
| --- | ---: |
| MPNST_1 ↔ MPNST_2 | 0.233 |
| MPNST_2 ↔ MPNST_4 | 0.211 |
| MPNST_1 ↔ MPNST_4 | 0.200 |
| MPNST_1 ↔ MPNST_3 † | 0.077 |
| MPNST_3 ↔ MPNST_4 † | 0.077 |
| MPNST_2 ↔ MPNST_3 † | 0.071 |

† MPNST_3 failed the immune sanity gate; its low similarity to the others is consistent with
its "tumour" set being immune cells.

Among the three reliable patients only **20–23% of broad events are shared by any two**, so CNV
architecture is mostly individual. The shared core is small but consistent: **chr18 loss 3/3,
chr2 gain 3/3, chr7 gain 3/3**, with chr3/5/6/11/15/19/22 events in 2/3. Only **1 of 59** broad
events appears in all four patients.

Wording rule (§50): these are **shared or recurrent CNA patterns across patients**, never "the
same clone". Clones are patient-scoped by construction.

## Q3 — Are the major CCC changes caused by one patient?

**The `Newly-supported` class is almost entirely single-patient and must be discounted.**

| change class | n | median patients | supported in 1 patient only | in ≥3 patients | fraction single-patient |
| --- | ---: | ---: | ---: | ---: | ---: |
| **Newly-supported** | 5,421 | 1 | **4,977** | **27** | **0.918** |
| Weakened | 8,446 | 2 | 3,619 | 1,031 | 0.428 |
| Strengthened | 2,746 | 2 | 284 | 1,031 | 0.103 |
| Lost | 3,963 | 0 | 0 | 0 | — |
| Sender-reassigned | 3,051 | 0 | 0 | 0 | — |

**91.8% of the 5,421 "newly supported" interactions rest on a single patient, and only 27 reach
≥3 patients.** So the apparent gain from refinement is weak evidence and is reported as such
rather than promoted. By contrast `Strengthened` is only 10.3% single-patient with 1,031
interactions at ≥3 patients — a much better-supported class.

## Q4 — Do the Phase 3 core interactions remain leave-one-patient-out stable?

**Yes.** Retention of the 34,893 supported refined interactions when each patient is dropped:

| patient dropped | interactions retained | retention |
| --- | ---: | ---: |
| MPNST_1 | 25,382 | **0.727** |
| MPNST_2 | 28,780 | 0.825 |
| MPNST_3 | 31,973 | 0.916 |
| MPNST_4 | 33,315 | 0.955 |

All well above 0.5, so **no single patient dominates the refined interaction set**. MPNST_1
contributes most (dropping it costs 27%), which is consistent with it having the largest
evaluable population set. This range (0.727–0.955) is comparable to Phase 3's 0.71–0.79.

## Q5 — Does the refined malignancy conclusion depend on one patient?

**The headline malignant-fraction expansion does — and this is the most important caveat in
Phase 4.**

| patient dropped | cells | refined Malignant | refined fraction | Phase 2 fraction | fold change |
| --- | ---: | ---: | ---: | ---: | ---: |
| **none** | 19,716 | 6,434 | **0.326** | 0.173 | **1.88×** |
| MPNST_1 | 12,101 | 4,548 | 0.376 | 0.082 | 4.58× |
| MPNST_2 | 17,432 | 5,783 | 0.332 | 0.191 | 1.74× |
| MPNST_3 | 16,776 | 6,222 | 0.371 | 0.190 | 1.95× |
| **MPNST_4** | 12,839 | 2,749 | **0.214** | **0.215** | **1.00×** |

**Without MPNST_4 the refined malignant fraction equals the Phase 2 fraction (0.214 vs 0.215).**
Not because refinement did nothing there, but because two opposing corrections cancel: roughly
1,149 fibroblasts promoted in MPNST_1 and MPNST_2 against roughly 1,933 `MPNST-Tumor` cells
demoted to `Ambiguous`.

So the claim "**17.35% → 32.63%**" is **not patient-robust** and must not be quoted as a cohort
property. The claim "**Phase 2's fibroblast compartment is predominantly malignant**" is a
different claim and *is* patient-robust (Q6).

## Q6 — Is the fibroblast conclusion consistent across patients?

| Phase 2 population | patients | patients majority-malignant after refinement | min | max | spread |
| --- | :--: | :--: | ---: | ---: | ---: |
| **Fibroblast** | 4 | **3** | 0.00 | 0.98 | 0.98 |
| MPNST-Tumor | 4 | 2 | 0.07 | — | — |
| Candidate-Malignant-Unresolved | 4 | 1 | 0.00 | 0.70 | 0.70 |
| Pericyte-VSMC | 4 | 1 | 0.00 | — | — |
| Uncertain | 4 | 1 | 0.00 | — | — |

**`Fibroblast` is the only disputed population that is majority-malignant in a majority of
patients** (MPNST_1 0.62, MPNST_2 0.55, MPNST_4 0.98; MPNST_3 0.08, the failed sample). The
`Candidate-Malignant-Unresolved` conclusion rests on MPNST_1 alone (0.70; 1,196 of its 1,231
cells are MPNST_1), and the `Pericyte-VSMC` and `Uncertain` malignant fractions rest on MPNST_4
alone. Those are stated as single-patient findings.

## Summary — what is robust and what is not

| claim | robust? | basis |
| --- | :--: | --- |
| Phase 2's `Fibroblast` compartment is predominantly malignant | **yes** | 3/4 patients; threshold-independent (0.15–0.40); reference-independent (agreement 0.966–0.997); clone composition |
| The Phase 3 CCC architecture survives refinement | **yes** | High concordance −2.6%; median axis change −3.1%; no axis lost; LOSO 0.727–0.955 |
| ECM→integrin interactions are substantially tumour-autocrine | **yes** | all four axes sender-reassigned; 1,565 fibroblast interactions reassigned |
| A small shared CNA core exists (chr18 loss, chr2/chr7 gain) | qualified | 3/3 reliable patients, but n = 3 and Jaccard only 0.20–0.23 overall |
| Refined malignant fraction is 32.63% | **no** | equals the Phase 2 fraction once MPNST_4 is dropped |
| `Candidate-Malignant-Unresolved` is malignant | qualified | resolved (0 non-malignant) but 97% one patient |
| Malignant transcriptional states are recurrent | **no** | 0 of 8; 97.2% of malignant cells in patient-private states |
| A tumour-state → TME model | **no** | only 12 of 176 evidence rows reach ≥3 patients, all in one 290-cell state |
| Newly-supported interactions after refinement | **no** | 91.8% single-patient |

## Tables and figures

`MALIGNANCY_LEAVE_ONE_PATIENT_OUT.tsv` · `MALIGNANCY_PATIENT_CONSISTENCY.tsv` ·
`CNV_PATTERN_PATIENT_SIMILARITY.tsv` · `CCC_CHANGE_PATIENT_DEPENDENCE.tsv` ·
`CCC_REFINED_LEAVE_ONE_PATIENT_OUT.tsv`. Figures `34_01_state_patient_representation` ·
`34_02_malignancy_leave_one_patient_out` · `34_03_refinement_patient_consistency` ·
`34_04_ccc_leave_one_patient_out`.
