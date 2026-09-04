# M33 — Tumour-State → TME Model

**Status: COMPLETE** · 2026-09-03 · SLURM 19897639 (chained with M31 and M34)

## Outcome — a negative result, stated as the headline

**The tumour-state → TME model asked for in §47 and §68 cannot be established in this cohort,
and the §68 structure was therefore NOT imposed on the data.**

M31 established that 7 of 8 malignant states are patient-specific or patient-dominated and that
97.2% of malignant cells sit in such states. Because `sample_id` = patient = dataset, a state
confined to one patient cannot have its signalling separated from that patient's identity.

| state | cells | evidence rows | max evaluable patients | rows ≥3 patients |
| --- | ---: | ---: | :--: | ---: |
| Cycling | 290 | 22 | **3** | **12** |
| Schwann_like | 179 | 22 | **2** | 0 |
| Mesenchymal_ECM-1…5 | 5,893 | 110 | **1** each | 0 |
| Interferon | 72 | 22 | **1** | 0 |

Of 176 evidence rows: **12 reach ≥3 patients, 34 reach ≥2, 121 rest on one patient, 21 on none.**
**All 12 of the ≥3-patient rows belong to `Cycling`**, which is itself 83% one patient. `Cycling`
leads all six programmes purely because it is the only state present in enough patients to be
evaluated — an evaluability artefact, not a biological preference. Reporting it as "cycling
tumour cells drive myeloid signalling" would be wrong, and it is not reported that way.

**The five `Mesenchymal_ECM` states hold 91.3% of the refined malignant compartment and have
exactly one evaluable patient each.** The populations that matter most are the ones whose
communication is least assessable here.

## What was done

Deliberately **targeted** (§45): LR frameworks were **not** run for every tumour state, which
with 8 states — most containing one patient's cells — would have manufactured combinatorial
noise. Instead:

* per-state, per-patient expression of the Phase 3-prioritized ligands (**40 of 40 present**)
  and receptors (**39 present**), at the same 10% detection floor and 10-cell stratum minimum
  Phase 3 used; 440 sender rows across 11 evaluable state × patient strata, 1,950 receiver rows;
* joint expression counted only when ligand **and** receptor both pass the floor **in the same
  patient**;
* cross-referenced against the M32 axis verdicts;
* Phase 3 NicheNet and LochNESS results **reused verbatim**, never recomputed.

**No composite score was invented.** Every column of `TUMOR_STATE_TME_EVIDENCE.tsv` is an
independent evidence stream and none are averaged (§47).

## §48 answers

**A/B/C/D — which states engage macrophages, NK/CD8, the vasculature, the perivascular niche:
not determinable.** Both evaluable states (`Cycling`, `Schwann_like`) show joint expression for
every prioritized programme; the five dominant ECM states cannot be compared. No preference can
be established.

**E — are previously fibroblast-like cells actually malignant? Yes, predominantly.** 4,036 of
5,064, majority-malignant in 3 of 4 patients, threshold- and reference-independent. The firmest
conclusion in Phase 4.

**F — does correcting malignant identity change the major Phase 3 conclusions? No.** High
concordance fell 2.6%, median axis change −3.1%, no axis lost. The substantive change is
**attributional**: 3,051 sender reassignments concentrated in the ECM→integrin axis.

## §46 — the Phase 3 distinction preserved

NicheNet is **not** forced to support APP–CD74. Phase 3 established that predicted APP–CD74
engagement and the CSF1-explained macrophage state are distinct findings; the myeloid
ligand programme (APP, CD99, ANXA1, HLA-F, THBS1) and the receiver-state cytokine programme
(CSF1, IL15, TGFB1) are carried as **separate** programmes and are never merged. LochNESS
enters only as receiver-lineage context and remains negative.

## Outputs

`TUMOR_STATE_TME_EVIDENCE.tsv` (176 rows) · `TUMOR_STATE_PROGRAMME_LEADERS.tsv` ·
`TUMOR_STATE_LIGAND_EXPRESSION.tsv.gz` · `RECEIVER_RECEPTOR_EXPRESSION.tsv.gz` ·
figures `33_01_tumor_state_specific_ccc`, `33_02_state_ligand_expression_by_patient`,
`33_03_programme_carriage_by_state` · `reports/phase4/TUMOR_STATE_COMMUNICATION_REPORT.md`

## Next milestone

**M34 — robustness.**
