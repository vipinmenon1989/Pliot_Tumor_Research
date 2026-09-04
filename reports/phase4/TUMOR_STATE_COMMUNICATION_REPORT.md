# Phase 4 — Tumour-State → TME Communication Report

**Milestone M33** · integrated as an **evidence matrix, not a composite score** (§47)

> scRNA-seq ligand–receptor analysis infers communication **potential**. It does not establish
> physical adjacency or direct signalling. Joint expression means the ligand *and* the receptor
> each pass a 10% detection floor **in the same patient**.

---

## 1. The headline: this model cannot be built in this cohort

**§48 asks which tumour states preferentially engage macrophages, NK/CD8 cells, the vasculature
and the perivascular niche. The honest answer is that this dataset cannot attribute a
signalling programme to a tumour state rather than to a patient.**

The reason is structural, and it comes straight from M31: **7 of 8 malignant states are
patient-specific or patient-dominated, and 97.2% of malignant cells sit in such states.** A
state that exists in one patient cannot have its signalling separated from that patient's
identity, because `sample_id` = patient = dataset.

Evaluability, per state (a state × patient stratum needs ≥10 cells to be assessed at all):

| state | cells | evidence rows | max evaluable patients | rows supported in ≥3 patients | in ≥2 |
| --- | ---: | ---: | :--: | ---: | ---: |
| Cycling | 290 | 22 | **3** | **12** | 22 |
| Schwann_like | 179 | 22 | **2** | 0 | 12 |
| Mesenchymal_ECM-1 | 1,851 | 22 | **1** | 0 | 0 |
| Mesenchymal_ECM-2 | 1,598 | 22 | **1** | 0 | 0 |
| Mesenchymal_ECM-3 | 1,291 | 22 | **1** | 0 | 0 |
| Mesenchymal_ECM-4 | 648 | 22 | **1** | 0 | 0 |
| Mesenchymal_ECM-5 | 505 | 22 | **1** | 0 | 0 |
| Interferon | 72 | 22 | **1** | 0 | 0 |

Of 176 evidence rows: **12 reach ≥3 patients, 34 reach ≥2, 121 rest on one patient, 21 on none.**

**Every one of the 12 rows supported in ≥3 patients belongs to `Cycling`** — a 290-cell state
that is itself 83% one patient. `Cycling` leads all six programmes not because cycling tumour
cells are the dominant senders, but because it is the only state present in enough patients to
be evaluated. That is an artefact of evaluability, and reading it as biology would be wrong.

**The five `Mesenchymal_ECM` states hold 91.3% of the refined malignant compartment and have
exactly one evaluable patient each.** The dominant malignant populations in this cohort are the
ones whose communication cannot be assessed across patients.

## 2. What can still be said

Within the two evaluable states, joint ligand–receptor expression is present for every
prioritized programme:

| programme | Cycling (≥3 patients) | Schwann_like (≥2 patients) |
| --- | --- | --- |
| Myeloid antigen-presentation / inhibitory | Macrophage, Monocyte, Dendritic (3/3 each); pDC, B-cell (2/2) | Macrophage, Monocyte, Dendritic (2/2) |
| NK / CD8 inhibitory and activating | CD4-T, T-cell-other (3/3); NK, CD8-T (2/2) | CD4-T, T-cell-other (2/2) |
| Angiogenic / vascular | Endothelial (3/3); Pericyte-VSMC (2/2) | Endothelial (2/2) |
| Notch perivascular | Endothelial (3/3) | — |
| Growth factor / guidance / ECM | Endothelial (3/3); Fibroblast, Pericyte-VSMC (2/2) | Endothelial (2/2) |
| Receiver-state cytokines | Macrophage, Monocyte, CD4-T, Endothelial (3/3 each) | — |

This is consistent with Phase 3 and with M32 — the tumour compartment does carry
myeloid-directed, vascular, Notch and immune-modulatory ligand programmes — but it locates
those programmes in **the malignant compartment as a whole**, which is what Phase 3 already
established and M32 confirmed survives refinement. It does **not** license a
state-to-programme assignment.

## 3. Answers to §48, stated as the evidence permits

**A. Which tumour states preferentially engage macrophages?**
Not determinable. `Cycling` and `Schwann_like` both show joint APP/CD99/ANXA1/HLA-F ↔
CD74/PILRA/FPR1/LILRB expression with Macrophage, Monocyte and Dendritic cells, in 3 and 2
patients respectively. The five ECM states cannot be compared. No preference can be established.

**B. Which states preferentially engage NK/CD8 populations?**
Not determinable, and the evaluable evidence is weaker than for myeloid receivers: NK and CD8-T
reach only 2 patients even for `Cycling`. Consistent with Phase 3's finding that the lymphoid
picture is mixed rather than uniformly suppressive.

**C. Which states drive vascular signalling?**
Not determinable by state. `Cycling` shows VEGFA-axis joint expression with Endothelial in 3/3
evaluable patients and `Schwann_like` in 2/2. M32 separately shows VEGFA→KDR *strengthened*
(+11.1%) and VEGFA→FLT1 *weakened* (−28.0%) after refinement, so the axis is real but its
sender composition changed.

**D. Which states participate in the JAG/NOTCH perivascular circuit?**
Not determinable by state. `Cycling` shows Notch joint expression with Endothelial in 3/3.
M32 shows the circuit itself is stable to strengthened (JAG1→NOTCH3 0.0%, JAG1→NOTCH2 +7.5%,
JAG2→NOTCH2 0.0%) with **JAG1→NOTCH3 sender-reassigned**, so refinement changed *who* sends
without dissolving the circuit.

**E. Are previously fibroblast-like cells actually malignant?**
**Yes, predominantly** — 4,036 of 5,064, majority-malignant in 3 of 4 patients, threshold- and
reference-independent, with all four MPNST_2 clones and seven of eight MPNST_4 clones
fibroblast-dominated. Answered in `AMBIGUOUS_POPULATION_MALIGNANCY_REPORT.md`, and it is the
firmest conclusion in Phase 4.

**F. Does correcting malignant identity change the major Phase 3 conclusions?**
**No — the architecture survives.** High-concordance interactions fell 2.6% (5,848 → 5,698)
despite 4,036 cells entering and 2,015 leaving the tumour compartment; the median change across
the 23 named axes is −3.1%; **no axis was lost**. The substantive change is **attributional**:
3,051 interactions were sender-reassigned, concentrated in the ECM→integrin axis, where all
four collagen/FN1 axes changed sender.

## 4. Receiver-program integration (§46) — the Phase 3 distinction preserved

Phase 3's NicheNet and LochNESS results are **reused verbatim, not recomputed**, and NicheNet
is **not** forced to support APP–CD74. Phase 3 established that predicted APP–CD74 engagement
and the CSF1-explained macrophage state are **distinct findings**, and Phase 4 does not collapse
them: the myeloid programme's ligand set (APP, CD99, ANXA1, HLA-F, THBS1) and the
receiver-state cytokine set (CSF1, IL15, TGFB1) are carried as **separate programmes** in the
evidence matrix, and both are present in the evaluable states without being merged.

LochNESS enters only as receiver-lineage context and remains a Phase 3 **negative** result: no
receiver-state structure associated with the tumour-derived APP context in any lineage (best
descriptive p = 0.333 against a 0.167 floor at n = 4), with the score itself unstable. High
LochNESS would not mean cell–cell communication, and its absence does not refute the LR
findings.

## 5. Design honesty

Every LR framework was **not** run for every tumour state (§45) — with 8 states that would have
manufactured combinatorial noise on strata that mostly contain one patient's cells. Instead the
evidence is per-state, per-patient ligand and receptor expression on the Phase 3-prioritized
axes (40/40 prioritized ligands and 39 receptors present), cross-referenced against the M32
sensitivity verdicts and the reused Phase 3 receiver-response results. **No composite score was
invented.** Each column of `TUMOR_STATE_TME_EVIDENCE.tsv` is an independent stream and none are
averaged.

## 6. What would make this answerable

More patients. The limitation is not depth, coverage or method choice — it is that with n = 4
and `sample_id` = patient = dataset, a tumour state confined to one patient is
indistinguishable from that patient. A cohort in which the mesenchymal/ECM states recur across
individuals would make the state→TME question tractable; this one does not.

## 7. Tables and figures

`TUMOR_STATE_TME_EVIDENCE.tsv` (176 rows) · `TUMOR_STATE_PROGRAMME_LEADERS.tsv` ·
`TUMOR_STATE_LIGAND_EXPRESSION.tsv.gz` · `RECEIVER_RECEPTOR_EXPRESSION.tsv.gz`.
Figures `33_01_tumor_state_specific_ccc` · `33_02_state_ligand_expression_by_patient` ·
`33_03_programme_carriage_by_state`.
