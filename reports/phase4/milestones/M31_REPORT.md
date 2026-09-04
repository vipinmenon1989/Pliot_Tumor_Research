# M31 — Malignant-Only Tumour-State Analysis

**Status: COMPLETE** · 2026-09-03 · SLURM **19897639** COMPLETED 00:02:22, MaxRSS 5.41 GiB
(first run **19897617** COMPLETED 00:02:24 and was superseded — its recurrence labelling
violated §36; see *Errors*)

## Input

`malignancy_refined == "Malignant"` — **6,434 cells**, not the broad Phase 2 `MPNST-Tumor`
label (§30). Expression from the M28-extracted RNA counts, LogNormalize 1e4; the 6 GB Phase 2
object was not reloaded and the four-model SCT assay was not used.

| Confidence | cells | | Patient | cells | | Phase 2 origin | cells |
| --- | ---: | --- | --- | ---: | --- | --- | ---: |
| High | 3,261 | | MPNST_1 | 1,886 | | Fibroblast | 4,036 |
| Moderate | 3,173 | | MPNST_2 | 651 | | MPNST-Tumor | 1,405 |
| | | | MPNST_3 | 212 | | Candidate-Malignant-Unresolved | 836 |
| | | | MPNST_4 | 3,685 | | Uncertain / Pericyte-VSMC | 157 |

## New reductions — Phase 2's are untouched (§31)

`malignant_pca` → `malignant_harmony` (batch variable `sample_id`, the same variable Phase 2
used) → `malignant_umap`. Dims 1:30, capturing **93.0%** of variance, matching the Phase 2
convention.

The frozen `pca`, `postint_harmony` and `postint_umap_harmony` were **never modified** —
asserted, not assumed: M35 confirms every Phase 2 embedding is numerically identical in the
final object.

**Was a new reduction justified?** Yes, and it was measured rather than asserted: kNN(k=20)
overlap between the Phase 2 Harmony embedding and the malignant-only Harmony embedding is
**0.1733**. The Phase 2 embedding was optimised to separate cell *types*; it does not resolve
states *within* the malignant compartment.

## Resolution sweep — no predetermined number of states (§32)

| resolution | clusters | <30 cells | multi-patient | frac multi-patient |
| ---: | ---: | ---: | ---: | ---: |
| **0.2** | **8** | **0** | 4 | **0.500** |
| 0.3 | 10 | 1 | 4 | 0.400 |
| 0.4 | 10 | 0 | 5 | 0.500 |
| 0.5 | 11 | 0 | 5 | 0.455 |
| 0.7 | 15 | 0 | 7 | 0.467 |
| 1.0 | 17 | 0 | 7 | 0.412 |

Selection rule, printed before it was applied: the largest resolution with no cluster under 30
cells *and* every cluster multi-patient; failing that, the largest fraction of multi-patient
clusters, then the smaller resolution. **No resolution reached full multi-patient coverage**,
so the fallback selected **0.2** (tied on 0.500 with 0.4, smaller resolution wins).

## The eight states

| state | cells | frac | patients | dominant patient share | clones | median CNA burden | defining programme evidence |
| --- | ---: | ---: | :--: | ---: | ---: | ---: | --- |
| Mesenchymal_ECM-1 | 1,851 | 0.287 | 1 | 1.00 | 7 | 0.357 | ECM 0.619, next Stress 0.163 |
| Mesenchymal_ECM-2 | 1,598 | 0.248 | 4 | **0.99** | 13 | 0.250 | ECM 0.956, next Interferon 0.044 |
| Mesenchymal_ECM-3 | 1,291 | 0.200 | 1 | 1.00 | 8 | 0.269 | ECM 1.095, next Immune-interacting 0.156 |
| Mesenchymal_ECM-4 | 648 | 0.100 | 3 | **1.00** | 5 | 0.310 | ECM 1.255, next Immune-interacting 0.659 |
| Mesenchymal_ECM-5 | 505 | 0.078 | 1 | 1.00 | 7 | 0.265 | ECM 0.887, next Interferon 0.102 |
| Cycling | 290 | 0.045 | 4 | 0.83 | 19 | 0.293 | Cycling 0.960, next ECM 0.698 |
| **Schwann_like** | 179 | 0.027 | 2 | **0.67** | 3 | 0.279 | Schwann 1.222, next Immune-interacting 0.966 |
| Interferon | 72 | 0.011 | 1 | 1.00 | 2 | 0.154 | Interferon 1.140, next Immune-interacting 0.619 |

State names were **not imposed**: a cluster takes its top programme's name only if that
programme's median score is positive and exceeds the runner-up by ≥0.02, otherwise it is
`Uncertain` (§34). No cluster required `Uncertain`.

## The headline negative result (§36)

```text
recurrent (>=3 patients, each contributing >=5% and >=10 cells) :  0 of 8
shared between 2 patients                                       :  1
patient-specific or patient-dominated                           :  7
fraction of malignant cells in patient-private states           :  0.972
```

**No malignant transcriptional state is recurrent across patients.** Only `Schwann_like`
(179 cells; MPNST_3 120, MPNST_4 59) is genuinely shared, and even that leans on the sample
whose CNV partition failed the sanity gate. **97.2% of malignant cells sit in states dominated
by a single patient.**

This is expected rather than surprising — `sample_id` = patient = dataset, n = 4, and the
malignant compartment ranges from 212 to 3,685 cells per patient, so residual patient structure
dominates. It is nonetheless a negative result about *tumour-state generality* and is reported
as one. Nothing here supports a claim about MPNST tumour states as a disease.

## The coherent biological finding underneath it

The state structure maps cleanly onto the malignancy refinement, and the mapping explains
Phase 2's error:

```text
Phase 2 MPNST-Tumor cells        ->  Cycling (248) · Schwann_like (175) · Interferon (72)
Phase 4 promoted fibroblasts     ->  Mesenchymal_ECM-1..5 (1,245 + 1,128 + 637 + 512 + 481)
```

**Five of eight states are mesenchymal/ECM-programme states, and they are built almost entirely
from cells Phase 2 called `Fibroblast`.** Phase 2 detected the marker-legible malignant cells —
Schwann-like, cycling, interferon-responsive — and missed the Mes-NC-like ECM cells precisely
because an ECM transcriptional programme is what a fibroblast looks like. That is why
copy-number evidence, and not a better marker panel, was the right instrument.

## Clone versus state (§33)

Kept as distinct concepts. `CLONE_VS_TUMOR_STATE.tsv` cross-tabulates them; the descriptive
Cramér's V is reported in the figure caption as an **association, not an equivalence**. Clone
counts per state range from 2 to 19, so a single transcriptional state routinely spans many CNV
clones — one-to-one correspondence is absent.

## Figures (§35)

`31_01_tumor_state_umap` · `31_02_tumor_state_patient_contribution` · `31_03_scevan_clone_umap`
· `31_04_clone_vs_tumor_state` · `31_05_tumor_state_marker_dotplot` ·
`31_06_tumor_state_program_heatmap` · `31_07_tumor_state_marker_heatmap` ·
`31_08_tumor_state_patient_composition`. PDF + PNG each.

## Tables

`TUMOR_STATE_ASSIGNMENTS.tsv` · `TUMOR_STATE_MARKERS.tsv` ·
`TUMOR_STATE_PATIENT_DISTRIBUTION.tsv` · `TUMOR_STATE_DEFINITIONS.tsv` ·
`TUMOR_STATE_PROGRAM_SCORES.tsv` · `TUMOR_STATE_RESOLUTION_SWEEP.tsv` ·
`CLONE_VS_TUMOR_STATE.tsv`. Object: `results/phase4/tumor_states/malignant_only_object.rds`.

## Errors

**One substantive error of mine, caught and corrected before it propagated.** The first run
(19897617) classified recurrence on `patients_represented` alone, which labelled
`Mesenchymal_ECM-2` "recurrent across ≥3 patients" when 1,587 of its 1,598 cells were MPNST_4
and the other three patients contributed 1, 3 and 7 cells. §36 forbids exactly that. The rule
was rewritten to require that no patient exceeds 80% **and** that ≥3 patients each contribute
≥5% and ≥10 cells, and the script now prints the negative result explicitly instead of leaving
it to be inferred from a table. M33 was halted mid-run and both were re-executed.

## Scientific decisions

1. Subset on `malignancy_refined`, never the broad Phase 2 label.
2. New reductions created; Phase 2's frozen ones untouched, and that is asserted in M35.
3. New reduction justified by a measured kNN overlap of 0.173, not by assertion.
4. Number of states not predetermined; selection rule printed before application.
5. Programme names earned by a stated margin, with `Uncertain` available.
6. Recurrence judged by patient *share*, not patient presence.
7. Clone and state kept separate throughout.

## Next milestone

**M33 — tumour-state → TME evidence model**, then **M34 — robustness**.
