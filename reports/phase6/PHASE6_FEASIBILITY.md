# Phase 6 Feasibility — Clonal, Regulatory and Plasticity Architecture

*Milestone M42 · 2026-09-03 · SLURM 19899748 (COMPLETED, 8 CPUs, ReqMem 128 G, MaxRSS 21.92 GiB)*

## 1. Input integrity

`results/phase5/phase5_final_object.rds` was verified against its own manifest by **both**
checksums before loading — md5 `839e5157bc7c3470ddf86746c2e719e1`, sha256
`fe99ecf51154046145cf21f9c6960664d10205b90736bb13c61cccf409c40f00`. 19,716 cells, 196 metadata
columns. The four frozen Phase 4 refined counts (6,434 / 9,078 / 3,766 / 438) were re-asserted inside
the Phase 5 object and all held.

## 2. Required fields

All present: `tumor_clone_phase4`, `scevan_clone`, `scevan_sample_reliable`, the four `cnv_*`
metrics, `sample_id`, `malignancy_confidence`, `malignancy_refined`, `tumor_state_phase4`, the
eight `program_P*_score` fields and the five Phase 5 dominant-program fields.

## 3. MPNST_3 exclusion — re-asserted, not assumed

The exclusion was not taken on trust from the phase brief. The **frozen Phase 4
`scevan_sample_reliable` flag** was read back per patient and independently reproduces it:

| patient | scevan_sample_reliable | clones | malignant cells | malignant with a clone |
| --- | :--: | ---: | ---: | ---: |
| MPNST_1 | TRUE | 7 | 1,886 | 1,886 |
| MPNST_2 | TRUE | 4 | 651 | 651 |
| **MPNST_3** | **FALSE** | 3 | 212 | 91 |
| MPNST_4 | TRUE | 8 | 3,685 | 3,628 |

MPNST_3 failed the Phase 4 SCEVAN immune sanity gate: its primary and sensitivity runs are inverted
(agreement 0.0864 against 0.9663–0.9974), six of nine canonical immune populations are called ~100%
malignant, and its three inferred clones are T/NK, plasma/pDC/B and plasma-dominated. **Its clone
structure is excluded from every clone-based inference.** Its Phase 5 transcriptional programs are
still described.

## 4. Clone evaluability

19 clones in the three reliable patients. At the declared minimum of 20 cells, **18 are evaluable**;
at 50 cells 18; at 100 cells 16. `MPNST_4_clone8` holds 5 malignant cells and is reported
**NOT EVALUABLE** rather than dropped.

**178 malignant cells carry no clone label** — Phase 4 assigned a clone only to cells its SCEVAN run
placed. Those cells are excluded from clone analyses and the number is stated, not hidden.
6,165 of the 6,256 clone-labelled malignant cells are in the three reliable patients.

## 5. Working artefacts

`phase6_cell_metadata.rds` (19,716 × 196) and `phase6_malignant_lognorm.rds` (31,764 genes × 6,434
malignant cells, the frozen RNA `data` layer) were exported so the 5.6 GB object is opened once.
Phase 4 and Phase 5 md5s were both re-verified unchanged after the export.

## 6. Feasibility verdict

**Phase 6 is feasible**, with four constraints carried into the design:

1. **Clone inference rests on 3 patients**, not 4.
2. **Clone labels are patient-scoped.** `MPNST_1_clone1` and `MPNST_4_clone1` are unrelated names;
   every association is computed inside one patient and never pooled as if clones matched.
3. **MPNST_2's malignant cells are entirely dominated by one program** (P3, 651/651), so its
   clone-program analysis has no variation at the hard-assignment level and must be read on the
   continuous scores.
4. **The CNA→expression analysis cannot be orthogonal validation.** SCEVAN inferred those events
   from expression, so agreement is expected and is an internal consistency check.
