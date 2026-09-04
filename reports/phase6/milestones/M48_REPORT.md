# M48 — Integrated Malignant Architecture and Per-Patient Portraits

**Date** 2026-09-04 · **Status** COMPLETE · **SLURM** 19899749

**Design.** Three candidate architectures were specified with numeric criteria **before** the
evidence was assembled, and whichever the data support is reported. The conclusion was not chosen in
advance.

```text
A clone-constrained   >= 70% of evaluable clones concentrated (>= 80% of cells in one program)
                      AND median eta^2 >= 0.30
                      AND between-clone divergence >= 50% of within-clone dispersion
B program-diverse     >= 70% of evaluable clones span >= 2 programs
                      AND <= 25% of program x patient pairs clone-associated
C mixed               otherwise
```

**Evidence.**

| quantity | value |
| --- | --- |
| evaluable clones (≥ 20 cells) | 18 |
| spanning ≥ 2 programs | **4 (22%)** |
| concentrated in one program | 11 (61%) |
| median η² | **0.059** → ~94% of program variance is WITHIN clones |
| max η² | 0.508 (P5 in MPNST_1) |
| clone-associated program × patient pairs | 7 of 24 (29%) |
| between-clone divergence / within-clone dispersion | **0.034 / 0.0014 / 0.0069** (median 0.0069) |

**Model A fails** on all three criteria — most decisively on the last: between-clone divergence is
**0.7% of within-clone dispersion**, so a patient's CNA-defined clones are transcriptionally
near-interchangeable. **Model B fails** because only 22% of clones span ≥ 2 programs, and all four
that do are in one patient.

## SELECTED: Model C — mixed architecture

> A minority of transcriptional programs associate detectably with specific CNA-defined clones while
> the great majority of program variance sits **within** clones. Median η² is 0.059, so roughly 94%
> of each program's variance is within-clone, and between-clone divergence is only 0.7% of
> within-clone dispersion — a patient's CNA-defined clones are transcriptionally
> near-interchangeable. How much program diversity a clone contains is itself patient-specific, so
> this is a mixed architecture rather than a single rule holding across the cohort.

## Per-patient portraits

| | MPNST_1 | MPNST_2 | MPNST_4 | MPNST_3 |
| --- | --- | --- | --- | --- |
| clone structure usable | yes | yes | yes | **NO** |
| malignant cells | 1,886 | 651 | 3,685 | 212 |
| clones / evaluable | 7 / 7 | 4 / 4 | 8 / 7 | 3 / — |
| dominant program | P2 Neuronal | P3 **Mesenchymal_ECM** | P1 Translation_ribosomal | P6 Schwann_like |
| effective programs per clone | **2.15** | 1.00 | 1.19 | — |
| median η² | 0.072 | 0.050 | 0.072 | — |
| clone-associated programs | 3 | 2 | 2 | — |
| broad CNA events evaluable | — | — | — | excluded |

**MPNST_3**: its Phase 5 transcriptional programs are described (dominant program P6, Schwann-like),
but **no clone-based interpretation is offered** and none was invented.

**Outputs.** `PATIENT_TUMOR_ARCHITECTURE.tsv` · `PHASE6_INTEGRATED_EVIDENCE.tsv`

**Next.** M49 — robustness.
