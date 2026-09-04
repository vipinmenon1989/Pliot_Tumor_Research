# Phase 6 Robustness Report

*Milestone M49 · 2026-09-04 · SLURM 19899749*

Five perturbations. Leave-one-patient-out is deliberately **not** applied to the clone analyses:
they are inherently within-patient, so withholding a patient deletes the analysis rather than
testing it. Patient-level direction consistency is tested directly instead.

## 1. Clone-size threshold

| minimum clone size | clones used | NOT EVALUABLE | median η² |
| ---: | ---: | ---: | ---: |
| ≥ 20 cells | 18 | 1 | **0.0586** |
| ≥ 50 cells | 18 | 1 | **0.0586** |
| ≥ 100 cells | 16 | 3 | **0.0586** |

**The clone-program coupling estimate is invariant to the threshold.** Clones failing a threshold
are reported NOT EVALUABLE, never hidden.

## 2. Dominant-program margin

| margin | clones | median effective programs | clones spanning ≥ 2 | median Mixed cells |
| ---: | ---: | ---: | ---: | ---: |
| ≥ 0.00 | 18 | 1.222 | 22% | 0% |
| ≥ 0.05 | 18 | 1.211 | 22% | ~1% |
| **≥ 0.10** | 18 | **1.209** | **22%** | ~1% |
| ≥ 0.15 | 18 | 1.190 | 22% | ~2% |
| ≥ 0.20 | 18 | 1.185 | 22% | ~2% |

The diversity estimate moves by 3% across a fourfold change in the hard-assignment threshold, and
the fraction of multi-program clones does not move at all. **The within-clone diversity result is not
a threshold artefact.**

## 3. High-confidence malignant cells only

3,261 of 6,165 clone-labelled cells in the reliable patients are High-confidence (MPNST_1 1,861;
MPNST_4 936; MPNST_2 464).

| quantity | agreement between all-malignant and High-only |
| --- | ---: |
| clone-program η² (24 patient × program pairs) | Spearman **0.901** |
| TF activity ↔ program ρ (5,520 pairs) | Spearman **0.940** |
| PROGENy ↔ program ρ (112 pairs) | Spearman **0.962** |

Per-pair η² changes are in the third decimal (largest Δ = 0.022, P3 in MPNST_2). **Every major
conclusion survives restriction to High-confidence cells.**

## 4. Patient-level direction consistency

| layer | pairs | strong pooled (|ρ| ≥ 0.20) | also concordant in ≥ 3 patients | fraction |
| --- | ---: | ---: | ---: | ---: |
| TF (CollecTRI) | 5,520 | 2,368 | 1,922 | **0.812** |
| Pathway (PROGENy) | 112 | 48 | 39 | **0.813** |
| Pathway (Hallmark) | 400 | 193 | 148 | **0.767** |

**Read this with the right prior.** With four patients, three agreeing in sign happens ≈ 31% of the
time by chance, so "concordant in ≥ 3 patients" is a weak bar and 77–81% is not as impressive as it
looks. It is reported because the alternative — accepting a strong pooled association on its own —
is worse: with `sample_id` = patient = dataset, a pooled cell-level correlation can be one patient's
signal entirely. **The named regulators in M45/M46 rest on their direction and magnitude, not on
this count.**

## 5. What survives

1. **Clone-program coupling is weak everywhere and invariant to how clones are filtered.** Median
   η² 0.059 at every clone-size threshold; ~94% of program variance is within-clone.
2. **Between-clone divergence is 0.1–3.4% of within-clone dispersion.** A patient's clones are
   transcriptionally near-interchangeable. This is the finding that rules out Model A.
3. **Within-clone diversity is real in MPNST_1 and essentially absent in MPNST_2 and MPNST_4**, and
   that difference is not a threshold effect.
4. **The regulatory and pathway architecture is stable to confidence filtering** (ρ 0.94 / 0.96).
5. **One confound is not dismissed.** Program dispersion tracks the fraction of High-confidence
   cells in a clone at ρ = 0.72 (M44). The effective-number metrics stay below the 0.7 bar
   (0.56–0.62), so the multi-program finding does not rest on the confounded metric — but any
   plasticity reading of *dispersion* specifically must carry this number.
