# M43 — Clone ↔ Program Coupling

**Date** 2026-09-03 · **Status** COMPLETE · **SLURM** 19899748

**Question.** Within each clone-reliable patient, how much of a malignant program's variation sits
*between* CNA-defined clones rather than *within* them?

**Design.** Every association is computed **inside one patient**. `MPNST_1_clone1` and
`MPNST_4_clone1` are unrelated labels and are never pooled. Strength is η² — the between-clone share
of a program's variance — with a null built by shuffling clone labels **inside** the patient
(1,000 permutations), which preserves clone sizes and the score distribution. Clones below 20 cells
are printed NOT EVALUABLE (`MPNST_4_clone8`, 5 cells).

**Result.**

| patient | clones used | cells | median η² | max η² | program at max | clone-associated programs |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| MPNST_1 | 7 | 1,886 | 0.072 | 0.508 | P5 | 3 of 8 |
| MPNST_2 | 4 | 651 | 0.050 | 0.154 | P2 | 2 of 8 |
| MPNST_4 | 7 | 3,623 | 0.072 | 0.208 | P2 | 2 of 8 |

**7 of 24 program × patient pairs are clone-associated** (η² ≥ 0.10 and permutation p < 0.01).
21 of 24 have permutation p < 0.001 — with 651–3,623 cells almost any non-zero η² is "significant",
which is exactly why the effect size and not the p-value carries the conclusion here.

**The number that matters: median η² ≈ 0.06, so roughly 94% of each program's variance sits WITHIN
clones.** Even the largest value, 0.51 for P5 in MPNST_1, leaves half the variance within clones.

**One program appears across essentially all clones of its patient.** MPNST_1's P2 is active in 7 of
7 clones (median usage 0.26–0.81), MPNST_2's P3 in 4 of 4 (0.84–0.93), MPNST_4's P1 in 7 of 7
(0.85–0.89). A patient's dominant program is not a property of one clone.

**Outputs.** `CLONE_PROGRAM_ACTIVITY.tsv` (152 rows) · `CLONE_PROGRAM_ASSOCIATION.tsv` (24) ·
`M43_CLONE_DOMINANT_PROGRAM.tsv` · `M43_ASSOCIATION_SUMMARY_BY_PATIENT.tsv` ·
`M43_PROGRAM_SPREAD_ACROSS_CLONES.tsv`

**Next.** M44 — within-clone program diversity.
