# M34 — Robustness

**Status: COMPLETE** · 2026-09-03 · SLURM 19897639 (chained with M31 and M33)

Full analysis in `reports/phase4/ROBUSTNESS_REPORT.md`. This is the milestone summary.

## The four §49 questions

| Question | Answer |
| --- | --- |
| Does each malignant state occur in multiple patients? | **No.** 0 of 8 recurrent; 1 shared between 2 patients; 97.2% of malignant cells in patient-private states. |
| Are the SCEVAN CNV patterns patient-specific? | **Predominantly yes.** Pairwise Jaccard 0.20–0.23 among reliable patients; only 1 of 59 broad events in all four. Small shared core: chr18 loss, chr2 gain, chr7 gain, each 3/3. |
| Are the major CCC changes caused by one patient? | **The `Newly-supported` class is: 91.8% single-patient**, only 27 of 5,421 at ≥3 patients. `Strengthened` is far better supported (10.3% single-patient, 1,031 at ≥3). |
| Do the Phase 3 core interactions remain leave-one-patient-out stable? | **Yes.** Retention 0.727–0.955; no patient dominates. |

## Two additional questions worth asking

**Q5 — does the refined malignant fraction depend on one patient? Yes, and this is Phase 4's
most important caveat.** Dropping MPNST_4 gives a refined fraction of **0.214 against a Phase 2
fraction of 0.215** — a fold change of 1.00. The headline 17.35% → 32.63% is **not
patient-robust**, because outside MPNST_4 the promotions (~1,149 fibroblasts) and the demotions
(~1,933 `MPNST-Tumor` → `Ambiguous`) roughly cancel.

**Q6 — is the fibroblast conclusion consistent across patients? Yes.** `Fibroblast` is the only
disputed population majority-malignant in a majority of patients (0.62, 0.55, 0.98; 0.08 in the
sanity-failed sample). By contrast `Candidate-Malignant-Unresolved` rests on MPNST_1 alone and
the `Pericyte-VSMC` / `Uncertain` malignant fractions on MPNST_4 alone.

## What is and is not robust

**Robust**: the fibroblast compartment is predominantly malignant · the Phase 3 CCC
architecture survives refinement · ECM→integrin interactions are substantially tumour-autocrine.

**Not robust**: the 32.63% refined malignant fraction · recurrence of malignant transcriptional
states (0 of 8) · a tumour-state → TME model · newly-supported interactions.

**Qualified**: the shared CNA core (3/3 reliable patients, but n = 3 and overall Jaccard only
0.20–0.23) · `Candidate-Malignant-Unresolved` malignancy (resolved, but 97% one patient).

## Outputs

`MALIGNANCY_LEAVE_ONE_PATIENT_OUT.tsv` · `MALIGNANCY_PATIENT_CONSISTENCY.tsv` ·
`CNV_PATTERN_PATIENT_SIMILARITY.tsv` · `CCC_CHANGE_PATIENT_DEPENDENCE.tsv` ·
`CCC_REFINED_LEAVE_ONE_PATIENT_OUT.tsv` · figures `34_01`–`34_04` ·
`reports/phase4/ROBUSTNESS_REPORT.md`

## Limitation carried forward

n = 4. Nothing in this milestone supports a population-level or epidemiological claim, and
`sample_id` = patient = dataset means biological and technical effects cannot be fully
separated.

## Next milestone

**M35 — freeze.**
