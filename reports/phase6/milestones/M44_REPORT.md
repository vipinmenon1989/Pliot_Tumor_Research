# M44 — Within-Clone Transcriptional-Program Diversity

**Date** 2026-09-03 · **Status** COMPLETE · **SLURM** 19899748

**Wording rule, enforced throughout.** This measures **diversity**, not switching. Cross-sectional
scRNA-seq cannot demonstrate a transition, so no transition rate, direction or trajectory is claimed
anywhere in this milestone.

**Design.** Diversity is measured on **both** the continuous program scores and a cautiously defined
hard assignment whose threshold is printed: a cell keeps its top program only when it leads the
runner-up by ≥ 0.10 relative usage, otherwise it is `Mixed`. Several metrics are used because none
is definitive: Shannon entropy, effective number of programs (exp of entropy), dominant-program
share, program dispersion around the clone centroid, and the fraction of Mixed cells.

**Result — diversity is real but strongly patient-specific.**

| patient | clones | effective programs per clone (hard) | dominant-program share | Mixed cells |
| --- | ---: | --- | --- | --- |
| MPNST_1 | 7 | **1.69 – 2.83** | 0.43 – 0.78 | 3 – 25% |
| MPNST_2 | 4 | **1.00** in all four | **1.000** in all four | 0% |
| MPNST_4 | 7 | 1.05 – 1.23 | 0.95 – 0.99 | 0 – 1% |

**Between-clone divergence is negligible next to within-clone spread.** Mean Jensen–Shannon
divergence between a patient's clone profiles against mean within-clone dispersion:

| patient | between-clone JSD | within-clone dispersion | ratio |
| --- | ---: | ---: | ---: |
| MPNST_1 | 0.0130 | 0.384 | **0.034** |
| MPNST_2 | 0.00018 | 0.129 | **0.0014** |
| MPNST_4 | 0.00085 | 0.123 | **0.0069** |

A patient's CNA-defined clones are **transcriptionally near-interchangeable**.

**Technical confounders, assessed before any plasticity wording.** One correlation reaches the
disqualifying |ρ| ≥ 0.7 bar and is reported rather than buried:

| covariate | continuous effective n | hard effective n | program dispersion |
| --- | ---: | ---: | ---: |
| fraction High-confidence | 0.62 | 0.56 | **0.72** |
| median cell-cycle program | 0.51 | 0.39 | 0.26 |
| median nFeature_RNA | 0.48 | 0.43 | 0.53 |
| median nCount_RNA | 0.18 | 0.07 | 0.23 |
| clone size | 0.03 | 0.17 | −0.12 |

**Program dispersion tracks the fraction of High-confidence malignant cells at ρ = 0.72.** That is
the correlation a reader should weigh before accepting any plasticity reading of the dispersion
metric. The two *effective-number* metrics stay below the bar (0.56–0.62), and clone size and
sequencing depth are not drivers.

**Clone-size sensitivity.** Median effective programs per clone is 1.22 / 1.22 / 1.22 at minimum
clone sizes of 20 / 50 / 100 cells — the estimate does not depend on the threshold.

**Verdict.** 4 of 18 evaluable clones span ≥ 2 programs, all four in MPNST_1. The permitted
statement is that **those clones contain cells spanning multiple malignant transcriptional
programs**; reading that as plasticity is possible only with the ρ = 0.72 confound in view.

**Outputs.** `CLONE_PROGRAM_DIVERSITY.tsv` · `PLASTICITY_METRICS.tsv` ·
`M44_BETWEEN_VS_WITHIN_CLONE.tsv` · `M44_DIVERSITY_CONFOUNDERS.tsv` ·
`M44_CLONE_SIZE_SENSITIVITY.tsv`

**Next.** M45 — regulatory architecture.
