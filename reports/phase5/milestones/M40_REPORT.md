# M40 — Phase 5 Program Robustness

**Date** 2026-09-03 · **Status** COMPLETE · **SLURM** 19899355 (LOO + high-confidence cNMF),
19899358 (technical-gene-free cNMF), 19899359 (robustness analysis)

**Scientific question.** Which of the eight programs survive perturbation of the rank, the patient
balance, the gene set, the cell-cycle genes, the confidence filter, and the removal of the patient
that dominates them?

## Six perturbations, plus two added because the results demanded them

| test | design |
| --- | --- |
| program-rank sensitivity | best-match cosine of each K = 8 program against the K = 7 and K = 9 solutions |
| patient-balanced | every patient downsampled to 651 malignant cells (dominant fraction 0.573 → 0.301) |
| gene-selection | 10 independent random 80% gene subsets, re-projected, Spearman vs the primary usage |
| cell-cycle | the 95 canonical Tirosh genes removed and the grid re-factorized (`nocc`) |
| malignancy confidence | recurrence recomputed on High-confidence cells only, plus a `highconf` cNMF run |
| leave-one-patient-out | four cNMF runs, each withholding one patient entirely |
| **technical covariates** *(added)* | within-patient Spearman of each program's usage against nCount_RNA, nFeature_RNA, percent.mt |
| **technical-gene-free universe** *(added)* | ribosomal, pseudogene/lncRNA and canonical myeloid genes removed; the full K grid re-run so this universe selects its own K |

The last two were added after M38 showed that **three of the eight programs are
technical-dominated** and that the only cross-patient program is one of them. They are declared
sensitivity analyses of the primary result, not a replacement for it, and both are reported in full
whichever way they come out.

The LOO and high-confidence runs factorize **only the selected K with 50 replicates** rather than
the whole grid at 100 — they exist to ask whether the primary programs are *recoverable*, not to
select a rank a second time. That is an evidence-based reduction in cost, not a shortcut past a
needed analysis.

## Declared pass thresholds

```text
patient-balanced best-match cosine        >= 0.60
gene-subset minimum Spearman rho          >= 0.80
leave-one-out worst-case cosine           >= 0.50
cosine when the DOMINANT patient is out   >= 0.60      <- the decisive one
recurrence class unchanged on High-only   required
max |within-patient rho| vs a technical covariate  < 0.50
```

A program failing any one of these is reported as **not robust** and kept in the table, never
dropped from it.

**Outputs.** `PHASE5_PROGRAM_ROBUSTNESS.tsv` · `PROGRAM_NOTECH_RECURRENCE.tsv` ·
`M40_LEAVE_ONE_PATIENT_OUT_RECOVERY.tsv` · `M40_PROGRAM_RANK_SENSITIVITY.tsv` ·
`M40_GENE_SUBSET_USAGE_CORRELATION.tsv` · `M40_PROGRAM_TECHNICAL_COVARIATES.tsv` ·
`reports/phase5/ROBUSTNESS_REPORT.md`

**Resources.** 19899355 COMPLETED 00:08:15, 8 CPUs, 48 G, MaxRSS 1.83 GiB (5 cNMF runs).
19899358 and 19899359 recorded in the handoff SLURM table.
