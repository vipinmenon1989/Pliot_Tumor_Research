# Milestone M26 — Robustness and Replication

**Phase 3 · MPNST** · *2026-09-03* · Status: **COMPLETE** → proceeded to M27.
Full report: [`reports/phase3/ROBUSTNESS_REPORT.md`](../ROBUSTNESS_REPORT.md)

## Script
`scripts/phase3/concordance/robustness.R` · SLURM **19896065**, COMPLETED, 00:00:24.

## Leave-one-patient-out
Retention of ≥2-patient support after dropping each tumour: MPNST_1 **71.1%**, MPNST_2 76.3%,
MPNST_3 79.2%, MPNST_4 79.4%. **No single patient dominates.**

## Single-patient dominance
**19,558 of 36,486 (53.6%)** supported interactions rest on one patient (MPNST_1 8,431;
MPNST_2 4,581; MPNST_4 3,394; MPNST_3 3,152). Retained and flagged; excluded from tiers P1–P4
by construction.

## Method disagreement — preserved
**5,485** discordant interactions (testable by ≥2 frameworks, supported by 1), against 5,848
high-concordance — the frameworks disagree about as often as they agree on jointly testable
pairs.

## Rare populations
Monocyte, B-cell, pDC, T-cell-other and NK all have a **median of 1–2 supporting patients**.
None was merged away; all are flagged. NK (106 cells) is least reliable.

## LochNESS robustness
LOSO Spearman 0.31–0.50; representation sensitivity 0.02–0.42 (CD8-T **−0.265**). The score is
not stable at this sample size, reinforcing the M24 negative.

Not a benchmarking study. Negative results and disagreement are kept in full.
