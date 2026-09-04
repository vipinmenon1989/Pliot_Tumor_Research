# Milestone M21 — CCC Concordance

**Phase 3 · MPNST** · *2026-09-03* · Status: **COMPLETE** → proceeded to M22.
Full report: [`reports/phase3/CCC_CONCORDANCE_REPORT.md`](../CCC_CONCORDANCE_REPORT.md)

## Scripts
`scripts/phase3/concordance/build_concordance.R`, `make_ccc_figures.R`

## SLURM
| Job | JobID | State | Elapsed | Note |
| --- | --- | --- | --- | --- |
| attempt 1 | 19896027 | **FAILED** | 00:04:00 | pandas wrote booleans as the strings `"True"`/`"False"`; `sum()` on a character column |
| attempt 2 | (pred.) | **FAILED** | — | `geom_curve` cannot draw autocrine self-loops (identical endpoints) |
| **final** | **19896062** | **COMPLETED** | **00:05:42** | |

Both failures were code defects fixed at source; no resource was increased.

## Result
**502,846 distinct interaction keys · 36,486 supported.**

| Class | n |
| --- | ---: |
| High concordance | 5,848 |
| Moderate concordance | 499 |
| Single-method | 24,654 |
| **Discordant/ambiguous** | **5,485** |

**Resource overlap governs interpretation:** testable by 3 frameworks 1,347; by 2 13,064;
by 1 **488,435**. 97% of keys exist in only one resource, so "single-method" usually reflects
resource non-overlap, not disagreement. Among the 1,347 jointly testable interactions, **763
(57%) were supported by all three**.

Raw scores were never averaged. NicheNet was not counted as a fourth LR framework.
Figures: `results/phase3/figures/M20/M20_01..05`, `M21/M21_06..10`.
