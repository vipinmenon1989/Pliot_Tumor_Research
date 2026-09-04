# Robustness and Replication Report

**Phase 3 · Milestone M26 · MPNST** · *Generated 2026-09-03 · SLURM JobID 19896065*

Not a benchmarking study. The question is whether Phase 3 conclusions survive dropping any
single patient, and where the frameworks genuinely disagree.

---

## 1. Leave-one-patient-out

Of the 16,928 interactions supported in ≥2 patients, how many still reach ≥2-patient support
after each patient is dropped:

| Patient dropped | Retained | % | Lost |
| --- | ---: | ---: | ---: |
| MPNST_1 | 12,036 / 16,928 | **71.1%** | 4,892 |
| MPNST_2 | 12,914 / 16,928 | 76.3% | 4,014 |
| MPNST_3 | 13,400 / 16,928 | 79.2% | 3,528 |
| MPNST_4 | 13,438 / 16,928 | 79.4% | 3,490 |

**No single patient dominates.** Retention is 71–79% whichever tumour is removed, and MPNST_1
— the largest sample at 5,905 CCC cells — has the greatest influence, as expected from cell
numbers rather than from any anomaly.

## 2. Single-patient dominance

**19,558 of 36,486 supported interactions (53.6%) rest on a single patient.**

| Sole supporting patient | Interactions |
| --- | ---: |
| MPNST_1 | 8,431 |
| MPNST_2 | 4,581 |
| MPNST_4 | 3,394 |
| MPNST_3 | 3,152 |

These are **retained and flagged**, never promoted. The prioritisation tiers exclude them from
P1–P4 by construction, since every tier above P5 requires ≥2-patient support. The tumour-centric
subset is listed in `CCC_TUMOUR_SINGLE_SAMPLE_INTERACTIONS.tsv`.

That more than half of all supported interactions are patient-private is the single most
important robustness statistic in Phase 3, and it is why recurrence is a prioritisation
criterion rather than an afterthought.

## 3. Method disagreement — preserved

**5,485 interactions are testable by ≥2 frameworks but supported by only one.** These are
classified `Discordant/ambiguous` and kept in full
(`CCC_METHOD_DISAGREEMENT.tsv`). Disagreement is a result, not noise to be averaged away.

Set against 5,848 `High concordance` interactions, the frameworks disagree about as often as
they agree on jointly testable pairs — a useful calibration on how much weight any single-tool
CCC result deserves.

## 4. Rare-population sensitivity

| Population | Interactions | ≥2 frameworks | ≥2 patients | Median patients |
| --- | ---: | ---: | ---: | ---: |
| Monocyte (190 cells) | 4,366 | 756 | 2,271 | **2** |
| B-cell (233) | 2,794 | 447 | 708 | 1 |
| Plasmacytoid-DC (287) | 2,753 | 507 | 999 | 1 |
| T-cell-other (258) | 2,362 | 373 | 1,172 | 1 |
| NK (106) | 2,436 | 454 | 1,151 | 1 |

Every rare population has a **median of 1–2 supporting patients**, versus the cohort-wide
picture where 46% reach ≥2. Interactions involving these populations are therefore
systematically less reproducible, and NK — the smallest CCC-ready population — is the least
reliable. None was merged away; they are flagged, per the Phase 2 rule against merging
biologically distinct immune types for size alone.

## 5. LochNESS robustness

Leave-one-reference-sample-out mean pairwise Spearman: Macrophage 0.497, Endothelial 0.448,
Fibroblast 0.347, CD8-T 0.310. Representation sensitivity (Harmony vs within-lineage PCA):
0.419, 0.219, 0.023, −0.265. **The LochNESS score is not stable** at this sample size, which
reinforces the negative result in `LOCHNESS_MPNST_REPORT.md` rather than contradicting it.

## 6. What survives everything

The interactions that survive method concordance, 4/4-patient recurrence, expression support
**and** receiver-response evidence are the 325 four-stream interactions, of which the
tumour-centric core is:

- **VEGFA → KDR / FLT1 / NRP1** (Endothelial)
- **HLA-E → KLRC1** (NK), **HLA-F → LILRB1/LILRB2** (Monocyte)
- **HLA-A / HLA-E → CD8A / CD8B** (CD8-T)
- **FGF2 → FGFR1** (Fibroblast)
- **COL1A2 → ITGA1/2/3/9_ITGB1** (Endothelial)

Notably **APP → CD74 is not in this set** — it has 3-framework, 4/4-patient and expression
support but lacks NicheNet receiver-response evidence (see `RECEIVER_RESPONSE_REPORT.md` §3).

## 7. Files

`results/phase3/tables/CCC_LEAVE_ONE_SAMPLE_OUT.tsv` ·
`CCC_SINGLE_SAMPLE_DOMINANCE.tsv` · `CCC_TUMOUR_SINGLE_SAMPLE_INTERACTIONS.tsv` ·
`CCC_METHOD_DISAGREEMENT.tsv` · `CCC_METHOD_DISAGREEMENT_SUMMARY.tsv` ·
`CCC_RARE_POPULATION_SENSITIVITY.tsv` · figures `results/phase3/figures/M26/`.
