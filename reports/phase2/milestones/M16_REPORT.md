# Milestone M16 — Annotation Refinement and Composition

**Phase 2 · MPNST single-cell integration** · *Generated 2026-09-02*
*Status: **COMPLETE** — proceeded automatically to M17.*

## 1. Scripts

Created: `scripts/R/phase2/refine_and_compose.R`, `scripts/shell/phase2/run_m16_composition.sh`,
`reports/phase2/milestones/M16_REPORT.md`.
Output object: `results/phase2/composition/phase2_harmony_refined.rds`.

## 2. SLURM

| JobID | State | ExitCode | Elapsed | CPUs | Mem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| **19886728** | **COMPLETED** | **0:0** | **00:07:34** | 8 | 128G | **8,890,820 K (8.48 GiB)** — 6.6% |

No failures or retries. R warnings: 0. All 6 round-trip validation checks TRUE.

## 3. Annotation refinement — computational review

Six canonical compartment panels (Malignant/tumour, T/NK, Myeloid, B/Plasma, Endothelial,
Fibroblast/Stromal) were scored per cell with `AddModuleScore` (seed 42) and averaged per
cluster after z-scaling. Each cluster's assigned Level 1 compartment was then compared with
the ranking of its own panel score under four documented rules:

| Rule | Condition | Action |
| --- | --- | --- |
| R1 | Assigned compartment is not one of the scored panels (`Uncertain`, `Other`) | No change — these are deliberate M15 non-assignments |
| R2 | Assigned compartment is the **top-scoring** panel | Confirmed |
| R3 | Assigned ranks 2nd **and** within 0.25 z of the top | Confirmed (adjacent compartments separate weakly) |
| R4 | Otherwise | **Confidence downgraded one step; the LABEL IS RETAINED** |

**Why labels are not overridden.** The M15 labels rest on specific marker genes and cited
literature. A coarse six-panel mean is weaker evidence than that, so it is used to flag
disagreement and lower confidence — never to replace a literature-grounded call.

### Outcome: 18 confirmed · 5 confidence-downgraded · 3 deliberately unscored · **0 labels changed**

| Cluster | Assigned | Top-scoring panel | Rank | Confidence | Interpretation |
| --- | --- | --- | ---: | --- | --- |
| C7 | Malignant/tumour | Fibroblast/Stromal | 6 | Moderate → **Low** | The classical Schwann panel does not support it. Consistent with reported MPNST dedifferentiation (SOX10/MPZ/PMP22 down-regulated), but honestly downgraded. |
| C13 | Malignant/tumour | Fibroblast/Stromal | 6 | Low → Low | Same; already Low. |
| C18 | Malignant/tumour | Fibroblast/Stromal | 6 | Low → Low | Same; already Low. |
| C20 | Cycling tumour-like | Fibroblast/Stromal | 6 | Moderate → **Low** | A dominant G2/M programme masks lineage; the downgrade is appropriate. |
| C21 | Myeloid (pDC) | B/Plasma | 4 | High → **Moderate** | A conservative artefact of panel overlap: pDCs express IGJ and MZB1, which load the B/Plasma panel. The marker-level evidence (LILRA4, CLEC4C, SPIB, GZMB) remains decisive, so the label stands with lowered confidence. |

C12, C15 and C23 fell under R1 (deliberate non-assignments) and were left untouched.

Confidence after refinement: **High 9,828 (49.8%) · Moderate 6,637 (33.7%) · Low 3,251 (16.5%)**
(before: High 10,115 · Moderate 7,650 · Low 1,951).

**Both versions are preserved.** `postint_celltype_level1/2/3_initial` and
`postint_annotation_confidence_initial` hold the M15 values;
`*_refined` hold the M16 values; `postint_annotation_review_rule` and
`postint_annotation_review_outcome` record why. Nothing was silently replaced.

## 4. Composition — descriptive only

> **No inferential condition-level or differential-abundance test was performed.** Cells
> are not independent biological replicates. Every composition table carries this note in a
> `note` column.

`sample_id` is simultaneously the dataset and the patient, so composition by sample,
by patient and by dataset are the same table; all three filenames are written for
convenience and are identical by construction. **No condition table exists** —
`cell_proportions_by_condition_NOT_APPLICABLE.txt` records that no biological-condition
field exists in the data and that fabricating one would be invalid.

**4 of 17 cell types are driven >80% by a single sample/patient**
(`celltype_sample_dependence.tsv`), which is exactly the caution this milestone was asked
to make obvious: those types must not be read as general MPNST biology.

Tables in `results/phase2/composition/`: `annotation_refinement_review.tsv`,
`cell_counts_by_annotation.tsv`, `cell_counts_by_annotation_level1.tsv`,
`cell_proportions_by_sample.tsv`, `cell_proportions_by_patient.tsv`,
`cell_proportions_by_dataset.tsv`, `cell_proportions_level1_by_{sample,patient,dataset}.tsv`,
`celltype_sample_dependence.tsv`, `m16_composition_record.json`, `prov_m16_composition.json`,
`figure_index_m16.tsv`, plus the condition not-applicable note.

## 5. Figures (`results/phase2/figures/M16/`, 10 figures × PDF + PNG)

Final broad-compartment UMAP · final detailed cell-type UMAP · annotation-confidence UMAP ·
broad composition by sample (%) · detailed composition by sample (%) · composition counts by
sample · cell-type abundance overview · **sample/patient contribution per cell type** (with
an 80% dominance line) · patient contribution per compartment · confidence by compartment
after review.

## 6. Next

Proceeded automatically to **M17 — Final Validation, Freeze and Handoff**.

---

# AMENDMENT ADDENDUM — CCC Annotation Audit

*Added 2026-09-03. The M16 refinement above is unchanged. The CCC layer was built in M15A
from this milestone's object (`results/phase2/composition/phase2_harmony_refined.rds`) and is
audited here alongside the detailed annotation.*

## Refined annotation (recap)

18 clusters confirmed, 5 confidence-downgraded, 3 deliberately unscored, **0 labels changed**.
Both versions preserved in `*_initial` and `*_refined`. Confidence after refinement:
High 9,828 (49.8%) · Moderate 6,637 (33.7%) · Low 3,251 (16.5%).

## CCC annotation audit

17 identities. `annotation_ccc` is derived from `postint_celltype_level2_refined`, so the M16
refinement propagates into the CCC layer by construction. The M16 confidence downgrades are
material: **C7 and C20, both downgraded to Low by the panel review, together contribute 38%
of `MPNST-Tumor`** — recorded in `MPNST_TUMOR_COMPOSITION.tsv` so that any communication
result can be re-tested with them excluded.

## Population sizes — `CCC_POPULATION_SIZE_AUDIT.tsv`

Thresholds: ≥100 cells total, present in ≥2 samples with ≥10 cells each.
**14 of 17 identities are READY; all 14 are present in all four samples.**

| Flagged | Cells | Reason |
| --- | ---: | --- |
| `Candidate-Malignant-Unresolved` | 1,231 | Malignant identity unestablished (Low confidence, patient-private evidence) |
| `Uncertain` | 720 | C12 hypoxia-vs-perineurial, C23 low-complexity |
| `Low-quality-excluded` | 438 | Mitochondrial technical artefact (C15) |

Thinnest READY population: `NK`, 106 cells, 11–43 per sample. **Not merged** — the amendment
forbids merging biologically distinct immune types merely because they are small.

## Sample and patient representation

`sample_id` is simultaneously the dataset and the patient, so the by-sample and by-patient
audits are one analysis. Every CCC-ready population appears in all four tumours:
`MPNST-Tumor` 97–2,428 cells per sample (median 448); `Macrophage` 361–1,133; `Fibroblast`
115–2,959. This is what makes sample-aware communication analysis feasible.
Tables: `ccc_counts_by_{sample,patient}.tsv`, `ccc_proportions_by_{sample,patient}.tsv`.

## Condition representation

**None.** No biological-condition field exists in the data.
`ccc_by_condition_NOT_APPLICABLE.txt` records why no condition table or figure was produced,
rather than fabricating one.

## MPNST-Tumor internal composition

SCP-like 1,680 (49.1%, C8+C9, Moderate) · Schwann-lineage tumour-like 1,011 (29.6%, C7, Low)
· NC-like 440 (12.9%, C14, Moderate) · Cycling tumour-like 289 (8.5%, C20, Low). All four
states retained in the final object.

## Potentially problematic populations

1. `Candidate-Malignant-Unresolved` — the largest excluded group; could be promoted on CNV
   evidence, which would raise the tumour fraction.
2. `Fibroblast` — the Mes-NC-like hypothesis means part of this 5,064-cell population might
   belong to the tumour compartment. Not collapsed, because the hypothesis is unresolved.
3. `NK` — thin, noisy in the smaller samples.
4. `CD4-T` / `T-cell-other` — soft boundary from sparse CD4 detection.
5. `Low-quality-excluded` — exclude from all downstream analysis.

## Figures and tables

CCC figure suite: `results/phase2/figures/CCC_annotation/` (14 figures × PDF + PNG).
Audit tables listed above plus `CCC_ANNOTATION_MAPPING.tsv`, `CCC_ANNOTATION_SUMMARY.tsv`,
`CCC_MARKER_SUMMARY.tsv`.

## SLURM accounting

CCC layer: JobID **19893067**, COMPLETED, 00:07:33, MaxRSS 8.44 GiB, 1 warning (CD4 gate
softness), no retries. This milestone's own job (19886728) is unchanged.

## Limitations

Descriptive only — no inferential differential-abundance or condition-level test was
performed anywhere in M16 or M15A. Cells are not independent biological replicates.
