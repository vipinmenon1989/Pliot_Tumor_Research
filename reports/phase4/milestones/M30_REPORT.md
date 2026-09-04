# M30 — Malignancy Integration

**Status: COMPLETE** · 2026-09-03 · SLURM **19897148** COMPLETED 00:01:00, MaxRSS 1.38 GiB
(first attempt **19897096** COMPLETED 00:01:06 under the pre-Amendment-A2 gate and was
superseded — see below)

## Input

`results/phase4/malignancy/phase4_cell_metadata.tsv.gz` (frozen Phase 2/3 metadata, M28) ·
the eight M29 SCEVAN classifications · the four per-sample CNA matrices ·
`results/phase4/scevan/by_sample/*/counts_raw.rds` for module scoring.
Decision rules from `reports/phase4/MALIGNANCY_DECISION_RULES.md`, **written before any SCEVAN
result was inspected**, plus Amendments A1 and A2, both dated and labelled post-hoc.

## Output

```text
results/phase4/malignancy/PHASE4_MALIGNANCY_CALLS.tsv        19,716 x 37, one row per cell
results/phase4/malignancy/phase4_malignancy_metadata.rds
results/phase4/malignancy/m30_malignancy_facts.json
results/phase4/tables/SCEVAN_VS_PHASE2_ANNOTATION.tsv        population x patient confusion
results/phase4/tables/AMBIGUOUS_POPULATION_MALIGNANCY_AUDIT.tsv
results/phase4/tables/SCEVAN_SAMPLE_RELIABILITY.tsv          both gates, per sample
results/phase4/tables/SCEVAN_POPULATION_MALIGNANT_FRACTION.tsv
results/phase4/tables/MALIGNANCY_THRESHOLD_SENSITIVITY.tsv
results/phase4/tables/MALIGNANCY_SUMMARY_BY_{CLUSTER,SAMPLE}.tsv
results/phase4/tables/LINEAGE_MARKER_EVIDENCE_BY_POPULATION.tsv
results/phase4/tables/SCEVAN_CELL_CLASSIFICATION{,_RAW}.tsv
results/phase4/tables/SCEVAN_CLONES.tsv · SCEVAN_CNV_SUMMARY.tsv
results/phase4/tables/SCEVAN_ARM_LEVEL_EVENTS.tsv · SCEVAN_RECURRENT_BROAD_EVENTS.tsv
results/phase4/tables/SCEVAN_RUN_ACCOUNTING.tsv · SCEVAN_NATIVE_FIGURE_INDEX.tsv
results/phase4/tables/CNA_PROFILES_BY_ANNOTATION.tsv.gz
reports/phase4/MALIGNANCY_REFINEMENT_REPORT.md
reports/phase4/AMBIGUOUS_POPULATION_MALIGNANCY_REPORT.md
```

## Scripts

`scripts/phase4/scevan/m29_aggregate.R` · `scripts/phase4/malignancy/m30_integrate_malignancy.R`
· `scripts/phase4/malignancy/m30_figures.R` · `scripts/shell/phase4/run_m30.sh`

## Headline

| Stage | Malignant | Fraction |
| --- | ---: | ---: |
| Phase 2 conservative | 3,420 | 17.35% |
| SCEVAN primary | 8,794 / 18,449 assessed | 47.7% of assessed |
| **Phase 4 refined** | **6,434** | **32.63%** |

Also 9,078 Non-malignant, **3,766 Ambiguous (19.10%)**, 438 Excluded. High-confidence
Malignant 3,261.

**Answer to §28**: yes, 17.35% was too conservative — but the correction is two-directional.
Phase 4 adds 4,036 fibroblasts, 836 candidate-malignant, 92 uncertain and 65 pericytes, and
simultaneously withdraws confidence from 2,015 of the 3,420 cells Phase 2 called malignant.
Phase 2 did not so much undercount as look in the wrong place.

## Figures (§27)

All eight required, PDF primary + PNG:
`01_phase2_vs_scevan_malignancy` · `02_scevan_malignancy_umap` · `03_refined_malignancy_umap` ·
`04_malignancy_by_phase2_cluster` · `05_malignancy_by_annotation` · `06_malignancy_by_patient` ·
`07_ambiguous_population_malignancy` · `08_cnv_profiles_by_annotation`.
Plus from M29 aggregation: `29_01_scevan_call_umap_by_patient` ·
`29_02_scevan_clones_umap_by_patient` · `29_03_clone_composition_by_annotation`.

## JobIDs / resources / runtime / MaxRSS

| JobID | State | Elapsed | ReqMem | MaxRSS | Note |
| --- | --- | --- | --- | --- | --- |
| 19897096 | COMPLETED | 00:01:06 | 64 G | 0.57 GiB | superseded by Amendment A2 |
| 19897097 | CANCELLED | 00:00:14 | 96 G | 0.35 GiB | downstream M32-prepare, cancelled with the chain |
| **19897148** | **COMPLETED** | **00:01:00** | 64 G | **1.38 GiB** | reported result |

Also cancelled without running: 19897098–19897103 and 19897104–19897111 (two duplicate
dependency chains submitted by a waiter and a monitor simultaneously; both stopped, no
scientific output from either).

**No failures and no resource increases.** 64 G was heavily over-provisioned for this stage
(1.38 GiB used) because it was sized for the CNA matrices; it was not reduced mid-phase, and
that over-provisioning is recorded rather than hidden.

## Errors

None during execution. Three latent defects were found by review **before** submission and
fixed: `aggregate()` with a formula could not see a variable outside `data=`; an empty marker
panel would have aborted `AddModuleScore`; and a fragile native-pipe lambda in the clone
summary. A fourth (`else` at top level in `m35_assemble.R`) was caught by a parse check.

## Warnings

* **`MPNST-Tumor` retention is the most fragile number in Phase 4.** MPNST_1's pop_frac is
  0.222, just under the a priori 0.25, so a 0.20 threshold would have retained 3,266 instead of
  1,405 cells. Disclosed in `MALIGNANCY_THRESHOLD_SENSITIVITY.tsv`, not discovered later.
* The fibroblast finding is by contrast **completely threshold-independent** (4,036 at every
  grid point from 0.15 to 0.40).
* `Fibroblast` falls to 1 cell in MPNST_4 and `Pericyte-VSMC` to 4, so those become NOT
  EVALUABLE there; `Plasmacytoid-DC` and `Uncertain` vanish from some patients. M32 must
  therefore distinguish "lost" from "no longer testable".

## Scientific decisions

1. **Amendment A1** (post-hoc, labelled): a per-sample immune sanity gate; a run calling >25%
   of canonical immune cells malignant cannot drive malignant **promotions**. Asymmetric on
   purpose — an unreliable CNV run may still support a `Non-malignant` call, because those
   always require concordant lineage evidence too.
2. **Amendment A2** (post-hoc, labelled, and prompted by the gate excluding the strongest
   result): plasma cells excluded from the sanity denominator on mechanistic grounds
   (immunoglobulin loci IGH 14q32 / IGK 2p11 / IGL 22q11 read as segmental gain — a documented
   artefact class for this method family); a **breadth** criterion added so a global inversion
   cannot dilute past the gate; and rule **R7p** so 1,125 plasma cells are called
   `Non-malignant` with the artefact named rather than dropped into `Ambiguous`. Both gates are
   reported for all four samples with a `flat_gate_would_have_failed` column.
3. Marker evidence can only **reduce** confidence, never raise it — letting the markers Phase 2
   already used raise confidence in a Phase 4 call would reintroduce circularity.
4. `Ambiguous` treated as a terminal outcome. 3,766 cells stay unresolved and are excluded from
   the refined CCC analysis as NOT EVALUABLE.
5. Threshold sensitivity computed and published rather than the single a priori value alone.

## Next milestone

**M32 — refined CCC input preparation and targeted sensitivity analysis** (chain 19897175 →
19897181), reusing the Phase 3 LIANA / CellChat / CellPhoneDB / concordance scripts
**completely unmodified**.
