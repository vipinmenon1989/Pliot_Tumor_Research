# Milestone M17 — Phase 2 Validation, Freeze and Handoff

**Phase 2 · MPNST single-cell integration** · *Generated 2026-09-03*
*Status: **COMPLETE — PHASE 2 FROZEN***

> Handoff document: [`reports/phase2/PHASE2_HANDOFF.md`](../PHASE2_HANDOFF.md)
> Manifest: `results/phase2/phase2_manifest.json`

## 1. Scripts

Created: `scripts/R/phase2/validate_and_freeze.R`, `scripts/shell/phase2/run_m17_freeze.sh`,
`reports/phase2/PHASE2_HANDOFF.md`, `reports/phase2/milestones/M17_REPORT.md`.
No new exploratory biology was performed — M17 validates, freezes and documents.

## 2. SLURM

| JobID | State | ExitCode | Elapsed | CPUs | Mem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| **19886743** | **COMPLETED** | **0:0** | **00:07:39** | 8 | 128G | **8,520,016 K (8.13 GiB)** — 6.4% |

No failures. **R warnings: 0.** No check returned FALSE.

## 3. Validation performed

**26 pre-save checks**, all TRUE: cell count 19,716 · no duplicate cell IDs · both assays
present · default assay `SCT` · all four reductions present (`pca`, `umap_preintegration`,
`postint_harmony`, `postint_umap_harmony`) · all four graphs present · all ten resolution
columns present · primary and alternative cluster columns present · all 14 annotation
columns present · Phase 1 metadata preserved · no missing annotations or cluster labels ·
Harmony reduction is 30-dimensional · **no non-finite values in any of the four
embeddings** · initial annotation retained.

**Scientific-scope guards**, all TRUE — the object contains no column matching
`cnv|infercnv|copykat`, `pseudobulk|deseq|edger|condition_de`, or
`pseudotime|velocity|monocle|slingshot`, confirming machine-checkably that no prohibited
Phase 3 analysis leaked into Phase 2.

**Phase 1 immutability**: `processed_mpnst.rds` (mtime 2026-07-06 11:29) and
`results/combined/pre_integration/combined_preintegration.rds` (mtime 2026-07-19 23:25)
untouched; the handoff object's md5 was recomputed and matches `88a442688f912d882f6c6da01820e329`.

**14 post-reload checks**, all TRUE. The final object was written, memory freed, then
**re-read from disk** and re-verified — cells, features, assays, default assay, reductions,
graphs, metadata columns, SCT model count, UMAP embedding digest, annotation digest, all ten
resolution columns, all annotation columns, no duplicate cells, no non-finite values.
`saveRDS()` returning without error was not treated as proof of validity.

## 4. Final object

| Item | Value |
| --- | --- |
| Path | `results/phase2/phase2_final_object.rds` |
| Size | 6,055,929,223 bytes (5.64 GiB) |
| md5 | `63146e84e43d8f036cca6fd6d1ef99b3` |
| sha256 | `51f0833f9b1e76a83c8a94046fd8bffcf647d842c59dae8f7b0705b0e9d9d233` |
| Cells | 19,716 |

Contents are enumerated in `PHASE2_HANDOFF.md` §11.

## 5. Deliverables assembled

- **`results/phase2/phase2_manifest.json`** — Phase 1 input and checksum, M11 object and
  checksum, intermediate objects, final object and checksums, Harmony variable/dimensions/
  parameters (explicit vs default), integration assessment and caveats, all ten clustering
  resolutions and the selected primary/alternative, marker parameters and table paths,
  annotation fields and map, composition scope, figure and table roots, all report paths,
  all script paths, environment and package versions, random seeds, every SLURM JobID with
  requested and actual resources including the four documented failures, git commit/branch/
  status, the full validation results, warnings, 11 known limitations, and the list of
  respected scientific prohibitions.
- **`results/phase2/tables/final/`** — 17 files.
- **`results/phase2/figures/final/`** — 28 curated figure files plus a not-applicable note.
- **`reports/FIGURE_INDEX.tsv`** — 438 rows with `phase`, `milestone` and `slurm_job_id`
  columns; 260 back-filled Phase 1 rows and 178 Phase 2 figure files.

## 6. Phase 2 completion criteria

| Criterion | Status |
| --- | --- |
| Default Harmony result | ✅ M11 |
| Pre/post assessment + quantitative diagnostics | ✅ M12 |
| Harmony neighbour graph, UMAP, resolution sweep | ✅ M13 |
| Primary and alternative resolution | ✅ 1.0 / 0.7 |
| All-cluster and top-marker tables, marker figures | ✅ M14 |
| Broad + detailed annotation, confidence, evidence table | ✅ M15 |
| Conservative MPNST annotation | ✅ graded labels, no CNV used |
| Sample / patient composition | ✅ M16 (condition N/A, documented) |
| Complete figure suite, PDF + key PNG, figure index | ✅ 178 Phase 2 files |
| Scripts, logs, environment, checksums, seeds, git, SLURM accounting | ✅ |
| M12–M17 reports, HARMONY_ASSESSMENT, CLUSTERING_ASSESSMENT, MARKER_REPORT, ANNOTATION_REPORT, PHASE2_HANDOFF, manifest | ✅ |
| PROGRESS.md and CHANGELOG.md current | ✅ |
| No Phase 1 object overwritten | ✅ verified by md5 and mtime |
| No prohibited Phase 3 analysis executed | ✅ machine-checked |

## 7. Unresolved biological questions

1. **Are the fibroblast-programme clusters Mes-NC-like malignant?** Malignant fraction is
   23.6% conservatively, up to ~49% otherwise. Needs CNV inference.
2. **Are C13, C17 and C18 malignant?** Patient-private neural programmes, Low confidence.
3. **What is C12?** Hypoxic malignant or perineurial — both defensible.
4. **Is the residual patient structure biology or batch?** Unanswerable by design.
5. Are the MPNST-G1 (SHH, C14) and G2 (SCP-like, C8/C9) subgroups genuinely present here?

## 8. STOP

**Phase 2 is complete and frozen. Phase 3 was not begun and requires separate
authorization.**

---

# AMENDMENT ADDENDUM — Re-freeze including `annotation_ccc`

*Added 2026-09-03 · SLURM JobID **19893267** (re-freeze) · initial freeze was 19886743*

Phase 2 was re-frozen after the approved CCC-oriented annotation amendment, so the final
object, manifest, table suite and figure suite all include the new layer.

## 1. `annotation_ccc` verification — all TRUE

Eleven CCC-specific pre-save guards were added to `scripts/R/phase2/validate_and_freeze.R`:

| Check | Result |
| --- | --- |
| `ccc_annotation_present` (all four CCC columns) | ✅ |
| `ccc_no_missing_labels` | ✅ |
| `ccc_mpnst_tumor_present` | ✅ |
| `ccc_immune_identities_kept` (≥4 immune identities) | ✅ (9 kept) |
| `ccc_stromal_kept_separate` (Fibroblast and Endothelial both present) | ✅ |
| `ccc_tumor_only_from_malignant_detail` | ✅ |
| `ccc_no_fibro_endo_in_tumor` | ✅ |
| `ccc_no_immune_in_tumor` | ✅ |
| `ccc_no_uncertain_in_tumor` | ✅ |
| `ccc_is_tumor_flag_consistent` | ✅ |
| `detailed_tumor_states_retained` (≥2 detailed states inside MPNST-Tumor) | ✅ (4) |

Plus four post-reload guards: `reload_ccc_present`, `reload_ccc_identical`,
`reload_ccc_mpnst_tumor`, `reload_detailed_and_ccc_coexist` — all TRUE.

**Total: 37 pre-save checks and 18 post-reload checks, none FALSE.** The scientific-scope
guards (no CNV, pseudobulk-DE or trajectory columns) still pass, and Phase 1 immutability was
re-verified by md5 and mtime.

## 2. Final object (re-frozen)

| Item | Value |
| --- | --- |
| Path | `results/phase2/phase2_final_object.rds` |
| Size | 5.64 GiB |
| **md5** | **`153d5f6acc70f9c05aa48cabc4f4ac2d`** |
| **sha256** | **`62e97524836309170c38be9335a10a4c76baa48a246b2963d2d440cebead0cf3`** |
| Cells | 19,716 · Features (default `SCT`) 29,113 · `RNA` 31,764 |
| Assays | `RNA` (Assay5, split layers), `SCT` (SCTAssay, 4 models) — default `SCT` |
| Reductions | `pca`, `umap_preintegration`, `postint_harmony`, `postint_umap_harmony` |
| Graphs | `SCT_nn`, `SCT_snn`, `postint_harmony_nn`, `postint_harmony_snn` |

*(The pre-amendment freeze was md5 `63146e84e43d8f036cca6fd6d1ef99b3`; the object now
additionally carries the four CCC columns.)*

## 3. CCC populations in the final object — 17 identities

`MPNST-Tumor` **3,420 cells (17.35%)**. Immune total **7,448 (37.8%)** across nine
identities: Macrophage 3,065 · Plasma-cell 1,709 · Dendritic 612 · CD4-T 606 · CD8-T 382 ·
Plasmacytoid-DC 287 · T-cell-other 258 · B-cell 233 · Monocyte 190 · NK 106.
Stromal/endothelial **6,459 (32.8%)**: Fibroblast 5,064 · Pericyte-VSMC 435 · Endothelial 960.
Uncertain/excluded **2,389 (12.1%)**: Candidate-Malignant-Unresolved 1,231 · Uncertain 720 ·
Low-quality-excluded 438.

## 4. Annotation metadata in the final object

Detailed (source of truth): `postint_celltype_level1/2/3`, `postint_annotation_confidence`,
`postint_celltype_level1/2/3_initial`, `postint_annotation_confidence_initial`,
`postint_celltype_level1/2/3_refined`, `postint_annotation_confidence_refined`,
`postint_annotation_initial`, `postint_annotation_refined`,
`postint_annotation_review_rule`, `postint_annotation_review_outcome`,
`postint_annotation_source_cluster`.
CCC layer: **`annotation_ccc`**, `annotation_ccc_compartment`, `annotation_ccc_is_tumor`,
`annotation_ccc_ccc_ready`.
Clustering: all ten `postint_harmony_clusters_res_*`, `postint_harmony_primary_cluster`,
`postint_harmony_alternative_cluster`. Plus all 93 original Phase 1 columns verbatim.

## 5. Final deliverables (updated)

- `results/phase2/phase2_manifest.json` — now carries a `ccc_annotation` section with the
  field names, mapping table and config, tumour-collapse criteria, per-population counts and
  proportions, `MPNST-Tumor` cell count and percentage, readiness report path, figures,
  script and record, plus three additional known limitations and the M15A/re-freeze JobIDs.
- `results/phase2/tables/final/` — **27 files**, adding `ccc_annotation_mapping.tsv`,
  `ccc_annotation_summary.tsv`, `ccc_population_size_audit.tsv`,
  `mpnst_tumor_composition.tsv`, `ccc_marker_summary.tsv`,
  `ccc_counts_by_{sample,patient}.tsv`, `ccc_proportions_by_{sample,patient}.tsv`.
- `results/phase2/figures/final/` — **34 files** renumbered to the amendment's requested
  scheme, plus `12_ccc_composition_by_condition_NOT_APPLICABLE.txt`.
- `reports/FIGURE_INDEX.tsv` — **472 rows**, extended with `figure_type`,
  `annotation_field`, `input_object` and `notes` columns; 28 M15A CCC rows and 34 final rows.
- **`reports/phase2/CCC_READINESS.md`** — verdict **`READY WITH CAVEATS`**.
- `PROJECT.md` — new Phase 2 architecture section (P2.1–P2.7).

## 6. Phase 3 readiness

**`READY WITH CAVEATS`.** All 14 CCC-ready populations are present in all four samples, so
sample-aware communication analysis is feasible. Caveats: exclude the three flagged
populations (12.1% of cells); the malignant fraction is 17.35% conservatively and could reach
~49% if the fibroblast clusters prove Mes-NC-like malignant, so **CNV inference should precede
substantive CCC conclusions**; two Low-confidence states contribute 38% of `MPNST-Tumor` and
warrant an exclusion sensitivity check; `NK` is thin; `CD4-T`/`T-cell-other` is a soft
boundary; no condition variable exists. Full detail in `CCC_READINESS.md`.

## 7. SLURM

| Job | JobID | State | ExitCode | Elapsed | CPUs | Mem | MaxRSS |
| --- | --- | --- | --- | --- | ---: | ---: | ---: |
| M15A CCC layer | 19893067 | COMPLETED | 0:0 | 00:07:33 | 8 | 128G | 8.44 GiB |
| M17 re-freeze | 19893267 | COMPLETED | 0:0 | 00:07:59 | 8 | 128G | 8.13 GiB |
| M17 initial freeze | 19886743 | COMPLETED | 0:0 | 00:07:39 | 8 | 128G | 8.13 GiB |

No failures or retries in the amendment. R warnings: 1 in M15A (CD4 gate softness), 0 in the
re-freeze.

## 8. STOP

**Phase 2 is complete and re-frozen with both annotation layers. Phase 3 was not begun and
requires separate authorization.**
