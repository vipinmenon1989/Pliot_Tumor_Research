# M35 — Phase 4 Freeze and Handoff

**Status: COMPLETE** · 2026-09-03 · SLURM **19897640** COMPLETED 00:12:16, ReqMem 96 G,
MaxRSS 13.38 GiB

No new exploratory biology (§56). This milestone validated and froze what M28–M34 produced.

## Final object

```text
path      results/phase4/phase4_final_object.rds
size      5.64 GB
md5       e85ba8486e456917e2483f2773bdbaf3
sha256    a658447744e5ee8621fc5a3127492ceb11bce37bac4d6f3f0f7255bdcb647ae4
cells     19,716     features  31,764 (RNA) / 29,113 (SCT)
assays    RNA, SCT   reductions  pca, umap_preintegration, postint_harmony, postint_umap_harmony
metadata  149 Phase 1-3 columns + 32 Phase 4 columns
```

## Validation — 23 / 23 passed

Preservation asserted **before** saving, column by column rather than by a whole-frame digest
(the lesson from Phase 2 M13, where a blunt digest fired on a change that had not happened):

```text
GUARD PASSED  all 149 pre-existing metadata columns are byte-identical
GUARD PASSED  reductions unchanged
GUARD PASSED  EVERY Phase 2 embedding numerically identical (PCA and Harmony included)
GUARD PASSED  assays and cell order unchanged
GUARD PASSED  annotation_ccc_phase3 is a verbatim copy of the frozen Phase 3 layer
GUARD PASSED  no non-Malignant cell carries the refined MPNST-Tumor label
GUARD PASSED  tumour states exist only on Malignant cells
```

Then reloaded under SLURM and re-checked: cell count, feature count, assays, RNA layers (12),
SCT layers (3), reductions, duplicate IDs (0), duplicate features (0), all Phase 2 metadata
retained (0 missing), all 32 Phase 4 fields present (0 missing), no missing
`malignancy_refined` / `malignancy_confidence` / `annotation_ccc_refined`, valid factor levels,
**0 non-finite values across all four embeddings** (1,656,144 values checked),
`annotation_ccc == annotation_ccc_phase3`, refined `MPNST-Tumor` implies Malignant, Ambiguous
never `MPNST-Tumor`, tumour states only on Malignant cells, and **save/load md5 stable**.

## Deliverables

| Item | Status |
| --- | --- |
| `results/phase4/phase4_manifest.json` | **35 sections**, `figures_missing: []`, `tables_missing: []`, 11 limitations, every JobID including failures |
| `reports/phase4/PHASE4_HANDOFF.md` | 27 sections, readable without conversation history |
| `results/phase4/figures/final/` | **16 of 16** required figures placed + 27 supporting = **86 files**, PDF + PNG |
| `results/phase4/tables/final/` | **47 tables**, all 15 §53-required present |
| `reports/FIGURE_INDEX.tsv` | 578 → **664 rows**, 86 Phase 4 rows registered |
| Milestone reports | M28–M35, eight files |
| Phase 4 reports | 8 files in `reports/phase4/` |

## The CNA heatmap (§20, §52)

`04_scevan_cna_heatmap` was built as a genuine ComplexHeatmap figure rather than substituted:
genome-ordered chromosomes 1–22 with four annotation tracks (frozen Phase 2 annotation, SCEVAN
call, refined call, SCEVAN clone), cells ordered by refined call then clone then Phase 2 label,
seeded 2,000-cell subsample per patient for legibility, and an explicit in-figure banner on
samples that failed the sanity gate. It is the figure that lets a reader see for themselves that
the promoted "fibroblasts" carry the same CNA architecture as the cells Phase 2 already called
malignant.

## SLURM

| JobID | State | Elapsed | ReqMem | MaxRSS |
| --- | --- | --- | --- | --- |
| **19897640** | COMPLETED | 00:12:16 | 96 G | 13.38 GiB |

Three stages, all exit 0: object build + validation, CNA heatmap, assembly. No failures, no
resource increases. Most of the runtime was PNG rasterisation of four large heatmaps without
the optional `magick` backend — an advisory ComplexHeatmap note, not an error.

## Warnings carried into the freeze

Eleven limitations are recorded in the manifest and §20 of the handoff. The three a reader
should see first:

1. **The 32.63% refined malignant fraction is not patient-robust** — dropping MPNST_4 gives
   0.214 against a Phase 2 fraction of 0.215.
2. **No malignant transcriptional state is recurrent** (0 of 8; 97.2% of malignant cells in
   patient-private states), so the tumour-state → TME model could not be built and the §68
   structure was **not** imposed.
3. **Amendment A2 was written after observing that the flat sanity gate excluded the sample
   carrying the strongest evidence.** Both gates are published; the fibroblast conclusion's
   direction survives rejecting A2, its magnitude does not.

## Phase 4 status

**COMPLETE.** Phase 5 not initiated; requires separate authorization.
