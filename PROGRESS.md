# MPNST Project Progress

**PHASE 1 = FROZEN**
**PHASE 2 = COMPLETE** — M10–M17 complete, plus milestone M15A (CCC-oriented annotation amendment). Re-frozen 2026-09-03. Phase 3 not begun and requires separate authorization.

---

# Phase 2 — Harmony Integration, Clustering, Markers and Annotation
---

## AMENDMENT — CCC-Oriented Annotation Layer (authorized 2026-09-03)

**Date/time authorized**: 2026-09-03
**Milestone status at authorization**: M10–M17 all **complete**; Phase 2 was frozen at
final object `results/phase2/phase2_final_object.rds` (md5 `63146e84e43d8f036cca6fd6d1ef99b3`).
**Amendment status**: **COMPLETE** (M15A JobID 19893067; Phase 2 re-frozen by M17 JobID 19893267)

### Rationale

The detailed M15/M16 annotation resolves the malignant compartment into five distinct
states across eight clusters. That resolution is scientifically correct, but for downstream
cell–cell communication analysis it fragments the tumour compartment into populations too
small and too patient-specific to serve as a communication source or receiver. A simplified
layer is therefore needed in which the evidence-supported malignant states are collapsed
into a single `MPNST-Tumor` population, while immune, stromal and endothelial populations
keep their biologically meaningful identities so that tumour↔immune, tumour↔stromal and
tumour↔endothelial axes can each be evaluated.

The detailed annotation **remains the source of truth and is not replaced or destroyed**.

### Work incorporated

Because M15–M17 were already complete, the amendment is applied as a new milestone
**M15A** consuming the validated M16 object, followed by a **re-run of M17** so that the
final object, manifest and handoff include the new layer. Per §36 of the amendment,
**Harmony (M11), the integration assessment (M12), clustering (M13) and marker discovery
(M14) are NOT rerun** — their outputs are valid and are reused.

### Planned metadata fields

`annotation_ccc` · `annotation_ccc_compartment` · `annotation_ccc_is_tumor` ·
`annotation_ccc_ccc_ready`. No existing annotation column is overwritten.

### Collapse rule (evidence-based, not "non-immune = tumour")

Applied to the detailed Level 2 label: `MPNST-like malignant (SCP-like)`,
`MPNST-like malignant (NC-like)`, `Schwann-lineage tumour-like` and `Cycling tumour-like`
collapse into `MPNST-Tumor`. `Candidate malignant` does **not** — it stays separate as
`Candidate-Malignant-Unresolved`, because the label is provisional, its confidence is Low
and its only evidence is patient-private neural expression. Uncertain and low-quality
clusters are never absorbed. Full rationale per cluster, including every exclusion, in
`results/phase2/annotation/CCC_ANNOTATION_MAPPING.tsv`.

### Expected outputs

- **Tables**: `CCC_ANNOTATION_MAPPING.tsv`, `CCC_ANNOTATION_SUMMARY.tsv`,
  `CCC_POPULATION_SIZE_AUDIT.tsv`, `MPNST_TUMOR_COMPOSITION.tsv`, `CCC_MARKER_SUMMARY.tsv`,
  `ccc_counts_by_{sample,patient}.tsv`, `ccc_proportions_by_{sample,patient}.tsv`
  (no condition table — no such field exists), updated `ANNOTATION_EVIDENCE.tsv`
- **Figures**: 15-figure CCC suite in `results/phase2/figures/CCC_annotation/` (PDF + PNG)
- **Scripts**: `scripts/R/phase2/build_ccc_annotation.R`,
  `scripts/shell/phase2/run_m15a_ccc_annotation.sh`, and an amended
  `scripts/R/phase2/validate_and_freeze.R`
- **Reports**: `reports/phase2/CCC_READINESS.md`; updated M15/M16/M17 reports, handoff,
  manifest, `PROJECT.md`, `reports/FIGURE_INDEX.tsv`
- **Object**: `results/phase2/annotation/phase2_harmony_ccc.rds`, then a re-frozen
  `results/phase2/phase2_final_object.rds`

### Expected downstream use

Phase 3 sample-aware cell–cell communication (LIANA consensus / NicheNet receiver-response),
computed per sample and only then compared across samples. **No CCC method is executed in
Phase 2.**

---


**Current Phase**: Phase 2
**Completed milestone**: **M17 — Phase 2 Validation, Freeze and Handoff** (2026-09-03; re-frozen after the M15A amendment)
**Next milestone**: none — Phase 2 is frozen. Phase 3 requires separate authorization.
**Final object**: `results/phase2/phase2_final_object.rds` (md5 `153d5f6acc70f9c05aa48cabc4f4ac2d`) — contains **both** the detailed annotation and `annotation_ccc`
**Handoff**: [reports/phase2/PHASE2_HANDOFF.md](reports/phase2/PHASE2_HANDOFF.md) · **CCC readiness**: [reports/phase2/CCC_READINESS.md](reports/phase2/CCC_READINESS.md) · **Manifest**: `results/phase2/phase2_manifest.json`
**Harmony assessment status**: **`ACCEPT DEFAULT HARMONY WITH CAVEATS`** — see
[reports/phase2/HARMONY_ASSESSMENT.md](reports/phase2/HARMONY_ASSESSMENT.md)
**STOP status**: **PHASE 2 FROZEN.** No Phase 3 analysis was executed.

| Milestone | Description | Status | Completion Date |
| --- | --- | --- | --- |
| **M10** | Phase 2 Reconstruction and Integration Planning | **Completed** | 2026-09-02 |
| **M11** | Default Harmony Integration | **Completed** | 2026-09-02 |
| **M12** | Pre/Post Harmony Evaluation | **Completed** | 2026-09-02 |
| **M13** | Post-Harmony Neighbours, UMAP and Clustering Sweep | **Completed** | 2026-09-02 |
| **M14** | Marker Discovery | **Completed** | 2026-09-02 |
| **M15** | Literature-Grounded Cell-Type Annotation | **Completed** | 2026-09-02 |
| **M16** | Annotation Refinement and Composition | **Completed** | 2026-09-02 |
| **M15A** | CCC-Oriented Annotation Layer (amendment) | **Completed** | 2026-09-03 |
| **M17** | Phase 2 Validation, Freeze and Handoff | **Completed** | 2026-09-03 |

## M15A — CCC-Oriented Annotation Layer (approved amendment)

- **Status**: Completed (2026-09-03). Reports: [M15_REPORT.md addendum](reports/phase2/milestones/M15_REPORT.md), [M16_REPORT.md addendum](reports/phase2/milestones/M16_REPORT.md), [CCC_READINESS.md](reports/phase2/CCC_READINESS.md)
- **Input**: `results/phase2/composition/phase2_harmony_refined.rds` (M16) · **Output**: `results/phase2/annotation/phase2_harmony_ccc.rds`
- **Harmony, M12, M13 and M14 were NOT rerun** — their outputs are valid and were reused.

### New metadata

`annotation_ccc` · `annotation_ccc_compartment` · `annotation_ccc_is_tumor` · `annotation_ccc_ccc_ready`. **No existing annotation column was overwritten** — every detailed column was verified byte-identical.

### `annotation_ccc` — 17 identities

| Identity | Cells | % | Compartment | CCC-ready |
| --- | ---: | ---: | --- | :---: |
| Fibroblast | 5,064 | 25.68% | Fibroblast-Stromal | ✅ |
| **MPNST-Tumor** | **3,420** | **17.35%** | MPNST-Tumor | ✅ |
| Macrophage | 3,065 | 15.55% | Immune | ✅ |
| Plasma-cell | 1,709 | 8.67% | Immune | ✅ |
| Candidate-Malignant-Unresolved | 1,231 | 6.24% | Uncertain | ❌ |
| Endothelial | 960 | 4.87% | Endothelial | ✅ |
| Uncertain | 720 | 3.65% | Uncertain | ❌ |
| Dendritic | 612 | 3.10% | Immune | ✅ |
| CD4-T | 606 | 3.07% | Immune | ✅ |
| Low-quality-excluded | 438 | 2.22% | Other | ❌ |
| Pericyte-VSMC | 435 | 2.21% | Fibroblast-Stromal | ✅ |
| CD8-T | 382 | 1.94% | Immune | ✅ |
| Plasmacytoid-DC | 287 | 1.46% | Immune | ✅ |
| T-cell-other | 258 | 1.31% | Immune | ✅ |
| B-cell | 233 | 1.18% | Immune | ✅ |
| Monocyte | 190 | 0.96% | Immune | ✅ |
| NK | 106 | 0.54% | Immune | ✅ |

### Scientific decisions

- **Collapse keyed on the detailed Level 2 label**, applied mechanically from `config/phase2/ccc_annotation_map.tsv`. Collapsed into `MPNST-Tumor`: C8, C9 (SCP-like, L1CAM/MPZ), C14 (NC-like, SHH), C7 (Schwann-lineage tumour-like), C20 (cycling tumour-like).
- **`Candidate malignant` (C13, C17, C18; 1,231 cells) deliberately NOT collapsed** — provisional label, Low confidence, patient-private evidence only. Retained as `Candidate-Malignant-Unresolved` so it is visible, excludable and promotable on future CNV evidence.
- **The prohibited inference `if (!immune) -> MPNST-Tumor` was never used.** M17 machine-verifies that no fibroblast, endothelial, immune, uncertain or low-quality cell entered `MPNST-Tumor`.
- **Immune compartment not collapsed**: nine identities retained. The single T/NK cluster C4 was sub-resolved into CD8-T (382), CD4-T (606), NK (106) and T-cell-other (258) by canonical marker gating restricted to that cluster; no unsupported subtype was invented.
- **Tumour heterogeneity preserved**: all four detailed malignant states remain in the final object, so tumour-state-specific CCC is possible later without re-annotating.
- **No marker rediscovery** — `CCC_MARKER_SUMMARY.tsv` reuses the M14 tables; figures use mean SCT expression, never Harmony coordinates.

### Audits

**14 of 17 identities are CCC-READY, and all 14 are present in all four samples** (`MPNST-Tumor` 97–2,428 cells per sample, median 448). Three flagged EXCLUDE (12.1% of cells). Thinnest ready population `NK` (106 cells) flagged, **not merged**. `MPNST-Tumor` internal composition: SCP-like 49.1%, Schwann-lineage tumour-like 29.6%, NC-like 12.9%, cycling 8.5%.

### Files

Tables: `CCC_ANNOTATION_MAPPING.tsv`, `CCC_ANNOTATION_SUMMARY.tsv`, `CCC_POPULATION_SIZE_AUDIT.tsv`, `MPNST_TUMOR_COMPOSITION.tsv`, `CCC_MARKER_SUMMARY.tsv`, `ccc_counts_by_{sample,patient}.tsv`, `ccc_proportions_by_{sample,patient}.tsv`, augmented `ANNOTATION_EVIDENCE.tsv`. Figures: 14 × PDF + PNG in `results/phase2/figures/CCC_annotation/` (figure 13 by condition N/A, documented). Scripts: `scripts/R/phase2/build_ccc_annotation.R`, `scripts/shell/phase2/run_m15a_ccc_annotation.sh`, `config/phase2/ccc_annotation_map.tsv`.

### SLURM

| JobID | State | ExitCode | Elapsed | CPUs | Mem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| **19893067** | **COMPLETED** | **0:0** | **00:07:33** | 8 | 128G | **8.44 GiB** (6.6%) |

No failures or retries. 1 R warning (CD4 transcript detection is sparse, so the CD4-T/T-cell-other boundary is soft). All 16 round-trip checks TRUE.

### Warnings and limitations

Malignant fraction is 17.35% conservatively and could reach ~49% if the fibroblast clusters prove Mes-NC-like malignant — CNV inference required. Two Low-confidence states (C7, C20) contribute 38% of `MPNST-Tumor`. `NK` is thin. `CD4-T`/`T-cell-other` boundary is soft.

---

## M17 — Phase 2 Validation, Freeze and Handoff

- **Status**: Completed (2026-09-03). Reports: [M17_REPORT.md](reports/phase2/milestones/M17_REPORT.md), [PHASE2_HANDOFF.md](reports/phase2/PHASE2_HANDOFF.md)
- **Final object**: `results/phase2/phase2_final_object.rds` — 6,055,929,223 bytes, md5 `63146e84e43d8f036cca6fd6d1ef99b3`, sha256 `51f0833f9b1e76a83c8a94046fd8bffcf647d842c59dae8f7b0705b0e9d9d233`, 19,716 cells
- **Manifest**: `results/phase2/phase2_manifest.json`

### Validation

**26 pre-save checks** and **14 post-reload checks**, all TRUE. The object was written, memory freed, then **re-read from disk** and re-verified — `saveRDS()` returning without error was not treated as proof of validity. No non-finite values in any of the four embeddings.

**Scientific-scope guards passed**: the object contains no column matching `cnv|infercnv|copykat`, `pseudobulk|deseq|edger|condition_de` or `pseudotime|velocity|monocle|slingshot`, confirming machine-checkably that no prohibited Phase 3 analysis leaked into Phase 2.

**Phase 1 immutability verified**: `processed_mpnst.rds` (mtime 2026-07-06 11:29) and the Phase 1 handoff object (mtime 2026-07-19 23:25, md5 re-verified as `88a442688f912d882f6c6da01820e329`) untouched.

### Deliverables

`results/phase2/phase2_manifest.json` (Phase 1 input, M11 object, final object, all checksums, Harmony parameters explicit vs default, integration assessment and caveats, clustering resolutions and selection, marker parameters, annotation fields, composition scope, figures, tables, reports, scripts, environment, seeds, every SLURM JobID with requested and actual resources including four documented failures, git state, full validation results, 11 known limitations, respected prohibitions) · `results/phase2/tables/final/` (17 files) · `results/phase2/figures/final/` (28 files + a not-applicable note) · `reports/FIGURE_INDEX.tsv` (438 rows).

### SLURM

| JobID | State | ExitCode | Elapsed | CPUs | Mem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| **19886743** | **COMPLETED** | **0:0** | **00:07:39** | 8 | 128G | **8.13 GiB** (6.4%) |

R warnings: 0. No check returned FALSE.

### Unresolved biological questions

1. Are the fibroblast-programme clusters Mes-NC-like malignant? (23.6% vs ~49% malignant fraction — needs CNV inference.)
2. Are C13/C17/C18 malignant? Patient-private neural programmes, Low confidence.
3. What is C12 — hypoxic malignant or perineurial? Both defensible.
4. Is the residual patient structure biology or batch? Unanswerable by design.
5. Are the MPNST-G1 (SHH, C14) and G2 (SCP-like, C8/C9) subgroups genuinely present?

### STOP

**PHASE 2 COMPLETE AND FROZEN.** Phase 3 was not begun and requires separate authorization.

---

## M16 — Annotation Refinement and Composition

- **Status**: Completed (2026-09-02). Report: [M16_REPORT.md](reports/phase2/milestones/M16_REPORT.md)
- **Object**: `results/phase2/composition/phase2_harmony_refined.rds`

### Refinement — 18 confirmed, 5 confidence-downgraded, 3 deliberately unscored, **0 labels changed**

Six canonical compartment panels were scored per cell and averaged per cluster; each cluster's assigned Level 1 was compared with its own panel ranking under four documented rules (R1 unscored/deliberate; R2 assigned = top-scoring → confirmed; R3 assigned ranks 2nd within 0.25 z → confirmed; R4 otherwise → **confidence downgraded one step, label retained**).

**Labels are never overridden by the panel score**: M15 labels rest on specific marker genes and cited literature, which is stronger evidence than a coarse six-panel mean.

| Cluster | Assigned | Top panel | Rank | Confidence | Note |
| --- | --- | --- | ---: | --- | --- |
| C7 | Malignant/tumour | Fibroblast/Stromal | 6 | Moderate → **Low** | Classical Schwann panel unsupportive; consistent with MPNST dedifferentiation but honestly downgraded |
| C13, C18 | Malignant/tumour | Fibroblast/Stromal | 6 | Low → Low | Same |
| C20 | Cycling tumour-like | Fibroblast/Stromal | 6 | Moderate → **Low** | G2/M programme masks lineage |
| C21 | Myeloid (pDC) | B/Plasma | 4 | High → **Moderate** | Conservative artefact of panel overlap (pDCs express IGJ/MZB1); marker evidence LILRA4/CLEC4C/SPIB/GZMB remains decisive |

Confidence after refinement: **High 9,828 (49.8%) · Moderate 6,637 (33.7%) · Low 3,251 (16.5%)**. Both versions preserved in `*_initial` and `*_refined` columns, with `postint_annotation_review_rule`/`_outcome` recording why.

### Composition — descriptive only

**No inferential condition-level or differential-abundance test was performed**; every table carries that note in a `note` column. `sample_id` is simultaneously dataset and patient, so those three groupings are one table (all filenames written for convenience). No condition table exists — `cell_proportions_by_condition_NOT_APPLICABLE.txt` records why.

**4 of 17 cell types are driven >80% by a single sample/patient** (`celltype_sample_dependence.tsv`) and must not be read as general MPNST biology.

### SLURM

| JobID | State | ExitCode | Elapsed | CPUs | Mem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| **19886728** | **COMPLETED** | **0:0** | **00:07:34** | 8 | 128G | **8.48 GiB** (6.6%) |

No failures. R warnings: 0. All 6 round-trip checks TRUE.

### Figures

`results/phase2/figures/M16/` — 10 figures × PDF + PNG, including the sample/patient contribution per cell type with an 80% dominance line.

---

## M15 — Literature-Grounded Cell-Type Annotation

- **Status**: Completed (2026-09-02). Reports: [M15_REPORT.md](reports/phase2/milestones/M15_REPORT.md), [ANNOTATION_REPORT.md](reports/phase2/ANNOTATION_REPORT.md)
- **Evidence table**: `results/phase2/annotation/ANNOTATION_EVIDENCE.tsv` · **Annotation map**: `config/phase2/annotation_map_M15.tsv`
- **Object**: `results/phase2/annotation/phase2_harmony_annotated.rds`

### Result

| Level 1 compartment | Cells | % |
| --- | ---: | ---: |
| Fibroblast/Stromal | 5,499 | 27.9% |
| **Malignant / tumour** | **4,651** | **23.6%** |
| Myeloid | 4,154 | 21.1% |
| B/Plasma | 1,942 | 9.8% |
| T/NK | 1,352 | 6.9% |
| Endothelial | 960 | 4.9% |
| Uncertain | 720 | 3.7% |
| Other (technical) | 438 | 2.2% |

17 Level 2 cell types. Confidence: High 10,115 (51.3%), Moderate 7,650 (38.8%), Low 1,951 (9.9%).

### Scientific decisions

- **Annotation is data, not code**: every label with its markers, conflicting evidence, confidence and citation lives in `config/phase2/annotation_map_M15.tsv`; the script applies the map and infers nothing.
- **No CNV inference** (outside Phase 2 scope). Malignant calls integrate Schwann/neural-crest markers, MPNST literature, absence of immune/stromal identity and provenance — never provenance alone, never "non-immune ⇒ tumour".
- **Graded malignant labels**: MPNST-like malignant SCP-like (1,680; C9 retains **MPZ**, C8 expresses **L1CAM**), MPNST-like malignant NC-like (440; C14 expresses **SHH**), Schwann-lineage tumour-like (1,011), Candidate malignant (1,231, Low confidence), Cycling tumour-like (289).
- **Largest open question stated explicitly**: the four fibroblast clusters (5,064 cells) may contain MPNST Mes-NC-like malignant cells. **Malignant fraction is 23.6% conservatively, up to ~49% if those clusters are malignant.** Resolving this needs CNV inference — the top Phase 3 priority.
- **Deliberate non-assignments retained**: C12 (hypoxia vs perineurial GLUT1/ITGB4 — both defensible) and C23 (ribosomal pseudogenes) left Uncertain; C15 labelled technical (mitochondrial-high, matching the Phase 1 MPNST_4 flag).
- **Marker evidence overrode a prior**: the coarse M12 programme score called C21 B/plasma; LILRA4/CLEC4C/SPIB/GZMB are decisive for plasmacytoid dendritic cells and the annotation follows the markers.
- Ten primary sources with DOIs/PMIDs triangulated across MPNST spatial transcriptomics, MPNST single-cell multiomics, NF1-PNST immunotyping, peripheral-nerve fibroblast biology and canonical immunology.

### SLURM

| JobID | State | ExitCode | Elapsed | CPUs | Mem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| **19886713** | **COMPLETED** | **0:0** | **00:08:35** | 8 | 128G | **8.55 GiB** (6.7%) |

No failures. R warnings: 0 (stderr holds only cosmetic ggplot colour-scale messages). All 7 round-trip checks TRUE.

### Figures

`results/phase2/figures/M15/` — broad/detailed/state/confidence UMAPs, labelled cluster UMAP, canonical marker dot plot and heatmap, eight per-compartment FeaturePlot panels, cluster→annotation map, confidence-by-compartment (PDF + PNG).

---

## M14 — Marker Discovery

- **Status**: Completed (2026-09-02). Reports: [M14_REPORT.md](reports/phase2/milestones/M14_REPORT.md), [MARKER_REPORT.md](reports/phase2/MARKER_REPORT.md)
- **Scope**: exploratory **cluster-characterisation** markers for annotation. **NOT condition-level differential expression**; cells are not independent biological replicates.

| Item | Value |
| --- | --- |
| Input | `results/phase2/clustering/phase2_harmony_clustered.rds` (md5 `b91f0eec…`) |
| Identities | `postint_harmony_primary_cluster` (resolution 1.0, 26 clusters); alternative 0.7 also run |
| Assay / layer | `SCT` / `data` — **`PrepSCTFindMarkers()` run** (4 SCT models, mandatory) |
| Parameters | `wilcox` (presto), `min.pct=0.25`, `logfc.threshold=0.25`, `only.pos=TRUE`, seed 42 — matching Phase 1 `config/config.yaml` |
| Harmony coordinates used for testing? | **No** — Harmony is embedding-level and carries no expression values |
| Results | 35,437 marker rows; 31,774 significant; every cluster returned markers (min 212, max 2,901) |
| Tables | 11 in `results/phase2/markers/` |
| Figures | 6 × PDF + PNG in `results/phase2/figures/M14/` |

### SLURM

| JobID | State | ExitCode | Elapsed | CPUs | Mem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| **19886699** | **COMPLETED** | **0:0** | **00:04:26** | 12 | 250G | **16.40 GiB** (6.6% efficiency) |

No failures or retries. R warnings: 0.

### Scientific findings

- Seurat v5 layer structure inspected explicitly: `RNA` layers are split four ways, `SCT` carries **4 models** → `PrepSCTFindMarkers()` mandatory and applied.
- **C9 retains MPZ** with GFRA3 and ABCB5 — strongest Schwann-lineage evidence in the dataset. **C8 expresses L1CAM** (reported SCP-like malignant marker). **C14 expresses SHH** with GAL3ST1 and KLK6.
- **C21 markers (LILRA4, CLEC4C, SPIB, GZMB) are decisive for pDC**, overriding the coarse M12 programme score that had called it B/plasma — marker evidence corrected a programme-score prior.
- **C15 is mitochondrial-dominated** (technical; matches the Phase 1 MPNST_4 flag) and **C23 ribosomal-pseudogene-dominated** — neither is a cell type.
- **C12 carries both a hypoxia programme and the reported perineurial GLUT1/ITGB4 signature** — genuinely ambiguous.

---

## M13 — Post-Harmony Neighbours, UMAP and Clustering Sweep

- **Status**: Completed (2026-09-02). Reports: [M13_REPORT.md](reports/phase2/milestones/M13_REPORT.md), [CLUSTERING_ASSESSMENT.md](reports/phase2/CLUSTERING_ASSESSMENT.md)
- **Recommendation**: **PRIMARY resolution 1.0 (26 clusters)** · **ALTERNATIVE resolution 0.7 (21 clusters)**

### Inputs / outputs

| Item | Value |
| --- | --- |
| Input | `results/phase2/harmony/phase2_harmony_integrated.rds` (md5 `cf63e843…`, verified) |
| Output object | `results/phase2/clustering/phase2_harmony_clustered.rds` — 5.64 GB, md5 `b91f0eecd73e202804ac7c6859db4c46` |
| Tables | 10 in `results/phase2/clustering/` |
| Figures | 22 figures × PDF + PNG in `results/phase2/figures/M13/` |

### Parameters

`FindNeighbors(reduction="postint_harmony", dims=1:30, k.param=20)` → graphs `postint_harmony_nn`/`postint_harmony_snn`; `RunUMAP(reduction="postint_harmony", dims=1:30, seed.use=42)` → `postint_umap_harmony`; `FindClusters(algorithm=1 Louvain, random.seed=42, resolution=0.1..1.0)`. All ten solutions preserved as `postint_harmony_clusters_res_*`. Cells 19,716 (in = out).

### SLURM

| JobID | State | ExitCode | Elapsed | MaxRSS | Note |
| --- | --- | --- | --- | ---: | --- |
| 19886683 | FAILED | 1:0 | 00:01:49 | 8.49 GiB | Preservation guard fired (whole-frame digest too blunt) |
| 19886685 | FAILED | 1:0 | 00:02:20 | 9.88 GiB | Unnamed-vector indexing defect |
| 19886687 | FAILED | 1:0 | 00:02:34 | 15.3 GiB | `sprintf("%d", median())` format defect |
| **19886690** | **COMPLETED** | **0:0** | **00:09:37** | **12.09 GiB** | Success (8 CPUs, 128G, 08:00:00 requested) |

All three failures were code defects fixed at source; **no resource was increased in response to any failure**. R warnings on the successful run: 0.

### Scientific decisions and findings

- 13 → 26 clusters over resolutions 0.1 → 1.0; **no tiny cluster at any resolution** (smallest 141 cells = 0.7%).
- **Every adjacent-resolution ARI ≥ 0.89** (0.892–0.976): higher resolution splits groups rather than reorganising them, so the resolution choice is not critical downstream.
- Primary selected by a documented composite rule (0.30 stability + 0.20 compactness + 0.15 compartment coverage + 0.15 program purity + 0.10 low fragmentation + 0.10 granularity); **UMAP appearance is not a criterion**.
- **Caveat**: resolution 1.0 sits at the sweep boundary, so its stability term is one-sided and inflated. It still leads on the two-sided criteria (highest silhouette 0.2078, full program coverage, no tiny clusters). Alternative 0.7 is interior and two-sided-validated.
- **11 of 26 primary clusters are >60% one sample** — expected, since M12 showed Harmony deliberately left the Schwann/neural-crest compartment patient-private.
- **M12 caveats carried forward as flags**: 3 clusters C1 (B/plasma >20%: C2, C21, C22), 6 clusters C2 (fibroblast >20%: C0, C3, C5, C7, C15, C25).
- Strongest Schwann-lineage/tumour-like candidates: C9 (713 cells, 79.8% MPNST_1, SchwannNC 0.986), C14 (440, 94.3%, 0.993), C8 (967, 91.0%, 0.799).

### Preservation

All Phase 1 / M11 state byte-identical after clustering and after a disk round-trip, including all 93 original metadata columns (value and R type). `FindClusters()` clobbers the legacy `seurat_clusters` column and active identities; both were snapshotted and restored verbatim and verified.

---

## M12 — Pre/Post Harmony Evaluation

- **Status**: Completed (2026-09-02). Reports: [reports/phase2/milestones/M12_REPORT.md](reports/phase2/milestones/M12_REPORT.md) and [reports/phase2/HARMONY_ASSESSMENT.md](reports/phase2/HARMONY_ASSESSMENT.md)
- **Harmony assessment status**: **`ACCEPT DEFAULT HARMONY WITH CAVEATS`**
- **Harmony re-run or re-tuned**: **No.** The M11 default result was evaluated as-is.

### Configuration

| Item | Value |
| --- | --- |
| Evaluated object | `results/phase2/harmony/phase2_harmony_integrated.rds` (md5 `cf63e843…`, sha256 `6085976f…`, both re-verified at runtime) |
| Cells evaluated | **19,716 — all. No subsampling anywhere in M12.** |
| PRE space | `pca` dims 1:30 (frozen Phase 1 M8 baseline) |
| POST space | `postint_harmony` dims 1:30 (M11) |
| Pre-Harmony UMAP | `umap_preintegration` — Phase 1, reused unmodified |
| Post-Harmony UMAP | `umap_harmony_m12` — identical `RunUMAP` call, only the input reduction changed (dims 1:30, seed 42, n.neighbors 30, min.dist 0.3, cosine — all Seurat defaults). Persisted as a 2-column embedding only; no third multi-GB object created. |
| k | 15 (primary, matches Phase 1) and 50 (robustness) |
| Methods | `RANN::nn2` exact kNN; `cluster::silhouette` on full pairwise distances |
| Unavailable packages | `lisi`, `kBET`, `clustree` — **not installed**; inverse Simpson computed natively (same estimator as `lisi`) and a dominance-ratio-vs-global-composition used in place of kBET |

### Major quantitative findings

**Pipeline validated against the frozen baseline.** M12 PRE-Harmony per-sample metrics reproduce Phase 1 M8 to 15 significant figures. Phase 1's headline `0.9613` is the *unweighted* mean of four per-sample means (cell-weighted = `0.9728`, the M12 value); Phase 1's entropy `0.0966` is in **log2** (= 0.0490 nats, the M12 value).

Because the samples are unequal in size, "fully mixed" means matching global composition: same-sample fraction **0.3065**, entropy **1.2683 nats**, inverse Simpson **3.2627** — not 0, log K and 4.

| Metric (k = 15, all cells) | PRE | POST | Fully-mixed ref. | Gap closed |
| --- | ---: | ---: | ---: | ---: |
| Same-sample neighbour fraction | 0.9728 | **0.7784** | 0.3065 | 29.2% |
| Neighbourhood entropy (nats) | 0.0490 | **0.3565** | 1.2683 | 25.2% |
| Inverse Simpson (iLISI-equiv.) | 1.0549 | **1.4408** | 3.2627 | 17.5% |
| Dominance ratio | 3.30 | **2.53** | 1.00 | — |

At k = 50: same-sample 0.9443 → 0.7053 (37.5% closed). Every sample improved — MPNST_2 most (0.910 → 0.504), MPNST_1 least (0.988 → 0.849).

**Mixing is compartment-specific — the decisive result.** Same-sample fraction PRE → POST: Panleukocyte 0.933 → **0.527**, Myeloid 0.958 → **0.611**, T/NK 0.953 → **0.689**, Endothelial 0.962 → **0.709**, Mural 0.951 → **0.724**, Fibroblast 0.990 → 0.839, **Schwann/neural-crest 0.990 → 0.900**, B/plasma 0.986 → 0.918. Harmony mixed the shared immune and vascular compartments hard while leaving the presumptive malignant Schwann-lineage compartment (86.5% MPNST_1) nearly as patient-private as it began — the desired behaviour, and the reason the *global* figure looks modest.

**A metric that failed.** Global silhouette of `sample_id` was −0.0104 **pre**-Harmony despite 97% same-sample neighbours, because each sample spans malignant, immune and stromal states so within- and between-sample mean distances nearly cancel. Batch structure here is local; only the kNN metrics detect it. Reported as a caveat, not as evidence.

### Biological-preservation findings

Legacy `orig.anno` was **deliberately excluded** as legacy integration-derived annotation. Two non-legacy references were used.

| Measure | Value |
| --- | ---: |
| Within-sample kNN retention (k = 15) | **0.8206** |
| Phase 1 independent-cluster coherence PRE → POST | 0.8644 → **0.8427** (−2.5% relative) |
| Global kNN retention | 0.6598 (the fall *is* the correction) |
| Canonical-program silhouette PRE → POST | 0.0627 → **0.0729** (+0.0103, improved) |

Six of ten canonical programs became *more* cohesive; largest gains T/NK +0.209, endothelial +0.218, Schwann/neural-crest +0.170, myeloid +0.104 — same-lineage cells from different patients were brought together. The canonical programs are a documented sanity check, **not** an annotation (26.3% of cells left `Unassigned`); they were never written into any saved object.

### Over-correction status — **present but bounded**

- **B/plasma cohesion fell 44%** (silhouette 0.437 → 0.243; n = 836, 93.5% MPNST_3). Same-sample fraction barely moved, so cells were not dispersed across samples — they became less separable from adjacent lymphoid/myeloid territory. Classic failure mode for a batch-private population with no counterpart to align to.
- **Fibroblast cohesion fell 44%** (0.174 → 0.097; n = 4,075, 70.6% MPNST_4), with the lowest within-sample retention of any program (0.795).
- **MPNST_4** shows the largest within-sample damage (retention 0.781, coherence −3.6%), consistent with the Phase 1 flag on its `percent.mt`-correlated cluster C09.
- Against widespread over-correction: overall program cohesion rose, within-sample retention 0.82, the 54 Phase 1 clusters remain discrete islands, and the Schwann/neural-crest compartment was neither dispersed nor collapsed.

### Under-correction status — **real in absolute terms, largely appropriate**

Only 29.2% of the achievable mixing gap closed; every dominance ratio remains above 2.2, highest for the two small samples (MPNST_2 4.35, MPNST_3 4.95) — uniform default `theta = 2` under a 3.3× size imbalance corrects minority samples proportionally less. But the residual sits in the compartments where patient-private structure is expected. The one genuine concern is the weaker relative correction of MPNST_2 and MPNST_3.

### Unresolved confounding

`sample_id` is simultaneously dataset, patient and the only batch proxy; no clinical or technical covariate exists. Residual structure cannot be apportioned between uncorrected batch effect and real between-patient tumour biology, and no metric can prove the removed variance was technical. The compartment-resolved result is a consistency argument, not proof. **The integrated embedding cannot support any between-tumour, between-patient, between-condition or differential-abundance claim.**

### Files created

- `scripts/R/phase2/evaluate_harmony.R`, `scripts/shell/phase2/run_m12_evaluate.sh`
- `reports/phase2/HARMONY_ASSESSMENT.md` (sections A–J), `reports/phase2/milestones/M12_REPORT.md`
- `results/phase2/harmony/evaluation/` — 12 machine-readable outputs (`pre_post_mixing_summary.tsv`, `neighborhood_mixing_metrics_by_sample.tsv`, `technical_silhouette_summary.tsv`, `biological_preservation_summary.tsv`, `sample_restricted_population_check.tsv`, `biological_program_gene_sets.tsv`, `program_by_sample_composition.tsv`, `harmony_evaluation_parameters.tsv`, `m12_headline_metrics.json`, `prov_m12_evaluation.json`, `umap_harmony_m12_embedding.{tsv,rds}`, `figure_index_m12.tsv`)
- `reports/phase2/figures/m12/` — 12 figures × PDF + PNG
- `logs/phase2/slurm/m12_evaluate_19886628.{out,err}`

### Files modified

- `reports/FIGURE_INDEX.tsv` — added `phase`, `milestone`, `slurm_job_id` columns; back-filled 260 Phase 1 rows; appended 24 M12 rows (284 total)
- `PROGRESS.md`, `CHANGELOG.md`

### Files deliberately not modified

`processed_mpnst.rds`, `results/datasets/**`, `results/combined/**`, `results/phase1_manifest.json`, and `results/phase2/harmony/phase2_harmony_integrated.rds` — the M11 object was opened read-only, its checksums re-verified, and its embedding digest re-checked after all computation.

### SLURM

| JobID | Node | CPUs | Mem | Walltime | State | ExitCode | Elapsed | MaxRSS | Mem eff. |
| --- | --- | ---: | ---: | --- | --- | --- | --- | ---: | ---: |
| **19886628** | ihc-grid-1-1-1 | 8 | 96G | 04:00:00 | **COMPLETED** | **0:0** | **00:04:09** | **8,338 M (8.14 GiB)** | 8.5% |

**R warnings: 0. R errors: none.** stderr holds one informational Seurat message about the `RunUMAP` backend (identical to the one Phase 1 used, preserving comparison fairness). Post-job inspection confirmed all 24 figures and 12 tables present and non-trivial, no NA/non-finite metric values, all four samples represented in every breakdown, no metric outside its valid range.

### Recommendation

**`ACCEPT DEFAULT HARMONY WITH CAVEATS`.** Both criteria met — mixing improved substantially and in the right compartments, and biological structure was preserved. Two caveats travel into M13–M17:

- **C1.** Any **B/plasma** cluster downstream may be distorted; cross-check against the preserved non-integrated baseline before annotating.
- **C2.** Same for **fibroblast/stromal** clusters and MPNST_4-derived structure.

A sensitivity analysis is **not** recommended and was not run: raising `theta` would push hardest on the already-fragile Schwann/neural-crest and B/plasma populations, and lowering it would undo the immune/vascular alignment that justifies integrating. If wanted anyway, the two candidate experiments are documented in `HARMONY_ASSESSMENT.md` §J; neither has been executed.

### M12 supplement (2026-09-02)

Under the continuous-execution authorization, the explicitly requested individually-named full-page pre/post panels were added in `results/phase2/figures/M12/` (`M12_01`–`M12_06`, `M12_09`–`M12_10`, PDF + PNG). `M12_07/08` (biological condition) were **deliberately not generated** — no condition field exists in the data; `M12_07_08_condition_figures_NOT_APPLICABLE.txt` records why. SLURM: 19886675 FAILED (`dpi = NA` rejected by ggplot2 4.0.1 — code defect, fixed at source), 19886682 COMPLETED 00:01:32, MaxRSS 7.18 GiB.

---

## M11 — Default Harmony Integration

- **Status**: Completed (2026-09-02). Full report: [reports/phase2/milestones/M11_REPORT.md](reports/phase2/milestones/M11_REPORT.md)
- **Scientific status**: **Harmony execution has completed, but integration quality has NOT yet been scientifically accepted. Pre/post-Harmony evaluation is reserved for M12.**

### Input

| Item | Value |
| --- | --- |
| Object | `results/combined/pre_integration/combined_preintegration.rds` |
| md5 / sha256 (both re-verified at runtime) | `88a442688f912d882f6c6da01820e329` / `c3fdce8b…dc66` |
| Cells | 19,716 · Default assay `SCT` (29,113 features) · `RNA` 31,764 features |
| PCA reduction consumed | `pca` (19,716 × 50 dims, assay `SCT`) |

### Harmony configuration

| Item | Value |
| --- | --- |
| Grouping variable | **`sample_id`** (4 levels: MPNST_1 7,615 · MPNST_2 2,284 · MPNST_3 2,940 · MPNST_4 6,877; identical to `orig.ident`) |
| Dimensions | **`1:30`** of `pca` |
| Output reduction | **`postint_harmony`** (key `postintharmony_`, assay `SCT`) |
| harmony version | **1.2.4** |
| Random seed | 42 |
| Scientific parameters | **all package defaults** — `theta` 2, `sigma` 0.1, `lambda` 1, `nclust` 100 (auto cap), `max_iter` 10, `early_stop` TRUE, `ncores` 1, `project.dim` TRUE, `harmony_options()` untouched |
| Explicit arguments | only `object`, `group.by.vars`, `reduction.use`, `dims.use`, `reduction.save`, `verbose` |

Exact call:
`harmony::RunHarmony(object = obj, group.by.vars = "sample_id", reduction.use = "pca", dims.use = 1:30, reduction.save = "postint_harmony", verbose = TRUE)`

Convergence: **converged after 9 of a maximum 10 iterations** (`early_stop` triggered); objective 753.089 → 353.656 over eight iterations, +0.35% at iteration 9. A matrix-level re-run with identical seed and defaults reproduced the embedding **bit-identically (max abs difference 0.000e+00)**, so the convergence history provably describes the primary run and the result is deterministic.

### Output

| Item | Value |
| --- | --- |
| Object | `results/phase2/harmony/phase2_harmony_integrated.rds` (6,048,637,630 bytes) |
| md5 | `cf63e84313a91de25fe2e41660f78f06` |
| sha256 | `6085976f788633b68c07c1ddff49e20e6f6f6727579cb7900854febd2ae6e69d` |
| Cells | 19,716 (in = out) |
| Reductions | `pca`, `umap_preintegration`, **`postint_harmony`** (19,716 × 30) |
| Graphs | `SCT_nn`, `SCT_snn` · Assays `RNA`, `SCT` (default `SCT`) · 93 metadata columns |

### SLURM

| JobID | Node | CPUs | Mem | Walltime | State | ExitCode | Elapsed | MaxRSS | Mem eff. |
| --- | --- | ---: | ---: | --- | --- | --- | --- | ---: | ---: |
| **19886486** | ihc-grid-1-1-1 | 8 | 64G | 04:00:00 | **COMPLETED** | **0:0** | **00:08:24** | **8,788,416 K (8.38 GiB)** | 13.1% |

Harmony itself took 10.1 s; the wall time is dominated by `saveRDS` (256.7 s) and checksumming (110 s) of a 5.6 GB object.

### Validation results — all TRUE

- **Phase 1 preservation proof**: `pca` embedding, `umap_preintegration` embedding, the full metadata frame, metadata column names, cell names, graphs, assays and default assay are all byte-identical before and after Harmony.
- **Round-trip validation**: the saved RDS was re-read from disk after freeing memory and re-verified — correct cell count, Phase 1 `pca`/`umap_preintegration` present and unchanged, `postint_harmony` present with 30 dims × 19,716 cells, serialised embedding identical to the in-memory one, no non-finite values, row names matching cells in order.
- **Pre-run guards**: output path new and outside every protected Phase 1 location; both input checksums matched; no duplicate barcodes; grouping variable matched the M10 structure exactly; requested dims available; target reduction name did not already exist.

### M10 Finding 2 — CLOSED

The zero-variance `percent.mt` regression for MPNST_1 produced **no** corruption: 0 non-finite values in `pca` (all 50 dims), 0 in `SCT@scale.data` (5,192 × 19,716, none in any sample), and 0 in the Harmony embedding. Incidental note for M14: `SCT@scale.data` spans **5,192** genes, not 3,000, because four SCT models contribute a union of variable features.

### Files created

- `scripts/R/phase2/run_harmony_integration.R` — M11 Harmony script with pre-run guards, preservation proof and round-trip validation; includes a `--validation-mode` flag used for synthetic smoke testing
- `scripts/shell/phase2/run_m11_harmony.sh` — SLURM launcher (8 CPUs, 64G, 04:00:00)
- `results/phase2/harmony/phase2_harmony_integrated.rds`
- `results/phase2/harmony/harmony_parameters.json`, `harmony_parameters.tsv`
- `results/phase2/harmony/harmony_convergence.tsv`, `harmony_kmeans_objective.tsv`
- `results/phase2/harmony/harmony_embedding_dimension_summary.tsv`, `harmony_grouping_composition.tsv`
- `results/phase2/harmony/prov_m11_harmony.json` — Phase 1-compatible provenance record
- `reports/phase2/milestones/M11_REPORT.md`
- `logs/phase2/slurm/m11_harmony_19886486.{out,err}`

### Files modified

- `PROGRESS.md`, `CHANGELOG.md`

### Files deliberately not modified

`processed_mpnst.rds`, `results/datasets/**`, `results/combined/**` (mtime still 2026-07-19 23:25), `results/phase1_manifest.json`, `config/config.yaml`, `workflow/Snakefile`, every Phase 1 report, and `reports/FIGURE_INDEX.tsv` (M11 generated no figures).

### Warnings and errors

- **R warnings: 0** (`options(warn = 1)` in force). **R errors: none.** ExitCode 0:0.
- `m11_harmony_19886486.err` contains only Harmony's own `message()` progress output for the primary and diagnostic runs, including `Harmony converged after 9 iterations`. No error or warning.

### Unresolved questions carried into M12

1. Perfect batch/biology confounding is now embedded in the corrected space; no between-tumour or differential-abundance claim can be supported from it.
2. Over-correction is untested. The three sample-restricted populations flagged at M10 (B cells ~95% MPNST_3, Malignant SCP-like ~97% MPNST_1, cycling tumour cells ~82% MPNST_4) must be checked explicitly in M12.
3. Group sizes are unbalanced 3.3× while default `theta = 2` is applied uniformly.
4. `nclust` defaulted to the package cap of 100 soft clusters.
5. Four SCT models remain; `PrepSCTFindMarkers()` stays mandatory for M14.
6. Phase 1 documentation defects D1–D9 remain uncorrected pending instruction.

### Recommendation

**Technically ready for M12.** A falling Harmony objective is an optimiser meeting its own criterion, not evidence of correct integration; the scientific assessment is M12's job.

### STOP

**STOPPED.** No clustering, neighbours, UMAP, markers, annotation or pre/post assessment were performed. M12 requires explicit researcher authorization.

---

## M10 — Phase 2 Reconstruction and Integration Planning

- **Status**: Completed (2026-09-02). Full report: [reports/phase2/milestones/M10_REPORT.md](reports/phase2/milestones/M10_REPORT.md)
- **Harmony executed**: **No.** M10 is planning and validation only.

### Files created

- `scripts/R/phase2/inspect_phase1_handoff.R` — read-only structural validation of the Phase 1 handoff object
- `scripts/shell/phase2/run_m10_inspect.sh` — SLURM launcher for the above
- `reports/phase2/PHASE2_PLAN.md` — modular Phase 2 plan, namespace contract, module map
- `reports/phase2/HARMONY_VARIABLE_DECISION.md` — Harmony grouping-variable decision (proposed, unapproved)
- `reports/phase2/milestones/M10_REPORT.md` — M10 milestone report
- `reports/phase2/PHASE2_ENVIRONMENT.tsv`, `reports/phase2/PHASE2_SESSIONINFO.txt` — environment records
- `results/phase2/handoff/` — 8 machine-readable inspection outputs (structure JSON, assays, reductions, metadata inventory, batch cross-tab, per-sample cells, per-sample QC, legacy annotation cross-tab)
- `logs/phase2/slurm/m10_inspect_19886411.{out,err}`
- Phase 2 directory skeleton under `results/phase2/`, `reports/phase2/`, `scripts/R/phase2/`, `scripts/shell/phase2/`, `logs/phase2/`, `benchmarks/phase2/`

### Files modified

- `PROGRESS.md`, `CHANGELOG.md`, `README.md` — Phase 2 status
- `scripts/python/validate_project_state.py` — the Phase 2 status check now requires an explicit `PHASE 2 = ...` declaration (NOT STARTED / IN PROGRESS / FROZEN / COMPLETE) instead of the hard-coded `NOT STARTED`, which became false at M10

### Files deliberately not modified

`processed_mpnst.rds`, `results/datasets/**`, `results/combined/**`,
`results/phase1_manifest.json`, `config/config.yaml`, `workflow/Snakefile`, and every
Phase 1 report. `reports/FIGURE_INDEX.tsv` is unchanged because M10 produced no figures.

### Jobs submitted

| JobID | Name | Partition/Node | CPUs | Mem | Walltime | State | ExitCode | Elapsed | MaxRSS |
| --- | --- | --- | ---: | ---: | --- | --- | --- | --- | ---: |
| 19886411 | `p2_m10_inspect` | ihc / ihc-grid-1-1-1 | 4 | 96G | 02:00:00 | COMPLETED | 0:0 | 00:02:17 | 4.73 GiB |

Memory efficiency 4.9% — over-provisioned. The M11 proposal (8 CPUs / 64G / 04:00:00) is
derived from this observation plus the Phase 1 M8 benchmark, not from the resource maximum.

### Scientific decisions

1. **Phase 2 input object**: `results/combined/pre_integration/combined_preintegration.rds`, verified against **both** Phase 1 checksums (md5 `88a4426…` from `prov_analysis.json`, sha256 `c3fdce8…` from `phase1_manifest.json`). Both match.
2. **Proposed Harmony grouping variable**: `sample_id` (4 levels), sole entry in `group.by.vars`. `orig.ident` is cell-for-cell identical and therefore redundant. Every other candidate is prohibited, degenerate, continuous, or non-existent.
3. **Proposed Harmony dimensions**: `1:30` of reduction `pca` — identical to the frozen Phase 1 baseline (91.1% cumulative variance), so M12 compares integration rather than dimensionality.
4. **Harmony parameters**: all scientific parameters left at `harmony 1.2.4` defaults (`theta` 2, `sigma` 0.1, `lambda` 1, `nclust` auto, `max_iter` 10, `early_stop` TRUE). Only technical arguments specified.
5. **Namespace**: Phase 1's mandated `postint_*` prefix reconciled with the Phase 2 specification names → `postint_harmony`, `postint_umap_harmony`, `postint_harmony_clusters_res_0.1`…`_res_1.0`, `postint_celltype_level1/2/3`, `postint_annotation_confidence`.
6. **Pre-Harmony baseline**: the existing `pca` and `umap_preintegration` reductions *are* the baseline and are preserved unmodified; no duplicate `pca_pre_harmony` copy is created.
7. **No Snakemake in Phase 2**: modular R scripts plus explicit `sbatch` scripts.
8. **Integration diagnostics computed natively** (inverse Simpson / entropy / same-batch fraction / silhouette / pre→post kNN retention) rather than installing `lisi` or `kBET`.
9. **Legacy `orig.anno` is provenance-only** and was not used to select any parameter.

### Warnings

- No R warnings or errors; the M10 job's stderr is empty.
- **Finding 1**: the Phase 1 `SCT` assay contains **four** SCTransform models (`model1`, `model1.1`, `model1.2`, `model1.3`), not the single global model described in `PRE_INTEGRATION_ASSESSMENT.md` §3.2 and `PHASE1_HANDOFF.md` §3.1. Seurat v5 ran SCTransform once per merged `RNA` layer. This does not invalidate integration, but `PrepSCTFindMarkers()` becomes mandatory in M14 (and `JoinLayers()` would be required if the `RNA` assay were used instead — its counts/data layers are split four ways).
- **Finding 2**: `percent.mt` is identically 0 in all 7,615 MPNST_1 cells, yet Phase 1 M8 passed `vars.to.regress = "percent.mt"` unconditionally, regressing a zero-variance covariate for that model. A non-finite-value integrity check is added to the front of M11.
- **Nine Phase 1 documentation defects (D1–D9)** are catalogued in `M10_REPORT.md` §7. None invalidates integration; all are prose-vs-record inconsistencies. Nothing was repaired. D8 (`config/config.yaml: input_rds` points at the deleted `/local/projects-t3/lilab/vmenon/Pilot_tumor/` root) is the only one with an operational consequence, and it does not affect Phase 2.

### Unresolved questions (researcher decision required)

1. Approve `sample_id` as the Harmony grouping variable, accepting that dataset = sample = patient = presumed technical batch is a single inseparable variable — and choose Option A / B / C in `HARMONY_VARIABLE_DECISION.md` §7.1.
2. Approve `dims.use = 1:30`.
3. Approve the `postint_*` naming contract.
4. May legacy `orig.anno` serve as a reference label for M12 biological-conservation diagnostics only? Default if unanswered: **no**.
5. Should Phase 1 documentation defects D1–D9 be corrected, now or at the Phase 2 freeze?
6. Should `lisi` / `kBET` / `clustree` be installed? Default: **no**.
7. Accept the four-SCT-model structure as-is with `PrepSCTFindMarkers()` in M14, rather than re-running M8? Recommended: **accept**.

---

# Phase 1 (frozen) — historical record


This document tracks the progress of the Phase 1 independent dataset analysis, reproducible workflow engineering, marker discovery, and pre-integration handoff.

| Milestone Status Summary

| Milestone | Description | Status | Target Date | Completion Date |
| --- | --- | --- | --- | --- |
| **M0** | Infrastructure, Environment, and SLURM Safety | **Completed** | 2026-07-06 | 2026-07-06 |
| **M1** | Real-Data Audit | **Completed** | 2026-07-07 | 2026-07-07 |
| **M2** | Dataset Extraction and Pre-Filter QC | **Completed** | 2026-07-08 | 2026-07-08 |
| **M3** | QC Filtering and Doublet Assessment | **Completed** | 2026-07-09 | 2026-07-09 |
| **M4** | Normalization and Variable Features | **Completed** | 2026-07-09 | 2026-07-09 |
| **M5** | PCA and PC Evaluation | **Completed** | 2026-07-10 | 2026-07-10 |
| **M6** | Clustering Resolution Sweep | **Completed** | 2026-07-18 | 2026-07-18 |
| **M7** | Marker Discovery and Dataset Recommendations | **Completed** | 2026-07-19 | 2026-07-19 |
| **M8** | Combined Pre-Integration Baseline | **Completed** | 2026-07-20 | 2026-07-20 |
| **M9** | Workflow Hardening, CI/CD, Provenance, and Phase 1 Freeze | **Completed** | 2026-07-20 | 2026-07-20 |

---

## Detailed Milestone Logs

### M0 — Infrastructure, Environment, and SLURM Safety
- **Status**: Completed (2026-07-06)
- **Conda Environment**: Successfully installed Snakemake without affecting critical packages. Generated a portable env specification `workflow/envs/R_env_portable.yaml` and a comparison diff.
- **Repository Setup**: Initialized directory skeleton. Moved verification checks to `scripts/shell/verify_milestone.sh`.
- **Configuration**: Created base config files and schema validation.
- **Synthetic Data**: Generated and verified synthetic data via Snakemake workflow execution from clean state.
- **SLURM Integration**: Submitted a trivial smoke job (JobID `19160128`) on `ihc` partition and verified successful execution log.
- **Whitespace / Quality Check**: Cleaned all trailing whitespaces. Verified that `git diff --cached --check` passes.
- **Git Commits**:
  - Initial M0 Implementation Commit: f6da0c537dc618925b6fed673a249e2601599d89
  - Documentation Commit: e1ab6d11c4234c498728db47c5b0e7962c6d5a94
  - M0 Cleanup and README Documentation Commit: 89ca4ea3aacae9d5a5da60f73c1a8e358fec7ae3

### M1 — Real-Data Audit
- **Status**: Completed & Extended to Complete Scientific Inventory (2026-07-07)
- **Codebase Update**: Developed `scripts/R/audit_object.R` to run read-only scientific audits, supporting both synthetic and split-layer Seurat v5 structures. Developed `scripts/R/inspect_additional.R` to retrieve advanced object statistics (memory size, duplicate genes/cells, and mitochondrial/ribosomal gene lists). Developed `tests/unit/test_audit.R` to validate audit outputs.
- **Synthetic Validation**: Ran synthetic data validation through local Snakemake execution. All unit tests passed successfully.
- **SLURM Production Run**: Submitted the real dataset audit via sbatch (JobID `19160528`) on partition `ihc` node `ihc-grid-1-1-1`. Run additional inspection via JobID `19160550`.
- **HPC Execution Metrics**:
  - State: COMPLETED (ExitCode 0:0)
  - Elapsed: 00:01:32 (SLURM job), 76.94 seconds (Snakemake rule)
  - MaxRSS: 9627268K (~9.18 GB)
  - CPU Usage: mean_load 54.72%, cpu_time 42.23s
  - Estimated object memory footprint: 10.08 GB
- **Key Scientific Findings**:
  - Object dimensions: 29,708 features x 22,661 cells in default SCT assay (RNA assay contains 31,764 features x 22,661 cells).
  - Raw counts: Verified split raw RNA counts exist in the Seurat object (split layers `counts.MPNST_1` through `counts.MPNST_4` under `RNA` assay). Reprocessing from raw UMI counts is scientifically feasible.
  - Dataset Identifier: Recommended using `sample_id` as the primary dataset splitter. Both `sample_id` and `orig.ident` have identical cardinality (4 samples: MPNST_1 (8338 cells), MPNST_2 (2830 cells), MPNST_3 (3682 cells), MPNST_4 (7811 cells)) and map 100% identically (0 mismatches).
  - Feature Statistics: 0 duplicated genes, 0 duplicated cell names. 13 mitochondrial genes matching `^MT-` and 433 ribosomal genes matching `^RP[SL]` identified.
  - Legacy Inventory: Successfully logged all legacy clusters (SCT, Harmony, CCA, RPCA, MNN) and dimensional reductions. Legacy columns have been marked for exclusion from Phase 1 processing to prevent data leakage.
- **Artifacts Generated / Expanded**:
  - [DATA_AUDIT.md](file://reports/DATA_AUDIT.md) (Expanded with full evidence tables and all 12 requested sections)
  - [DATASET_STRUCTURE.tsv](file://reports/DATASET_STRUCTURE.tsv)
  - [object_inventory.json](file://reports/object_inventory.json) (Expanded with all machine-readable metrics)
  - [object_inventory.tsv](file://reports/object_inventory.tsv) (Expanded with recommended usages)
  - [M1_REPORT.md](file://reports/milestones/M1_REPORT.md) (Revised milestone report)
  - [m1_audit_provenance.json](file://reports/audits/m1_audit_provenance.json)

### M2 — Dataset Extraction and Pre-Filter QC
- **Status**: Completed (2026-07-08)
- **Codebase Update**: Developed `scripts/R/extract_dataset.R` to split the Seurat object by `sample_id` and generate validation reports/manifests, and updated it with a robust counts checker. Developed `scripts/R/generate_qc_plots.R` to produce 4 types of publication-quality diagnostic plots (PDF and PNG) per sample and compile local figure indices. Developed `scripts/R/generate_qc_recommendation.R` to assess QC metrics against thresholds and write markdown reports. Added unit test `tests/unit/test_extraction.R` to validate extracted object integrity.
- **Workflow Integration**: Updated `workflow/Snakefile` with rules for dataset extraction, plotting, recommendations, figure index merging, and extraction validation tests. Enabled sample list configuration in `config.yaml` and `config.test.yaml`.
- **Synthetic Validation**: Successfully executed synthetic tests locally, verifying that all rules and tests run to completion and pass.
- **SLURM Production Run**: Submitted the production run to SLURM (JobID `19161116`) on partition `ihc` node `ihc-grid-1-1-1`.
- **HPC Execution Metrics**:
  - State: COMPLETED (ExitCode 0:0)
  - Elapsed: 00:06:35
  - MaxRSS: 66935856K (~63.83 GB) — Efficiency of 99.7% of the 64 GB limit.
- **Key Scientific Findings**:
  - Successfully extracted four datasets (`MPNST_1` to `MPNST_4`).
  - Discovered that `MPNST_1` has 0.00% mitochondrial transcripts across all 8,338 cells.
  - Calculated percent.ribo across all cells using `^RP[SL]`.
  - Assessed expected filtering impact of default thresholds (combined filter excludes 16.21% for MPNST_1, 8.02% for MPNST_2, 5.43% for MPNST_3, and 4.07% for MPNST_4).
- **Artifacts Generated**:
  - `results/datasets/MPNST_*/` (Split Seurat objects and extraction provenance logs)
  - `reports/datasets/MPNST_*/DATASET_VALIDATION.md`
  - `reports/datasets/MPNST_*/QC_RECOMMENDATION.md`
  - `reports/datasets/MPNST_*/manifest.json`
  - `reports/datasets/MPNST_*/` QC plots (violins, scatter, histograms, metrics in PDF and PNG) and metrics summary TSVs
  - `reports/FIGURE_INDEX.tsv`
  - `reports/milestones/M2_REPORT.md`

### M3 — QC Filtering and Doublet Assessment
- **Status**: Completed (2026-07-09)
- **Codebase Update**:
  - Developed `scripts/R/filter_and_detect_doublets.R` to run scDblFinder doublet detection and apply QC filters, supporting a custom output suffix to enable dual-strategy evaluation.
  - Developed `scripts/R/qc_optimization_review.R` to run threshold sensitivity simulations and grid intersection overlap plots.
  - Developed `scripts/R/compare_qc_strategies.R` to quantitatively and visually compare global vs dataset-specific strategies.
  - Created unit tests `tests/unit/test_filtering.R` and `tests/unit/test_filtering_specific.R` to check compliance for both strategies.
- **Workflow Integration**: Updated `workflow/Snakefile` with rules for both Global (Strategy A) and Dataset-Specific (Strategy B) filtering, optimization review, comparison analysis, unit tests, and global figure index merging. Added wildcard constraints to prevent filename ambiguities.
- **Synthetic Validation**: Successfully executed synthetic tests locally, verifying that all rules and tests run to completion and pass.
- **SLURM Production Runs**:
  - Global Run: Submitted to SLURM (JobID `19161422` on partition `ihc` node `ihc-grid-1-1-1`, COMPLETED, 3m 24s, MaxRSS 36.88 GB).
  - Specific & Comparison Run: Submitted to SLURM (JobID `19161426` on partition `ihc` node `ihc-grid-1-1-1`, COMPLETED, 9m 40s, MaxRSS 63.67 GB).
- **Key Scientific Findings**:
  - Flat global thresholds lead to a catastrophic cell loss of **46.9%** in `MPNST_2`, **48.3%** in `MPNST_3`, and **56.8%** in `MPNST_4`.
  - Under proposed dataset-specific thresholds, cell retention increases to **80.7%** (`MPNST_2`), **79.8%** (`MPNST_3`), and **88.0%** (`MPNST_4`), avoiding artificial truncation of biological expression profiles (such as ribosomal and mitochondrial fractions naturally elevated in sarcomas) while successfully removing doublets (~8.0-9.2% of cells) and low-complexity cells.
- **Artifacts Generated**:
  - `results/datasets/MPNST_*/MPNST_*_filtered.rds` & `MPNST_*_filtered_specific.rds` (Filtered Seurat objects for both strategies)
  - `reports/datasets/MPNST_*/FILTER_REPORT.md` & `FILTER_REPORT_SPECIFIC.md` (Dataset filter reports)
  - `reports/datasets/MPNST_*/manifest_filtered.json` & `manifest_filtered_specific.json`
  - `reports/datasets/MPNST_*/` post-filter QC plots (violins, scatter, density, histograms, doublet summaries, filtering summaries for both strategies)
  - `reports/qc_optimization/` plots and statistics (distribution comparisons, sensitivity curves, and grid intersection overlaps)
  - `reports/qc_comparison/` plots and tables (retention barplot, mt/ribo violin comparison, composite distributions, and strategy comparison summaries)
  - `reports/QC_OPTIMIZATION_REPORT.md` and `reports/QC_COMPARISON_REPORT.md`
  - `reports/FIGURE_INDEX.tsv` (Authoritative merged index of all 32 generated figure paths)
  - `reports/milestones/M3_REPORT.md`

### M4 — Normalization and Variable Features
- **Status**: Completed (2026-07-09)
- **Codebase Update**:
  - Developed `scripts/R/normalize_and_find_features.R` to run Seurat v5 SCTransform v2 or LogNormalize dynamically, handling zero-mitochondrial libraries (like `MPNST_1`) and setting parallelization global limits options for the `future` package.
  - Developed `scripts/python/generate_m4_report.py` to compile normalization runtimes and HVFs into a consolidated summary report.
  - Created unit test `tests/unit/test_normalization.R` to validate normalized objects and diagnostic outputs.
- **Workflow Integration**: Extended `workflow/Snakefile` with rules `normalize_and_find_features`, `test_normalization`, and `generate_m4_report` using dataset-specific QC inputs (`*_filtered_specific.rds`). Updated global figure index rules.
- **Synthetic Validation**: Successfully validated the normalization workflow on synthetic datasets locally.
- **SLURM Production Run**: Submitted the production run to SLURM (JobID `19174235`) on partition `ihc` node `ihc-grid-1-1-1`.
- **HPC Execution Metrics**:
  - State: COMPLETED (ExitCode 0:0)
  - Elapsed: 00:04:17
  - MaxRSS: 25287072K (~24.12 GB) — Highly efficient memory utilization under SLURM.
- **Key Scientific Findings**:
  - Standardized SCTransform v2 normalization decoupled sequencing depth covariance across MPNST constituent libraries.
  - `MPNST_1` had no mitochondrial transcripts, and the workflow dynamically bypassed mitochondrial regression, completing successfully.
  - Highly variable features (HVFs) across libraries successfully captured sarcoma-related biology, including extracellular matrix elements (collagens, APOD), cell-cycle markers, and macrophage-associated chemokine markers (CCL3, CCL4).
- **Artifacts Generated**:
  - `results/datasets/MPNST_*/MPNST_*_normalized.rds` (Normalized Seurat objects)
  - `reports/datasets/MPNST_*/NORM_REPORT.md` (Dataset-specific normalization reports)
  - `reports/datasets/MPNST_*/variable_features.tsv` (List of 3,000 highly variable features with metrics)
  - `reports/datasets/MPNST_*/` diagnostic plots (scatter, distribution, and top 6 expression violins in PDF and PNG)
  - `results/datasets/MPNST_*/normalization_provenance.json` (Digests and session specifications)
  - `reports/milestones/M4_REPORT.md` (Consolidated Milestone 4 report)

### M5 — PCA and PC Evaluation
- **Status**: Completed (2026-07-10)
- **Codebase Update**:
  - Developed `scripts/R/run_pca_and_evaluation.R` to execute PCA on normalized SCT assays, calculate variance explained relative to z-scored residuals, test cell score correlations with technical covariates, find PC selection ranges (conservative, recommended, maximum) using a geometric elbow knee-point detector, and save detailed results.
  - Developed `scripts/python/generate_m5_report.py` to aggregate dataset-specific metrics, generate comparison tables, and write `reports/PCA_RECOMMENDATIONS.tsv` and `reports/milestones/M5_REPORT.md`.
  - Created validation unit tests in `tests/unit/test_pca.R` to verify PCA reductions, coordinate dimensions, loadings, reports, and provenance.
- **Workflow Integration**: Extended `workflow/Snakefile` with rules `run_pca_and_evaluation`, `test_pca`, and `generate_m5_report`. Updated the figure index merging rule.
- **Synthetic Validation**: Validated local Snakemake execution on synthetic datasets, checking test parameter parsing.
- **SLURM Production Run**: Submitted workflow execution to SLURM (JobID `19176730`) on partition `ihc` node `ihc-grid-1-1-1`.
- **HPC Execution Metrics**:
  - State: COMPLETED (ExitCode 0:0)
  - Elapsed: 00:03:51
  - MaxRSS: 11410872K (~10.88 GB)
- **Key Scientific Findings**:
  - PCA geometric elbows identified optimal recommended PC cutoffs: `PC8` (MPNST_1), `PC6` (MPNST_2), `PC9` (MPNST_3), and `PC5` (MPNST_4).
  - SCTransform normalization decoupled sequencing depth covariance, resulting in low correlations (R < 0.2) in leading PCs.
  - Leading PCs are heavily dominated by biological programs: antigen presentation/immune response (CD74, HLA-DRA, HLA-DRB1) and extracellular matrix structure/remodeling (COL1A1, COL1A2, COL3A1, DCN, SFRP2), representing core MPNST biology.
  - JackStraw analysis was omitted due to lack of statistical validity on regularized negative binomial SCT z-scored residuals and high computational footprint.
- **Artifacts Generated**:
  - `results/datasets/MPNST_*/MPNST_*_pca.rds` (PCA-embedded Seurat objects)
  - `reports/datasets/MPNST_*/PCA_REPORT.md` (Dataset-specific PCA reports)
  - `reports/datasets/MPNST_*/top_loading_genes.tsv` (Top loading gene lists)
  - `reports/datasets/MPNST_*/pc_technical_correlations.tsv` (Correlation metrics)
  - `reports/datasets/MPNST_*/pca_variance_explained.tsv` (Raw variance scores)
  - `reports/datasets/MPNST_*/` diagnostic plots (elbow, cumulative variance, PC loadings, PC heatmaps, and technical correlations in PDF and PNG)
  - `reports/PCA_RECOMMENDATIONS.tsv` (Machine-readable recommendations)
  - `reports/milestones/M5_REPORT.md` (Consolidated Milestone 5 report)
  - `reports/FIGURE_INDEX.tsv` (Updated global figure index with all 40 PCA plots)

### M6 — Clustering Resolution Sweep
- **Status**: Completed (2026-07-18)
- **Codebase Update**:
  - Developed `scripts/R/run_clustering_sweep.R` to run SNN graph construction using dataset-specific PCA recommendations, execute a Louvain clustering resolution sweep from 0.1 to 1.0, run subsampling stability bootstrapping (5 rounds of 80% cells) with cell-order-aligned ARI calculations, calculate linear regression $R^2$ of technical covariates (nCount_RNA, nFeature_RNA, percent.mt, percent.ribo) with handling for zero-variance covariates (like `percent.mt` in `MPNST_1`), select recommended resolutions using a biologically-relevant trade-off heuristic, and generate diagnostic plots (umap grid, recommended umap, metrics vs resolution, stability/R2 vs resolution, and clustering transition tree).
  - Developed `scripts/python/generate_m6_report.py` to compile dataset sweep metrics into a consolidated cross-dataset recommendations summary.
  - Created validation unit tests in `tests/unit/test_clustering.R` to verify clustered objects, coordinate embeddings, and stats.
  - Created SLURM script `scripts/shell/run_m6_workflow.sh` to encapsulate HPC execution configurations.
- **Workflow Integration**: Extended `workflow/Snakefile` with rules `run_clustering_sweep`, `test_clustering`, and `generate_m6_report`. Integrated into the figure index merging rule.
- **SLURM Production Run**: Submitted workflow execution to SLURM (JobID `19399140`) on partition `ihc` node `ihc-grid-1-1-1`.
- **HPC Execution Metrics**:
  - State: COMPLETED (ExitCode 0:0)
  - Elapsed: 00:07:58
  - MaxRSS: 29593148K (~28.22 GB)
- **Key Scientific Findings**:
  - Graph construction and clustering sweep successfully executed using dataset-specific PCs.
  - Subsampling stability checks identified optimal biologically-relevant recommended resolutions: `0.6` (MPNST_1, 18 clusters, stability ARI = 0.919), `0.3` (MPNST_2, 9 clusters, stability ARI = 0.933), `0.6` (MPNST_3, 13 clusters, stability ARI = 0.910), and `0.7` (MPNST_4, 14 clusters, stability ARI = 0.740).
  - `MPNST_4` flagged key technical covariate correlation with mitochondrial percentage (`R2_percent_mt = 0.47` at recommended resolution 0.7), signifying potential MT bias that needs to be monitored in downstream analysis.
- **Artifacts Generated**:
  - `results/datasets/MPNST_*/MPNST_*_clustered.rds` (Clustered Seurat objects)
  - `reports/datasets/MPNST_*/CLUSTERING_REPORT.md` (Dataset-specific clustering reports)
  - `reports/datasets/MPNST_*/*.tsv` (Sweep stats, recommendation summaries, and stability metrics)
  - `reports/datasets/MPNST_*/` diagnostic plots (umap grid, recommended umap, metrics vs resolution, stability/R2 vs resolution, and clustering transition tree in PDF and PNG)
  - `reports/CLUSTERING_RECOMMENDATIONS.tsv` (Consolidated recommendations)
  - `reports/CLUSTERING_SWEEP_SUMMARY.tsv` (Consolidated sweep summary)
  - `reports/milestones/M6_REPORT.md` (Consolidated Milestone 6 report)
  - `reports/FIGURE_INDEX.tsv` (Updated global figure index with all 20 clustering sweep plots)

### M7 — Marker Discovery and Dataset Recommendations
- **Status**: Completed (2026-07-19)
- **HPC Execution Metrics**:
  - **Marker Discovery Sweep**: SLURM JobID `19399784` (COMPLETED, Elapsed: 11m 8s, MaxRSS: ~28.43 GB)
  - **Visualization Audit & Index Merge**: SLURM JobID `19399799` (COMPLETED, Elapsed: 2m 21s, MaxRSS: ~5.22 GB)
- **Key Scientific & Engineering Accomplishments**:
  - Executed Wilcoxon rank-sum marker discovery across 4 datasets × 10 resolutions = 40 combinations using the project-approved `PrepSCTFindMarkers` workflow.
  - Recommended resolutions confirmed to possess robust, high-quality marker support (median 1446 markers for MPNST_1 at 0.6, 816 markers for MPNST_2 at 0.3, 816 markers for MPNST_3 at 0.6, 794 markers for MPNST_4 at 0.7) with zero weak or small clusters.
  - Validated MPNST_4 mitochondrial bias concern at resolution 0.7 (correlation with `percent.mt` $R^2 = 0.47$) and proposed resolution 0.5 (reducing correlation to $R^2 = 0.24$ and merging stress-response clusters) as the primary alternative baseline.
  - Standardized all recommended-resolution figures (top 5 marker heatmap, dot plot, and a 2x3 UMAP FeaturePlot panel showing 5–6 representative marker programs) and generated local figure indices (`figure_index_m7.tsv`).
  - Hardened the Snakemake workflow by explicitly tracking all 6 visual outputs (PDF/PNG format for heatmap, dot plot, and representative FeaturePlots) in `visualize_markers` and introducing a strict DAG dependency on `reports/FIGURE_INDEX.tsv` in `test_markers` to eliminate parallel race conditions.
  - Expanded unit tests in `tests/unit/test_markers.R` to programmatically validate all 24 required visual files and their global index registration.
- **Artifacts Generated**:
  - `reports/datasets/MPNST_*/markers/resolution_*/markers_all.tsv`, `markers_filtered.tsv`, `top_markers.tsv`, and `marker_summary.tsv` for all 40 combinations.
  - `reports/datasets/MPNST_*/markers/resolution_{rec}/figures/` (Top 5 marker heatmap, dot plot, and 2x3 representative marker UMAP panel in PNG and PDF, total of 24 visual files).
  - `reports/datasets/MPNST_*/markers/resolution_{rec}/figure_index_m7.tsv` (Local M7 figure index tables).
  - `reports/datasets/MPNST_*/ANALYSIS_RECOMMENDATION.md` (Dataset-specific analysis recommendation reports).
  - `reports/milestones/M7_REPORT.md` (Consolidated Milestone 7 report).
  - `reports/FIGURE_INDEX.tsv` (Consolidated global figure index with all 12 recommended-resolution marker figures).

### M8 — Combined Pre-Integration Baseline
- **Status**: Completed (2026-07-20)
- **HPC Execution Metrics**:
  - **Pre-Integration Baseline Run**: SLURM JobID `19403199` (COMPLETED, Elapsed: 17m 17s, MaxRSS: ~32.00 GB)
- **Key Scientific & Engineering Accomplishments**:
  - Merged the four constituent clustered Seurat objects, verifying cell counts (expected 19,716, actual 19,716) and cell name uniqueness (0 duplicates).
  - Implemented namespaced clustering assignments (`preint_MPNST_{ds}_res_{resolution}`) for all 10 resolutions (0.1 to 1.0) and generated globally unique recommended/alternative cluster labels (e.g. `MPNST_1_C00` ... `MPNST_4_C13`).
  - Constructed a mathematically valid shared non-integrated expression space by running a unified SCTransform on the merged raw counts, regressing `percent.mt` and selecting the top 3,000 variable features.
  - Computed a new shared PCA and a new shared UMAP embedding (`umap_preintegration`) using PCs 1-30, showing complete spatial segregation of the four datasets.
  - Implemented a block-based nearest neighbor algorithm in pure R (requiring zero external package dependencies) to calculate dataset-mixing diagnostics: mean same-dataset neighbor fraction is extremely high (>98% across all datasets), confirming massive batch/patient-specific segregation.
  - Performed a clinical metadata audit and documented complete confounding of patient/sample identity with dataset identity.
  - Generated and saved 10 UMAP and PCA baseline visualizations (PDF and PNG, total of 20 files) and registered them in the global figure index (`reports/FIGURE_INDEX.tsv`).
  - Generated machine-readable metadata inventory (`metadata_inventory.tsv`) and dictionary (`metadata_dictionary.tsv`).
  - Created a test suite `tests/unit/test_preintegration.R` to programmatically validate M8 integrity.
- **Artifacts Generated**:
  - `results/combined/pre_integration/combined_preintegration.rds` (Combined pre-integration Seurat object, 3.35 GB)
  - `reports/combined/pre_integration/metadata_dictionary.tsv` (Metadata dictionary)
  - `reports/combined/pre_integration/metadata_inventory.tsv` (Metadata inventory)
  - `reports/combined/pre_integration/composition_dataset.tsv` and `composition_cluster.tsv` (Composition analysis TSVs)
  - `reports/combined/pre_integration/confounding_summary.tsv` (Confounding summary TSV)
  - `reports/combined/pre_integration/neighborhood_mixing_summary.tsv` and `pca_variance_explained.tsv` (Diagnostics TSVs)
  - `reports/combined/pre_integration/figure_index_m8.tsv` (Local M8 figure index table)
  - `reports/combined/pre_integration/pca/` and `umap/` directories (10 figures in PDF and PNG format, total 20 visual files)
  - `reports/PRE_INTEGRATION_ASSESSMENT.md` (Pre-integration baseline assessment report)
  - `reports/INTEGRATION_PREPARATION.md` (Integration preparation design report)
  - `reports/milestones/M8_REPORT.md` (Consolidated Milestone 8 report)
  - `reports/FIGURE_INDEX.tsv` (Consolidated global figure index with all 10 M8 figures)

### M9 — Workflow Hardening, CI/CD, Provenance, and Phase 1 Freeze
- **Status**: Completed (2026-07-20)
- **Workflow Hardening**:
  - Implemented top-level canonical target `phase1_complete` in `workflow/Snakefile` that runs all validation and manifest steps.
  - Hardened configuration validation against schemas, path portability inside analysis scripts, and rerun safety.
  - Implemented failure propagation tests to ensure corrupted inputs halt execution immediately.
- **CI/CD Pipeline**:
  - Configured GitHub Actions CI workflow in `.github/workflows/ci.yml` that performs linting, config validation, synthetic data generation, and clean-room synthetic workflow runs.
- **Static Checks & Audits**:
  - Implemented Python project-state consistency checker `validate_project_state.py` to prevent documentation drift.
  - Compiled detailed three-layer verification audit for Milestones 0-8 in `reports/audits/M0_M8_RECONCILIATION.md` and `M0_M8_RECONCILIATION.tsv`.
  - Built preflight checker `preflight_checker.py` to verify system requirements.
- **Phase 1 Handoff**:
  - Generated machine-readable contract `results/phase1_manifest.json` and human-readable final handoff report `reports/PHASE1_HANDOFF.md`.
- **Artifacts Generated**:
  - `.github/workflows/ci.yml` (GitHub Actions workflow file)
  - `reports/audits/M0_M8_RECONCILIATION.md` and `M0_M8_RECONCILIATION.tsv` (Verification matrices)
  - `scripts/python/preflight_checker.py` (Environment checker)
  - `scripts/python/validate_project_state.py` (Consistency checker)
  - `scripts/python/generate_phase1_manifest.py` (Manifest generator script)
  - `results/phase1_manifest.json` (Machine-readable handoff contract)
  - `reports/PHASE1_HANDOFF.md` (Human-readable handoff document)
  - `reports/milestones/M9_REPORT.md` (Milestone 9 execution report)

---

# PHASE 2 STATUS: COMPLETE

---

# Phase 3 — Tumour–Immune Cell–Cell Communication (IN PROGRESS, authorized 2026-09-03)

**PHASE 3 STATUS: COMPLETE** — frozen 2026-09-03. M18–M27 all complete. Phase 4 (spatial validation) not begun and requires separate authorization.

| Milestone | Description | Status | Date |
| --- | --- | --- | --- |
| **M18** | Phase 3 reconstruction and method audit | **Completed** | 2026-09-03 |
| **M19** | Sample-aware CCC input preparation | **Completed** | 2026-09-03 |
| **M20** | Multi-method CCC execution (LIANA · CellChat · CellPhoneDB) | **Completed** | 2026-09-03 |
| **M21** | CCC concordance | **Completed** | 2026-09-03 |
| **M22** | Biological prioritisation | **Completed** | 2026-09-03 |
| **M23** | Receiver-response / NicheNet | **Completed** | 2026-09-03 |
| **M24** | LochNESS adaptation | **Completed** | 2026-09-03 |
| **M25** | CCC + LochNESS integration | **Completed** | 2026-09-03 |
| **M26** | Robustness and replication | **Completed** | 2026-09-03 |
| **M27** | Final freeze and handoff | **Completed** | 2026-09-03 |


All of M10, M11, M12, M13, M14, M15, M15A, M16 and M17 are complete and validated. The final
object carries both the detailed literature-supported annotation and the CCC-oriented
`annotation_ccc` layer. Phase 3 was not begun and requires separate authorization.

### Phase 3 summary

**Objective**: how MPNST tumour cells interact with immune and other TME populations; which
findings are concordant across independent frameworks; whether receiver cells show the
expected transcriptional response.

**Input**: `results/phase2/phase2_final_object.rds` (md5 `153d5f6a…`), `annotation_ccc`.
**CCC object**: `results/phase3/ccc/ccc_input_object.rds` (md5 `aeb02c59…`), **14 populations,
17,327 cells, all populations and all 182 directed pairs evaluable in all 4 samples**.
**Expression basis**: RNA assay joined + LogNormalize (SCT has 4 models); Harmony coordinates
never used as expression.

**Methods**: LIANA 0.1.14 · CellChat 2.2.0.9001 · CellPhoneDB 5.0.1 · NicheNet 2.2.1.1 ·
LochNESS (official MMCA formulation, R). Raw scores never averaged; each LR framework
contributes a boolean flag. LIANA/CellPhoneDB partial non-independence declared.

**Concordance**: 502,846 keys, 36,486 supported — High 5,848 · Moderate 499 · Single-method
24,654 · **Discordant 5,485**. Testable by 3 frameworks: only 1,347, of which **763 (57%)
were supported by all three**.

**Prioritisation**: lexicographic tiers, no composite score. P1 1,587 · P2 1,435 · P3 163.
**547 tumour-centric interactions reach ≥3 evidence streams; 325 reach all 4.**

**Key findings**: tumour→myeloid APP→CD74 (5 receivers, 3 frameworks, 4/4 patients),
CD99→PILRA, ANXA1→FPR1, HLA-F→LILRB1/2 · tumour→vasculature VEGFA→KDR/FLT1/NRP1 · reciprocal
perivascular **Notch** circuit (tumour JAG1→pericyte NOTCH3; pericyte/endothelial
JAG1/JAG2/DLL4→tumour NOTCH2) · mixed lymphoid picture (class I retained and IL15 best explains
CD8/NK/CD4 programmes, but HLA-E→KLRC1 and HLA-F→LILRB are inhibitory) · NicheNet myeloid
receivers converge on **CSF1**, lymphoid on **IL15**.

**Negative results retained**: LochNESS found no receiver-state structure associated with the
APP context in any lineage (best descriptive p = 0.333 against a 0.167 floor) and is itself
unstable (representation ρ 0.02–0.42, LOSO ρ 0.31–0.50) · NicheNet did **not** corroborate APP
for macrophages · fibroblast receiver-response is a clear negative (best AUPR 0.022) ·
TGFB1→TGFBR only weakly supported · 53.6% of supported interactions rest on one patient.

**LochNESS implementation comparison**: R (MMCA) vs Python with equivalent inputs →
**Pearson 1.0000, max abs diff 0**; vs the perturb-seq formulation as written → 0.3764, because
that formulation omits same-sample exclusion.

**Robustness**: leave-one-patient-out retention 71–79%; no single patient dominates.

**SLURM**: 19894719, 19894761, 19895124, 19895459, 19894939, 19895159 (after 19895133 failed),
19895438 (after 19895344 failed), 19895738–41 (after 19895564 was cancelled for wall clock),
19895947 (after two truncated downloads), 19896062 (after two failures), 19896063, 19896064,
19896065. Peak memory 27.43 GiB (6% of the envelope). **No failure was addressed by increasing
resources, and no Phase 2 package was downgraded.**

**Environment**: Seurat/SeuratObject/harmony/Matrix unchanged from Phase 2. One recorded
change: igraph 2.2.1 → 2.1.4 (conda solver; Phase 2 outputs are frozen artefacts, unaffected).
CellPhoneDB isolated in `cpdb_env`.

**Deliverables**: `results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv` ·
`results/phase3/phase3_manifest.json` · `reports/phase3/PHASE3_HANDOFF.md` · 36 final figures ·
29 final tables · 10 reports · 10 milestone reports · FIGURE_INDEX.tsv at 578 rows.

**STOP**: Phase 4 (spatial validation) not begun; requires separate authorization.


---

# PHASE 3 STATUS: COMPLETE

---

# PHASE 4 — SCEVAN MALIGNANCY REFINEMENT, TUMOUR-STATE RESOLUTION, TARGETED CCC REASSESSMENT

*Authorized 2026-09-03. Milestones M28 → M35, continuous execution.*

**PHASE 4 STATUS: COMPLETE** — frozen 2026-09-03. M28–M35 all complete. Phase 5 not begun and
requires separate authorization.

**Central question**: Phase 2 defined the malignant compartment conservatively — `MPNST-Tumor`
= 3,420 cells (17.35%) — because the only alternative available was the prohibited
non-immune-equals-tumour inference. That leaves 5,064 `Fibroblast`, 1,231
`Candidate-Malignant-Unresolved`, 720 `Uncertain` and 435 `Pericyte-VSMC` cells whose malignant
status **marker expression cannot resolve**, and that uncertainty propagates into every
tumour-centric Phase 3 interaction. Phase 4 brings orthogonal evidence — inferred copy number
via SCEVAN — and then asks whether the Phase 3 biological architecture survives a better
tumour definition.

**Phase 4 does not invalidate Phase 2 or Phase 3.** Phase 2 annotation is historical
biological annotation; Phase 4 malignancy is additional orthogonal evidence. Phase 4 only adds
metadata fields. Details in `PROJECT.md` §P4.

| Milestone | Description | Status | Date |
| --- | --- | --- | --- |
| **M28** | Phase 4 reconstruction and SCEVAN feasibility | **Completed** | 2026-09-03 |
| **M29** | Per-patient SCEVAN execution | **Completed** | 2026-09-03 |
| **M30** | Malignancy integration and refined calls | **Completed** | 2026-09-03 |
| **M31** | Malignant-only tumour-state analysis | **Completed** | 2026-09-03 |
| **M32** | Targeted Phase 3 CCC sensitivity analysis | **Completed** | 2026-09-03 |
| **M33** | Tumour-state → TME model | **Completed** | 2026-09-03 |
| **M34** | Robustness | **Completed** | 2026-09-03 |
| **M35** | Phase 4 freeze and handoff | **Completed** | 2026-09-03 |

### M28 — Phase 4 Reconstruction and SCEVAN Feasibility

**Status** Completed · 2026-09-03 · **Next** M29

**Input** `results/phase2/phase2_final_object.rds`, md5 `153d5f6acc70f9c05aa48cabc4f4ac2d`
(verified equal to the Phase 2 manifest). Read-only.

**Output** per-sample sparse integer RNA counts in
`results/phase4/scevan/by_sample/*/counts_raw.rds` · `m28_feasibility_facts.json` ·
`phase4_cell_metadata.tsv.gz` (19,716 x 155) · frozen `phase2_harmony_embedding.rds` ·
`reports/phase4/SCEVAN_FEASIBILITY.md` · `reports/phase4/MALIGNANCY_DECISION_RULES.md` ·
`reports/phase4/environment/PHASE4_INSTALL_LOG.md`

**Scripts** `scripts/phase4/scevan/m28_feasibility.R`,
`scripts/shell/phase4/run_m28_feasibility.sh`,
`scripts/shell/phase4/install_umap_learn.sh`, `scripts/phase4/utils/phase4_plot_utils.R`

**JobIDs / resources / runtime / MaxRSS**

| JobID | State | Elapsed | ReqMem | MaxRSS |
| --- | --- | --- | --- | --- |
| 19896576 | FAILED | 00:02:15 | 96 G | 27.60 GiB |
| 19896623 | FAILED | 00:03:14 | 96 G | 27.60 GiB |
| 19896654 | FAILED | 00:00:11 | 16 G | 0 |
| 19896672 | COMPLETED | 00:05:39 | 16 G | 0.48 GiB |
| **19896711** | **COMPLETED** | **00:05:07** | 96 G | **26.90 GiB** |

No resource was raised in response to any failure.

**Figures** none (M28 is an inspection milestone). **Tables** the facts JSON and the per-sample
count audit within it.

**Key numbers** 19,716 cells · 31,764 RNA features · **28,340 (89.2%) map to SCEVAN's chr1-22
annotation** · gene identifiers are **symbols** · all four `counts` layers integer-verified ·
0 duplicate cell IDs · round-trip `setequal` TRUE · 7,448 high-confidence immune reference
candidates · **11,308 disputed cells never used as a reference** · smoke test 600 cells ->
363 tumour / 214 normal / 23 filtered / 5 subclones in 3.86 min.

**Errors** three, all root-caused: (1) SCEVAN 1.0.3 ignores `output_dir` in its read-back and
plot helpers -> each run now sets its own working directory, package unpatched; (2)
`subcloneAnalysisPipeline` -> `plotTSNE` needs Python umap-learn, and it runs *before* the
clone labels are written, so it is not optional -> umap-learn 0.5.12 in an isolated
`p4_umap_env` via `RETICULATE_PYTHON`, `R_env` verified unchanged; (3) our own `set -u` versus
conda's `qt-main` activation hook -> `set -o pipefail` without `-u`.

**Warnings** SCEVAN covers chr1-22 only, so X/Y events are not assessable; SCEVAN removes
cell-cycle genes **and all `HLA-*` genes**, so the CNV analysis is structurally blind to the
HLA-E/HLA-F loci Phase 3 highlighted.

**Scientific decisions** per-patient SCEVAN is primary (`sample_id` = patient = dataset) ·
primary run uses `norm_cell = NULL` so confident-normal detection sees no Phase 2 label ·
**`FIXED_NORMAL_CELLS = TRUE` prohibited in this project** because the source forces every
non-reference cell to malignant · Fibroblast/Pericyte-VSMC/Candidate-Malignant-Unresolved/
Uncertain never used as references, and Endothelial excluded so its call stays independent ·
**decision rules committed before any result was inspected** · SCT (four models) not used.

### M29 — Per-Patient SCEVAN Execution

**Status** Completed · 2026-09-03 · **Next** M30

**Input** `results/phase4/scevan/by_sample/*/counts_raw.rds` (M28, sparse integer RNA counts)
and `phase4_cell_metadata.tsv.gz`. The 6 GB Phase 2 object was never reloaded.

**Output** 8 runs (4 patients x primary/sensitivity), each with `classDf`, CNA matrices,
clonal and subclonal `.seg` profiles, onco-heatmaps and native CNA heatmaps under
`results/phase4/scevan/by_sample/<patient>/<run>/output/`; plus `m29_<patient>_summary.json`.

**Scripts** `scripts/phase4/scevan/m29_run_scevan.R`, `scripts/shell/phase4/run_m29_scevan.sh`

**JobIDs / resources / runtime / MaxRSS** — array 19896712, 8 CPUs / 64 G / 12 h per task

| JobID | Patient | State | Elapsed | MaxRSS |
| --- | --- | --- | --- | --- |
| 19896712_1 | MPNST_1 | COMPLETED | 00:34:20 | **60.87 GiB (95% of request)** |
| 19896712_2 | MPNST_2 | COMPLETED | 00:08:13 | 6.24 GiB |
| 19896712_3 | MPNST_3 | COMPLETED | 00:11:20 | 18.10 GiB |
| 19896712_4 | MPNST_4 | COMPLETED | 00:42:07 | 50.49 GiB |

No failures, no retries, no resource increases. All four concurrent: 32 CPUs, 256 G of 450 G.

**Key numbers** 18,449 of 19,716 cells assessed (1,267 SCEVAN-filtered = **not assessed**,
never "non-malignant") · **8,794 SCEVAN-malignant** · 22 subclones (7/4/3/8) ·
primary/sensitivity agreement **0.9951 / 0.9663 / 0.0864 / 0.9974**.

**Malignant fraction of assessed, by Phase 2 population** — Fibroblast 0.628 / 0.559 / 0.078 /
**1.000**; Candidate-Malignant-Unresolved 0.733 / 0.167 / 0.789 / 0.429; MPNST-Tumor 0.225 /
0.089 / 0.429 / 0.912; Pericyte-VSMC ≤0.081 except MPNST_4 0.942; T/NK/B/myeloid 0.000 in
MPNST_1 and MPNST_2.

**Errors** none fatal. Every run emitted a caught `plotCloneTree` error (ggtree calls
`ggplot2:::is.waive()`, removed in ggplot2 4.x); SCEVAN catches it internally, so only the clone
**phylogeny plot** is missing and all clone assignments, CN profiles, CNA matrices and
onco-heatmaps are intact. **ggplot2 was deliberately not downgraded** — that would destabilise
the Phase 3 figure suite. Limitation K.

**Warnings** MPNST_1 peaked at 95% of its 64 G request; the M28-derived estimate was closer
than intended and a larger sample would need more. **MPNST_3 failed the immune sanity check**:
its primary run calls the lymphoid compartment malignant and the myeloid compartment normal
while the sensitivity run inverts it (agreement 0.0864, only 25 confident normal cells, lowest
depth of the four at median 1,594 genes/cell); 6 of 9 immune populations individually exceed
25% malignant, and its three "clones" are pure immune lineages. **MPNST_4's high overall immune
rate (0.293) is driven almost entirely by Plasma-cell (98.2%)** — an immunoglobulin-locus
artefact — with 0 of 9 populations failing and the highest agreement of the four (0.9974).

**Scientific decisions** per-patient runs (§14) · primary run non-circular (`norm_cell = NULL`)
· `FIXED_NORMAL_CELLS = TRUE` prohibited project-wide · both reference strategies run so
reference-sensitivity is measurable, which is the only reason MPNST_3's failure was detectable
· disagreement with Phase 2 preserved in both directions · clone labels patient-scoped, and
cross-patient similarity described as a shared/recurrent CNA pattern, never the same clone.

### M30 — Malignancy Integration

**Status** Completed · 2026-09-03 · **Next** M32 (then M31/M33/M34)

**Input** the eight M29 classifications, four CNA matrices, frozen Phase 2/3 metadata, and
`MALIGNANCY_DECISION_RULES.md` (thresholds fixed a priori; Amendments A1/A2 post-hoc and
labelled as such).

**Output** `PHASE4_MALIGNANCY_CALLS.tsv` (19,716 x 37, per-cell rule and printed reason) ·
`SCEVAN_VS_PHASE2_ANNOTATION.tsv` · `AMBIGUOUS_POPULATION_MALIGNANCY_AUDIT.tsv` ·
`SCEVAN_SAMPLE_RELIABILITY.tsv` · `MALIGNANCY_THRESHOLD_SENSITIVITY.tsv` ·
`MALIGNANCY_SUMMARY_BY_{CLUSTER,SAMPLE}.tsv` · `SCEVAN_{CLONES,CNV_SUMMARY,ARM_LEVEL_EVENTS,RECURRENT_BROAD_EVENTS,RUN_ACCOUNTING}.tsv`
· `reports/phase4/MALIGNANCY_REFINEMENT_REPORT.md` ·
`reports/phase4/AMBIGUOUS_POPULATION_MALIGNANCY_REPORT.md`

**Scripts** `m29_aggregate.R`, `m30_integrate_malignancy.R`, `m30_figures.R`, `run_m30.sh`

**JobIDs** 19897096 COMPLETED 00:01:06 / 0.57 GiB (superseded by Amendment A2) ·
**19897148 COMPLETED 00:01:00 / 1.38 GiB** (reported) · 19897097 CANCELLED with the chain ·
19897098–19897111 cancelled without running (two duplicate dependency chains fired
simultaneously by a background waiter and a monitor; both stopped, no scientific output).
No failures, no resource increases; 64 G was over-provisioned at 1.38 GiB used.

**Figures** all eight §27 figures plus three M29 aggregation figures, PDF + PNG.

**Key numbers** Phase 2 **3,420 (17.35%)** -> SCEVAN **8,794 of 18,449 assessed (47.7%)** ->
**refined 6,434 (32.63%)**, of which **3,261 High confidence**. Non-malignant 9,078 ·
**Ambiguous 3,766 (19.10%)** · Excluded 438. **R99 fail-safe = 0 cells**, so rule coverage is
complete. Fibroblast **4,036 / 5,064 Malignant** · Candidate-Malignant-Unresolved **836 / 1,231
Malignant, 0 Non-malignant** · MPNST-Tumor **1,405 / 3,420 retained, 2,015 -> Ambiguous** ·
Pericyte-VSMC 348 / 435 Non-malignant · Uncertain 476 / 720 Non-malignant. Recurrent broad
CNA events in the three reliable patients: **chr18 loss 3/3, chr2 gain 3/3, chr7 gain 3/3**.

**Warnings** the `MPNST-Tumor` retention is the most fragile number in Phase 4 (MPNST_1's
pop_frac 0.222 sits just below the a priori 0.25; a 0.20 threshold would retain 3,266 rather
than 1,405) whereas the **fibroblast finding is completely threshold-independent** (4,036 at
every grid point 0.15-0.40). `Fibroblast` falls to 1 cell in MPNST_4 and `Pericyte-VSMC` to 4,
so those become NOT EVALUABLE there.

**Scientific decisions** Amendments A1 and A2 (both post-hoc, both dated, A2 explicitly noting
it was prompted by the flat gate excluding the strongest result, with both gates reported and a
`flat_gate_would_have_failed` column) · rule **R7p** for the plasma-cell immunoglobulin-locus
artefact · marker evidence may only **reduce** confidence, never raise it, to avoid
reintroducing circularity · `Ambiguous` is terminal and never becomes tumour · threshold
sensitivity published rather than a single a priori value.

### M31 — Malignant-Only Tumour-State Analysis

**Status** Completed · 2026-09-03 · **Next** M33

**Input** `malignancy_refined == "Malignant"`, **6,434 cells** (3,261 High / 3,173 Moderate),
built from the M28 RNA counts + LogNormalize; the 6 GB Phase 2 object was not reloaded and the
four-model SCT assay was not used.

**Output** `results/phase4/tumor_states/malignant_only_object.rds` ·
`TUMOR_STATE_{ASSIGNMENTS,MARKERS,PATIENT_DISTRIBUTION,DEFINITIONS,PROGRAM_SCORES,RESOLUTION_SWEEP}.tsv`
· `CLONE_VS_TUMOR_STATE.tsv` · `m31_tumor_state_facts.json` · 8 figures (`31_01`–`31_08`)

**Scripts** `scripts/phase4/tumor_states/m31_tumor_states.R`,
`scripts/shell/phase4/run_m31_m33_m34.sh`

**JobIDs** 19897617 COMPLETED 00:02:24 / 5.16 GiB (**superseded** — recurrence labelling
violated §36) · **19897639 COMPLETED 00:02:22 / 5.41 GiB** (reported). No resource increases.

**New reductions** `malignant_pca` → `malignant_harmony` (batch `sample_id`) → `malignant_umap`,
dims 1:30 capturing **93.0%** of variance. Phase 2's `pca`/`postint_harmony`/
`postint_umap_harmony` were never modified, and M35 asserts every Phase 2 embedding is
numerically identical in the final object. The new reduction was **justified by measurement**:
kNN(k=20) overlap with the Phase 2 Harmony embedding is only **0.1733**.

**Resolution sweep** 0.2→8, 0.3→10, 0.4→10, 0.5→11, 0.7→15, 1.0→17 clusters. No resolution
achieved full multi-patient coverage, so the printed fallback rule selected **0.2** (8 states).

**Key numbers — 8 states**: Mesenchymal_ECM-1 1,851 · ECM-2 1,598 · ECM-3 1,291 · ECM-4 648 ·
ECM-5 505 · Cycling 290 · Schwann_like 179 · Interferon 72. No cluster needed the `Uncertain`
label. Clone counts per state 2–19, so a state routinely spans many CNV clones — clone and
state are not equivalent (§33).

**HEADLINE NEGATIVE RESULT (§36)** `recurrent (>=3 patients each >=5%) = 0 of 8` ·
`shared between 2 patients = 1` (Schwann_like) · `patient-specific or dominated = 7` ·
**97.2% of malignant cells sit in patient-private states**. Expected with n = 4 and
`sample_id` = patient = dataset, but reported as a negative result about tumour-state generality.

**Coherent biological finding** Phase 2 `MPNST-Tumor` cells map to Cycling (248) /
Schwann_like (175) / Interferon (72); the promoted fibroblasts map to the five Mesenchymal_ECM
states (1,245 + 1,128 + 637 + 512 + 481). **Phase 2 detected the marker-legible malignant cells
and missed the Mes-NC-like ECM cells, because an ECM programme is what a fibroblast looks like**
— which is why copy-number evidence rather than a better marker panel was the right instrument.

**Errors** one substantive error of mine, caught before it propagated: the first run classified
recurrence on patient *presence*, labelling `Mesenchymal_ECM-2` "recurrent across ≥3 patients"
when 1,587 of 1,598 cells were MPNST_4 and the others contributed 1, 3 and 7. §36 forbids this.
The rule was rewritten to require no patient above 80% **and** ≥3 patients each contributing
≥5% and ≥10 cells, and the script now prints the negative result explicitly. M33 was halted
mid-run and both were re-executed.

**Scientific decisions** subset on `malignancy_refined`, not the broad Phase 2 label · new
reductions only, Phase 2's untouched and asserted so · state count not predetermined and the
selection rule printed before application · programme names earned by a stated ≥0.02 margin with
`Uncertain` available · recurrence judged by patient share, not presence · clone and state kept
separate.

### M32 — Targeted Phase 3 CCC Sensitivity Analysis

**Status** Completed · 2026-09-03 · **Next** M33/M34

**Input** `PHASE4_MALIGNANCY_CALLS.tsv` + the frozen Phase 2 object → refined CCC object
`ccc_input_object_refined.rds` (**15,036 cells, 14 populations**, md5 `99fd641f…`). Phase 3
baseline `MPNST_CCC_MASTER_TABLE.tsv` and `CCC_CONCORDANCE.tsv`.

**Method reuse** `scripts/phase3/ccc/{run_liana.R,run_cellchat.R,run_cellphonedb.py}` and
`scripts/phase3/concordance/build_concordance.R` reused **completely unmodified** at Phase 3
versions (liana 0.1.14, CellChat 2.2.0.9001, CellPhoneDB 5.0.1 in `cpdb_env`), same resource,
thresholds, min-cells 10, seed 42 and expression basis. Only `ccc_label` differs — which is
what makes this a label-sensitivity analysis rather than a tool-drift confound. Nothing was
reinstalled.

**Output** `PHASE3_VS_PHASE4_CCC.tsv` (506,605 rows) · `PHASE3_AXIS_SENSITIVITY.tsv` ·
`FIBROBLAST_INTERACTION_REASSIGNMENT.tsv` · `CCC_LABEL_PHASE3_VS_PHASE4.tsv` ·
`CCC_REFINED_SAMPLE_CELL_COUNTS.tsv` · refined `CCC_CONCORDANCE.tsv` · figures `32_01`, `32_02`
and the `build_concordance` suite in `results/phase4/figures/M32/`

**Scripts** `prepare_refined_ccc_inputs.R`, `m32_ccc_sensitivity.R`, `run_m32_{prepare,liana,cellchat,cellphonedb,concordance,sensitivity}.sh`

**JobIDs / resources / runtime / MaxRSS**

| JobID | Stage | State | Elapsed | ReqMem | MaxRSS |
| --- | --- | --- | --- | --- | --- |
| 19897175 | refined CCC input | COMPLETED | 00:07:04 | 96 G | 27.78 GiB |
| 19897176 | LIANA | COMPLETED | 00:34:36 | 64 G | 38.87 GiB |
| 19897177_1 | CellChat MPNST_1 | COMPLETED | 00:39:48 | 48 G | 18.33 GiB |
| 19897177_2 | CellChat MPNST_2 | COMPLETED | 00:20:28 | 48 G | 14.49 GiB |
| 19897177_3 | CellChat MPNST_3 | COMPLETED | 00:19:40 | 48 G | 14.58 GiB |
| 19897177_4 | CellChat MPNST_4 | COMPLETED | 00:23:32 | 48 G | 18.60 GiB |
| 19897178 | CellPhoneDB | COMPLETED | 00:02:31 | 48 G | 1.64 GiB |
| **19897179** | concordance + sensitivity | **FAILED** | 00:03:39 | 48 G | 1.34 GiB |
| 19897616 | sensitivity rerun | COMPLETED | ~00:01 | 48 G | — |

**Errors** 19897179's *concordance* step exited 0 (outputs intact); only my sensitivity script
failed, from two defects: `fib3` derived from `t3` had no `rk` column so `NULL %in% reasg`
returned `logical(0)`; and `canon()` stripping punctuation collapsed an upstream resource's
`"CD8 RECEPTOR"` and `"CD8_RECEPTOR"` into one key, leaving 4 duplicate keys per table and a
many-to-many join. Fixed by computing `rk` explicitly and collapsing duplicate keys while
keeping the strongest support **and logging how many were affected**, rather than silencing the
warning. Only the failed script was rerun — the 96-minute CCC computation was not repeated. No
resource increases. CellChat was parallelised per patient from the outset, applying the Phase 3
lesson.

**Key numbers** supported **36,486 → 34,893 (−4.4%)** · **High concordance 5,848 → 5,698
(−2.6%)** · Moderate 499 → 513 · Single-method 24,654 → 23,302 · Discordant 5,485 → 5,380 ·
testable-by-three 1,347 → 1,297. Change classes over 506,605 union keys: Stable 18,276 ·
Weakened 8,446 · Newly-supported 5,421 · Lost 3,963 · **Sender-reassigned 3,051** ·
Strengthened 2,746. **Median change across the 23 named axes −3.1%; no axis lost.**

**Named-axis verdicts** APP→CD74 −3.1% (stable) · ANXA1→FPR1 0.0% · HLA-F→LILRB1 −1.4% ·
HLA-F→LILRB2 −3.1% · HLA-E→KLRC1 −2.0% · CD58→CD2 −1.4% · JAG1→NOTCH3 0.0% · JAG2→NOTCH2 0.0% ·
FGF2→FGFR1 −4.3% · **strengthened** FN1→ITGAV_ITGB8 +33.3%, CD99→PILRA +14.0%, VEGFA→KDR
+11.1%, JAG1→NOTCH2 +7.5% · **weakened >20%** COL1A1→ITGAV_ITGB8 −44.4%, SLIT2→ROBO1 −38.9%,
DLL4→NOTCH2 −31.6%, COL6A2→ITGAV_ITGB8 −30.0%, VEGFA→FLT1 −28.0%.

**The fibroblast question (§43)** of 9,483 Phase 3 supported interactions involving
`Fibroblast`: **6,547 (69.0%) retained**, **1,565 (16.5%) reassigned to the refined tumour
compartment**, 1,371 (14.5%) lost or no longer testable. **All four collagen/FN1 → ITGAV_ITGB8
axes changed sender** — the ECM→integrin axis Phase 3 read as stroma-to-tumour is substantially
tumour-autocrine.

**Warnings** refinement emptied populations in some patients (`Fibroblast` 1 cell in MPNST_4,
`Pericyte-VSMC` 4), so evaluable populations fell to 10 in MPNST_3 and 12 in MPNST_4. The change
classes therefore distinguish `Ambiguous` (no longer testable) from `Lost` (testable, no longer
supported) — absence of evidence is not evidence of absence.

### M33 — Tumour-State → TME Model

**Status** Completed · 2026-09-03 · **Next** M34

**Output** `TUMOR_STATE_TME_EVIDENCE.tsv` (176 rows) · `TUMOR_STATE_PROGRAMME_LEADERS.tsv` ·
`TUMOR_STATE_LIGAND_EXPRESSION.tsv.gz` · `RECEIVER_RECEPTOR_EXPRESSION.tsv.gz` · figures
`33_01`–`33_03` · `reports/phase4/TUMOR_STATE_COMMUNICATION_REPORT.md`

**Scripts** `scripts/phase4/tumor_states/m33_state_tme_model.R`

**JobID** 19897639 (chained with M31/M34), COMPLETED 00:02:22, MaxRSS 5.41 GiB

**HEADLINE NEGATIVE RESULT** the tumour-state → TME model **cannot be established in this
cohort**, and the §68 structure was therefore NOT imposed on the data. Of 176 evidence rows,
**12 reach ≥3 patients, 34 reach ≥2, 121 rest on one patient, 21 on none** — and **all 12 of the
≥3-patient rows belong to `Cycling`**, a 290-cell state that is itself 83% one patient.
`Cycling` leads all six programmes purely because it is the only state present in enough
patients to be evaluated: an evaluability artefact, not a biological preference. The five
`Mesenchymal_ECM` states hold **91.3%** of the malignant compartment and have **exactly one
evaluable patient each**.

**Design** deliberately targeted (§45): LR frameworks were **not** run per state. Evidence is
per-state, per-patient expression of the prioritized ligands (**40/40 present**) and receptors
(39 present) at the Phase 3 10% detection floor and 10-cell stratum minimum, joint only when
both sides pass **in the same patient**; 440 sender rows over 11 evaluable strata, 1,950
receiver rows. **No composite score invented.** Phase 3 NicheNet and LochNESS reused verbatim.

**§46 preserved** NicheNet is not forced to support APP–CD74; the myeloid ligand programme and
the receiver-state cytokine programme (CSF1/IL15/TGFB1) are carried as separate streams and
never merged. LochNESS enters only as receiver-lineage context and remains negative.

### M34 — Robustness

**Status** Completed · 2026-09-03 · **Next** M35

**Output** `MALIGNANCY_LEAVE_ONE_PATIENT_OUT.tsv` · `MALIGNANCY_PATIENT_CONSISTENCY.tsv` ·
`CNV_PATTERN_PATIENT_SIMILARITY.tsv` · `CCC_CHANGE_PATIENT_DEPENDENCE.tsv` ·
`CCC_REFINED_LEAVE_ONE_PATIENT_OUT.tsv` · figures `34_01`–`34_04` ·
`reports/phase4/ROBUSTNESS_REPORT.md`

**Scripts** `scripts/phase4/malignancy/m34_robustness.R` · **JobID** 19897639, MaxRSS 5.41 GiB

**§49 answers** states multi-patient? **No** — 0/8 recurrent, 97.2% patient-private. CNV
patterns patient-specific? **Predominantly** — pairwise Jaccard 0.20–0.23 among reliable
patients, only 1 of 59 broad events in all four, shared core chr18 loss / chr2 gain / chr7 gain
each 3/3. CCC changes driven by one patient? **The `Newly-supported` class is — 91.8%
single-patient**, only 27 of 5,421 at ≥3 patients, whereas `Strengthened` is 10.3%
single-patient with 1,031 at ≥3. Core interactions LOSO-stable? **Yes**, retention
0.727 / 0.825 / 0.916 / 0.955.

**MOST IMPORTANT CAVEAT (Q5)** dropping MPNST_4 gives a refined malignant fraction of **0.214
against a Phase 2 fraction of 0.215** — fold change 1.00. **The headline 17.35% → 32.63% is NOT
patient-robust**: outside MPNST_4 the ~1,149 promoted fibroblasts and ~1,933 demoted
`MPNST-Tumor` cells roughly cancel. The *fibroblast* conclusion is a different claim and **is**
robust (Q6: majority-malignant in 3 of 4 patients).

**Q6** `Fibroblast` is the only disputed population majority-malignant in a majority of
patients (0.62 / 0.55 / 0.98; 0.08 in the sanity-failed sample).
`Candidate-Malignant-Unresolved` rests on MPNST_1 alone, `Pericyte-VSMC` and `Uncertain` on
MPNST_4 alone.

### M35 — Phase 4 Freeze and Handoff

**Status** Completed · 2026-09-03 · **Next** none; Phase 5 requires separate authorization

**Input** every M28–M34 artefact plus the frozen Phase 2 object (md5 re-verified).

**Output** `results/phase4/phase4_final_object.rds` (5.64 GB, md5
`e85ba8486e456917e2483f2773bdbaf3`, sha256 `a658447744e5ee8621fc5a3127492ceb11bce37bac4d6f3f0f7255bdcb647ae4`)
· `phase4_final_object_validation.json` · `results/phase4/phase4_manifest.json` (**35 sections**,
`figures_missing: []`, `tables_missing: []`, 11 limitations) · `reports/phase4/PHASE4_HANDOFF.md`
(27 sections) · **86 files in `results/phase4/figures/final/`** (16/16 required + 27 supporting,
PDF + PNG) · **47 tables in `results/phase4/tables/final/`** (all 15 §53-required present) ·
`reports/FIGURE_INDEX.tsv` 578 → **664 rows**

**Scripts** `m35_build_final_object.R`, `m35_cna_heatmap.R`, `m35_assemble.R`, `run_m35.sh`

**JobID / resources / runtime / MaxRSS** **19897640 COMPLETED 00:12:16, ReqMem 96 G, MaxRSS
13.38 GiB.** Three stages all exit 0. No failures, no resource increases.

**Validation 23/23 passed.** Preservation asserted column by column before saving: all 149
pre-existing metadata columns byte-identical · reductions unchanged · **every Phase 2 embedding
numerically identical (PCA and Harmony included)** · assays and cell order unchanged ·
`annotation_ccc_phase3` a verbatim copy of the frozen Phase 3 layer · no non-Malignant cell
carries the refined `MPNST-Tumor` label · tumour states only on Malignant cells. After reload:
0 duplicate IDs, 0 missing Phase 2 columns, 0 missing Phase 4 fields, **0 non-finite values
across 1,656,144 embedding values**, save/load md5 stable.

**Figures** `04_scevan_cna_heatmap` was built as a real ComplexHeatmap figure (genome-ordered
chr1-22, four annotation tracks, seeded 2,000-cell subsample per patient, in-figure banner on
sanity-failed samples) rather than substituting the profile-line figure for it.

**Warnings** three carried into the freeze: the 32.63% refined fraction is **not
patient-robust**; **no malignant state is recurrent** (0 of 8), so the tumour-state → TME model
could not be built and the §68 structure was not imposed; and **Amendment A2 was written after
observing that the flat gate excluded the sample carrying the strongest evidence** — both gates
published, direction robust, magnitude not.

**Scientific decisions** no new exploratory biology (§56); validate and freeze only. Phase 4
only ever added metadata, and that is asserted rather than assumed.

---

## M35A — SCEVAN Figure Consolidation and Evidence Visualization

*2026-09-03. Visualization milestone appended after the Phase 4 freeze. **Nothing scientific was
re-run or changed**: SCEVAN, the malignancy calls, the decision rules, amendments A1/A2, the clone
assignments, the tumour states, the CCC results and every threshold were re-used exactly as frozen.
`results/phase4/phase4_final_object.rds` was neither loaded nor modified (md5
`e85ba8486e456917e2483f2773bdbaf3` unchanged).*

**Why.** The SCEVAN CNA evidence existed but was not visible. 267 native SCEVAN files sat behind
`SCEVAN_NATIVE_FIGURE_INDEX.tsv`, and no final figure showed a native CNA heatmap — so the claim
that carries Phase 4 was defensible in the tables and only partly visible in the figures.

**Input** `results/phase4/malignancy/phase4_malignancy_metadata.rds` (5 MB, md5
`694b49bca0a52b0ea00266544276abaa`, 19,716 × 190) · the frozen Phase 4 TSVs · the stored native
SCEVAN PNG and `.RData` outputs. **The 6 GB final object was never loaded.**

**Output** 18 files in `results/phase4/figures/final/` — `17_scevan_native_cna_MPNST1` ·
`18_..._MPNST2` · `19_..._MPNST4` · `20_scevan_cna_reliable_patients` ·
`21_scevan_clone_composition_phase2` · `22_fibroblast_cna_burden_by_patient` ·
`23_fibroblast_malignant_vs_nonmalignant_cna_profile` ·
`24_phase2_to_phase4_malignancy_transition` · `25_MPNST3_scevan_failure_qc` ·
`26_fibroblast_malignancy_threshold_robustness` · `27_broad_cna_recurrence` ·
`28_phase4_scevan_evidence_summary` (PNG alongside for 21, 22, 24, 25, 26, 28). Three new tables:
`SCEVAN_CLONE_COMPOSITION_VISUALIZATION.tsv` (176 rows), `M35A_FIBROBLAST_CNA_EFFECT_SIZES.tsv`,
`M35A_FIBROBLAST_CNA_PROFILE_CORRELATION.tsv`. Reports
`reports/phase4/M35A_SCEVAN_FIGURE_AUDIT.md` and `PHASE4_HANDOFF.md` §28.
`FIGURE_INDEX.tsv` 664 → **682 rows**; `phase4_manifest.json` 35 → **36 sections**.

**Scripts** `scripts/phase4/figures/m35a_{common,native_figures,clone_composition,fibroblast_cna,transition,mpnst3_qc,threshold,broad_cna,summary,finalize}.R`
· `scripts/shell/phase4/run_m35a_figures.sh`, `run_m35a_finalize.sh`

**JobIDs / resources** 19899301, 19899307, 19899308, 19899311, **19899312** (figures; each
COMPLETED, ~00:04:25, ReqMem 24 G, peak MaxRSS **2.97 GiB**) and 19899313, **19899315** (finalize;
COMPLETED 00:00:09, ReqMem 8 G). The figure reruns are quality iterations found by inspecting the
rendered output, not failures. **No plotting problem was solved by increasing memory.**

**Native evidence surfaced, never redrawn.** The native SCEVAN CNA heatmaps are embedded verbatim
as rasters — the CNA matrix, its cell ordering and its subclone track are SCEVAN's own output — with
the Phase 2 / clone / malignancy annotation they cannot carry supplied as an aligned companion
panel built from the same frozen clone assignments. **No CNA value was redrawn.** The cytoband
onco-heatmaps were deliberately not surfaced (unreadable at page scale, and they invite gene-level
over-reading). **All 8 `CloneTree.png` files are blank** — the documented ggtree/ggplot2 4.x defect;
ggplot2 was not downgraded and phylogeny plotting was not forced.

**Validation** 18 frozen-value assertions run at the head of every figure script, all passed:
19,716 cells · Fibroblast 5,064 → 4,036 / 908 / 120 · Candidate-Malignant-Unresolved 836 / 395 / 0 ·
MPNST-Tumor 1,405 / 2,015 · refined 6,434 / 9,078 / 3,766 / 438 · clones 7 / 4 / 3 / 8 · agreement
0.9951 / 0.9663 / 0.0864 / 0.9974.

**New quantitative results** (all from stored output, none newly computed by any analysis method):
Cliff's delta for Fibroblast Malignant vs Non-malignant CNA metrics is **0.834–0.917 in MPNST_1** and
**0.781–0.890 in MPNST_2**; genome-wide CNA profile correlation r(Fib→Mal, MPNST-Tumor→Mal) =
**0.934** in MPNST_1 and **0.972** in MPNST_4 against r(Fib→Mal, Fib→Non-mal) = **0.184** in MPNST_1.

**Two things reported honestly rather than smoothed.** (1) The malignant-vs-non-malignant fibroblast
contrast **exists in 2 of 4 patients only** — MPNST_4 retains a single non-malignant fibroblast and
MPNST_3 contributes no malignant fibroblasts — and every figure prints NOT EVALUABLE where that is
the case. (2) **A wording discrepancy was found in `PHASE4_HANDOFF.md` §12**: "Every MPNST_1 clone
mixes Candidate-Malignant-Unresolved, Fibroblast and MPNST-Tumor" holds for **6 of 7** clones;
`MPNST_1_clone6` (176 cells) carries no Fibroblast cells. No call, count or conclusion changes; §12
is left as written and the correction is recorded in §28.7.

**A defect caught by not trusting an exit code.** The first finalize run (19899313) exited 0 and
silently rounded every frozen decimal in `phase4_manifest.json` (`jsonlite::toJSON` defaults to
`digits = 4`; agreement `0.99514117 → 0.9951`). It was caught by diffing against a pre-write backup.
The manifest was restored and regenerated with `digits = NA`, and `m35a_finalize.R` now asserts,
after writing, that all 33 pre-existing sections are byte-identical to a pre-write snapshot.

---

# PHASE 5 — MALIGNANT TRANSCRIPTIONAL PROGRAMS

*Opened 2026-09-03. Phases 1–4 remain complete and frozen and are not modified. Phase 5 reads
`results/phase4/phase4_final_object.rds` (md5 `e85ba8486e456917e2483f2773bdbaf3`, sha256
`a658447744e5ee8621fc5a3127492ceb11bce37bac4d6f3f0f7255bdcb647ae4`) read-only and writes only
under `results/phase5/`, `reports/phase5/`, `scripts/phase5/` and `logs/phase5/`.*

**Scientific question.** Phase 4 reported a negative result: 0 of 8 discrete malignant
transcriptional states were recurrent across ≥3 patients and 97.2% of malignant cells sat in
patient-private states. Phase 5 does **not** re-run clustering to rescue that. It asks whether
recurrent **continuous** transcriptional programs exist across patients even when discrete
clusters do not, and what separates malignant ECM-like MPNST cells from genuine fibroblasts.

**Milestones.** M36 reconstruction/feasibility · M37 continuous program discovery (cNMF) ·
M38 program annotation, pathways and relationship to the Phase 4 discrete states · M39 malignant
ECM-like cells vs true fibroblasts · M40 robustness · M41 freeze.

**Architecture decisions taken before production work.** Program discovery uses the 6,434
`Malignant` cells only, on RNA expression (never Harmony/UMAP/PCA/SCEVAN-CNA); `Ambiguous` cells
are projected afterwards and never reclassified; the historical annotation field is **verified**
against the frozen fibroblast split rather than assumed, and `annotation_ccc_refined` is barred
from that role because it encodes the Phase 4 conclusion; patient imbalance is handled with a
mandatory patient-balanced sensitivity analysis; recurrence criteria are declared before biology
is inspected; cell cycle is retained as a candidate program with a declared exclusion sensitivity
run. `R_env` is not modified — cNMF lives in an isolated `p5_cnmf_env` and gene sets are
downloaded as versioned GMT files with checksums.

---

# PHASE 6 — CLONAL, REGULATORY AND PLASTICITY ARCHITECTURE

*Planned; begins only after M41 completes successfully.*

**Scientific question.** Is malignant transcriptional heterogeneity primarily associated with
distinct CNA-defined clones, with substantial program diversity within clones, or a mixture — and
which TF and pathway activities distinguish the Phase 5 programs?

**Milestones.** M42 reconstruction · M43 clone↔program coupling · M44 within-clone program
diversity · M45 regulatory architecture (decoupleR/CollecTRI) · M46 pathway activity
(PROGENy/Hallmark) · M47 broad CNA → transcriptional consequences · M48 integrated architecture ·
M49 robustness · M50 freeze.

**Constraints taken before production work.** MPNST_3's SCEVAN clone structure is unreliable and
is excluded from all clone-based inference (QC panels only). Clone labels are patient-scoped and
never treated as homologous across patients. Within-clone diversity is described as
*transcriptional-program diversity*, interpretable as *consistent with plasticity* only after
technical confounders are assessed — never as observed switching or a trajectory. CNA→expression
associations are an internal consistency analysis, **not** orthogonal validation, because SCEVAN
infers copy number from expression in the first place.

---

## M36 — Phase 5 Reconstruction and Feasibility

*2026-09-03 · COMPLETE · SLURM 19899335 (00:02:05, 4 CPUs, 96 G, MaxRSS 29.28 GiB)*

**Question** Can the frozen Phase 4 malignant compartment support a search for continuous programs,
and which frozen field carries the pre-Phase-4 annotation? **Input** `phase4_final_object.rds`,
verified by md5 AND sha256, opened read-only. **Output** `PHASE5_FEASIBILITY.md`, five working
artefacts, patient-distribution and state-by-patient tables.

13 frozen Phase 4 counts re-asserted, all passed. The historical annotation field is
**`annotation_ccc_phase3`, VERIFIED** by reproducing the frozen fibroblast split 5,064 →
4,036/908/120 — `annotation_ccc_refined` was excluded by design because it encodes the Phase 4
conclusion. **Patient imbalance is 17.4×** (MPNST_4 3,685 = 57.3%; MPNST_3 212 = 3.3%) and
confounded with depth. 178 malignant cells carry no clone label; MPNST_3 has **0 High-confidence**
malignant cells. **Warnings** both facts constrain every later per-patient conclusion.
**Next** M37.

## M37 — Continuous Program Discovery (cNMF)

*2026-09-03 · COMPLETE · SLURM 19899333, 19899339, 19899353, 19899354*

**Question** Do continuous programs exist in the 6,434 malignant cells, and at what rank?
**Method** cNMF 1.7.1 in the **isolated `p5_cnmf_env`**; `R_env` untouched. Declared before any
factor was seen: malignant cells only, RNA raw counts, genes in ≥ 0.5% of malignant cells with
`^MT-` removed (19,663), K 4–15, 100 replicates, seed 42, dt 0.10. Cell-cycle, ECM, HLA, interferon,
Schwann, neural-crest and angiogenesis genes all retained.

Three runs declared in advance (`primary`, `balanced`, `nocc`); three more added later as declared
sensitivity analyses (`loo_×4`, `highconf`, `nortech`). Balanced design: cap 651 cells (MPNST_2, the
smallest patient with ≥ 500), cutting the dominant fraction 0.573 → 0.301.

**K-selection rule declared in code before it ran** — largest K with max program cosine ≤ 0.75, no
dead program, stability ≥ the median of qualifying K. Eligible {4,5,6,7,8}; **K = 8**.
**Resources** MaxRSS **1.92 GiB**; the cost is 1,200 NMF fits per run, not RAM. **Next** M38.

## M38 — Program Annotation, Pathways, Phase 4 State Comparison

*2026-09-03 · COMPLETE · SLURM 19899354, re-run 19899356*

**8 programs. Three are technical-dominated and are named for it**: P1 and P5
`Translation_ribosomal` (54% ribosomal / 66% pseudogene in their top 50) and P8
`Myeloid_ambient_like` (18 myeloid markers). The first labelling pass called P1 `Mesenchymal_ECM` on
a 0.20 family overlap while its top 15 genes were ribosomal proteins; an explicit technical-content
check was added and M38 re-run. **Only the naming changed** — not the factorization, K, or the
recurrence rule. The five biological programs are P2 Neuronal, P3 **Mesenchymal_ECM**, P4
Hypoxia_Angio, P6 Schwann_like, P7 **Cycling**.

**Recurrence, criteria declared before any label: 0 recurrent · 1 shared-limited · 7
patient-private**, and the one shared-limited program is technical. **Not one biologically
interpretable program is carried by even two patients at the declared threshold.** At a relaxed 0.10
activity threshold, **P3 (ECM) reaches three patients** — the one biological program that comes close.

**The five Mesenchymal_ECM states are not one program**: they map to three. Within MPNST_4, three
separately-clustered ECM states (3,383 cells) collapse onto one; across patients they do not.
**Every program is recovered in the balanced run** (cosine 0.50–0.98). **Next** M39.

## M39 — Malignant ECM-like Cells versus True Fibroblasts

*2026-09-03 · COMPLETE · SLURM 19899354*

Groups from the verified historical field: A 4,036 / B 908 / C 120 (projection only). **Evaluability
checked before any test**: only MPNST_1 (512 vs 303) and MPNST_2 (637 vs 503) qualify — MPNST_4 has
**1** non-malignant fibroblast and MPNST_3 **0** malignant ones. **A four-patient paired test was not
constructed.**

**Cross-patient correlation of the pseudobulk log2FC is only 0.091** — stated first, because it
bounds the whole comparison. **92-gene signature** under declared criteria applied in every evaluable
patient, with **no classifier and no train/test split**: 61 up in malignant ECM-like (CA12, IGFBP3,
COL11A1, COL14A1, GJA1, SPOCK1 …), 31 up in true fibroblast (**CDH19, APOD, SCN7A, ABCA6/8/9/10,
VIT, SPARCL1** — markers of nerve-associated / endoneurial stroma). The 120 Ambiguous fibroblasts
were projected for description and **not reclassified**. **Next** M40.

## M40 — Phase 5 Robustness

*2026-09-03/04 · COMPLETE · SLURM 19899355, 19899358, 19899359, 19899742*

Eight perturbations; two added after M38 (technical covariates, technical-gene-free universe).
Programs are stable to rank (median cosine 0.998 at K ± 1), patient balance (0.940), cell-cycle
removal (**1.000**) and confidence filtering, and no program's usage tracks a within-patient
technical covariate above **ρ = 0.37**.

**The decisive sensitivity result**: a fourth full cNMF run with ribosomal, pseudogene/lncRNA and
canonical myeloid genes removed, across the **whole K grid**, gives **0 recurrent programs at every
K from 5 to 15**. The patient-private result is **not** an artefact of those gene classes.
**Every patient-private program vanishes when its own patient is withheld** (LOO 0.17–0.35).

Two verdicts are published because the strict one is uninformative alone: `robust_overall` is FALSE
for all eight (five fail only the test a patient-private program cannot pass), while
`robust_excluding_patient_scope` is TRUE for P2, P3, P4 and P7. **Next** M41.

## M41 — Phase 5 Freeze

*2026-09-03/04 · COMPLETE · SLURM 19899359, 19899742, 19899744, 19899747, 19899750*

**Object** `results/phase5/phase5_final_object.rds`, 6,058,471,006 bytes, md5
`839e5157bc7c3470ddf86746c2e719e1`, sha256 `fe99ecf51154046145cf21f9c6960664d10205b90736bb13c61cccf409c40f00`,
19,716 cells, 196 metadata columns (183 Phase 1–4 + 13 Phase 5). **10 preservation guards + 8 reload
validations, all passed**, asserted column by column. Cells outside the malignant compartment carry
**projected** scores flagged `program_score_source == "projected"`.

**12/12 figures** (PDF + PNG) · **16 tables** · manifest 30 sections · `FIGURE_INDEX.tsv` 682 → **706**.
`R_env_PRE_PHASE5.yml` and `R_env_POST_PHASE5.yml` byte-identical.

**Failures, both root-caused and fixed at source, both preserved.** *19899359* — the figures script
read UMAP coordinates from `meta.data`, but they live in the object's **reduction**; fixed by reading
the frozen Phase 4 embedding table. The object build in that job had already passed all validations
and was not repeated. *19899742* — a new column named `sig` shadowed the `sig` data frame inside the
same `mutate()`; renamed `is_sig`. **Visual QC** found and fixed three further figure defects.

---

## M42 — Phase 6 Reconstruction

*2026-09-03 · COMPLETE · SLURM 19899748 (00:02:42, 8 CPUs, 128 G, MaxRSS 21.92 GiB)*

Phase 5 object verified against its manifest by md5 and sha256. **The MPNST_3 exclusion was
re-derived** from the frozen `scevan_sample_reliable` flag rather than taken on trust — it marks
MPNST_3 and only MPNST_3. 18 of 19 clones clear the 20-cell minimum; `MPNST_4_clone8` (5 cells) is
reported NOT EVALUABLE. **Next** M43.

## M43 — Clone ↔ Program Coupling

*2026-09-03 · COMPLETE · SLURM 19899748*

Within-patient only; permutation null shuffles clone labels **inside** the patient (1,000 perms).
**Median η² 0.050–0.072 — roughly 94% of each program's variance sits WITHIN clones.** 7 of 24
program × patient pairs are clone-associated. 21 of 24 reach p < 0.001, which is why effect size and
not the p-value carries the conclusion. **A patient's dominant program is active in essentially
every one of its clones.** **Next** M44.

## M44 — Within-Clone Program Diversity

*2026-09-03 · COMPLETE · SLURM 19899748*

Measured on both continuous scores and a hard assignment whose 0.10 margin is printed. **Diversity
is real and strongly patient-specific**: MPNST_1 clones carry 1.69–2.83 effective programs; MPNST_2's
four clones carry **exactly 1.00**; MPNST_4's carry 1.05–1.23. **Between-clone divergence is
0.1–3.4% of within-clone dispersion** — a patient's clones are transcriptionally near-interchangeable.

**Warning reported, not buried**: program dispersion tracks the fraction of High-confidence cells in
a clone at **ρ = 0.72**, above the declared disqualifying bar. The effective-number metrics stay
below it (0.56–0.62). Clone size and depth are not drivers. **Next** M45.

## M45 — Regulatory Architecture

*2026-09-03/04 · COMPLETE · SLURM 19899749*

decoupleR `run_ulm` over CollecTRI (41,674 edges, 1,201 TFs) on the frozen RNA log-normalised layer.
**P7 Cycling → E2F4 0.43\*, MYC 0.34\*, E2F1 0.24\*** and **P5 Translation → HIF1A 0.45\*, ATF4
0.43\*, HSF1 0.45\*** — textbook, and recovered without being imposed.

**A dependency defect worked around without moving a frozen package**: `decoupleR::get_collectri()`
and `OmnipathR::collectri()` both error in OmnipathR 3.14.0 (`unnest_evidences`). **OmnipathR was
not upgraded**; the identical data was taken from OmniPath's documented REST endpoint and cached with
provenance. **Caution stated**: with n = 4, three-of-four sign agreement occurs ≈ 31% of the time by
chance, so the concordance count is a weak filter and the named TFs rest on direction and magnitude.
**Next** M46.

## M46 — Pathway Activity

*2026-09-04 · COMPLETE · SLURM 19899749*

PROGENy (14 pathways) and Hallmark (50 sets), **side by side, never merged**. The pathway layer
independently supports labels derived from an entirely separate source: **P4 `Hypoxia_Angio` → PROGENy
Hypoxia 0.47\*** · **P7 `Cycling` → Hallmark E2F_TARGETS 0.43\*, MYC_TARGETS_V1 0.44\*** · **P3
`Mesenchymal_ECM` → Hallmark EMT 0.33\*, PROGENy EGFR 0.57\***. **Next** M47.

## M47 — Broad CNA → Transcriptional Consequences

*2026-09-04 · COMPLETE · SLURM 19899749*

Within-patient carrier vs non-carrier clone comparison on the mean expression of all genes in each
broad segment. **299 events tested, 263 evaluable, direction matches in 180 (68%)**; largest effects
are MPNST_1's chr8 gains at Cohen's d 2.65–3.23. **This is internal consistency, NOT validation** —
SCEVAN inferred those events from expression. NF1/NF2 segments are tabulated with the permitted and
prohibited wording printed next to each row. **Next** M48.

## M48 — Integrated Architecture

*2026-09-04 · COMPLETE · SLURM 19899749*

Criteria specified before the evidence was assembled. **Model A fails** on all three criteria, most
decisively because between-clone divergence is **0.7%** of within-clone dispersion. **Model B fails**
because only 22% of clones span ≥ 2 programs. **SELECTED: Model C, mixed architecture** — a minority
of programs associate detectably with clones while ~94% of program variance is within-clone, and a
patient's clones are transcriptionally near-interchangeable. MPNST_3's programs are described; **no
clone-based interpretation was invented for it**. **Next** M49.

## M49 — Phase 6 Robustness

*2026-09-04 · COMPLETE · SLURM 19899749*

Median η² is **0.0586 at every clone-size threshold**; effective programs per clone moves only
1.222 → 1.185 across a fourfold margin change; High-confidence restriction reproduces η² at Spearman
**0.901**, TF at **0.940**, PROGENy at **0.962**. LOO deliberately not applied to inherently
within-patient analyses. **Next** M50.

## M50 — Phase 6 Freeze

*2026-09-04 · COMPLETE · SLURM 19899749, 19899752*

**Object** `results/phase6/phase6_final_object.rds`, 6,060,757,871 bytes, md5
`b8c01dd01755dda10b2be4e0f5ef7ef7`, sha256
`bde592d461d22253b646da1426db8545404b9d31b7f1b948e23deed03aa643f6`, 19,716 cells, **241 metadata
columns** (196 Phase 1–5 + 46 Phase 6). **8 preservation guards + 8 reload validations, all passed.**
The full 690-regulon TF and 50-set Hallmark matrices stay as separate artefacts rather than bloating
the object.

**13/13 figures** · **12 tables** · manifest 31 sections · `FIGURE_INDEX.tsv` 706 → **732**.
`R_env_PRE_PHASE6.yml` and `R_env_POST_PHASE6.yml` byte-identical.

**Failures, both root-caused and fixed at source.** *19899746* — `m43` grouped by a column that
exists only after a rename. *19899748* — `decoupleR::get_collectri()` errored inside OmnipathR;
**OmnipathR was not upgraded**, the REST endpoint was used instead. **Visual QC** fixed one figure
defect (shared colour scale across rows spanning three orders of magnitude).

**Phase 7 is NOT initiated and requires separate authorization.**

---

# PHASE 4 STATUS: COMPLETE — SCEVAN EVIDENCE VISUALLY CONSOLIDATED
# PHASE 5 STATUS: COMPLETE
# PHASE 6 STATUS: COMPLETE
