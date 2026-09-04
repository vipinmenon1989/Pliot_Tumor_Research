# Phase 2 Handoff — MPNST Harmony Integration, Clustering, Markers and Annotation

*Frozen 2026-09-03 · re-frozen 2026-09-03 after the approved CCC-oriented annotation amendment · Phase 2 status: **COMPLETE***

This document is self-contained: a new researcher or agent can reconstruct Phase 2 from it
without any conversation history.

---

## 1. Scientific objective

Take the four independently processed MPNST tumours frozen at the end of Phase 1, integrate
them at the embedding level with Harmony, and produce a single joint space in which cell
types can be clustered, characterised by markers and annotated across all four tumours —
while proving that the non-integrated Phase 1 baseline is preserved and that integration
did not destroy real biology.

---

## 2. Phase 1 input

| Item | Value |
| --- | --- |
| Object | `results/combined/pre_integration/combined_preintegration.rds` |
| md5 / sha256 | `88a442688f912d882f6c6da01820e329` / `c3fdce8b61602725977f989d8bbf10030ffe48ef6140372a13d2c59903a8dc66` |
| Cells / features | 19,716 · `SCT` 29,113, `RNA` 31,764 |
| Constituent tumours | MPNST_1 (7,615) · MPNST_2 (2,284) · MPNST_3 (2,940) · MPNST_4 (6,877) |
| Pre-Harmony reductions | `pca` (50 dims, assay SCT), `umap_preintegration` |
| Machine-readable manifest | `results/phase1_manifest.json` |

**Verified unchanged at the Phase 2 freeze** (md5 re-checked in M17; mtime still
2026-07-19 23:25). `processed_mpnst.rds` also untouched (mtime 2026-07-06 11:29).

Two Phase 1 structural findings, discovered in M10 and carried throughout:

- The `SCT` assay carries **four SCTransform models** — Phase 1 M8 ran SCTransform
  layer-wise, not once globally as its prose states. This does not invalidate integration,
  but it made `PrepSCTFindMarkers()` mandatory in M14.
- `percent.mt` is identically zero in all 7,615 MPNST_1 cells while M8 regressed it
  unconditionally. M11 verified this produced **zero** non-finite values anywhere. Closed.

---

## 3. Harmony decision (M10)

**Grouping variable: `sample_id` (4 levels).** It is the only field carrying batch
structure; `orig.ident` is cell-for-cell identical; every other candidate is prohibited
(legacy integration output), degenerate (zero-variance QC flags), continuous, or
non-existent.

**The central constraint of this dataset:** dataset = sample = patient = presumed technical
batch is a **single 4-level variable**. There is no clinical or technical covariate of any
kind. Correcting batch is mathematically indistinguishable from erasing between-patient
tumour biology. Full analysis: `reports/phase2/HARMONY_VARIABLE_DECISION.md`.

---

## 4. Harmony implementation (M11)

```r
set.seed(42)
harmony::RunHarmony(object = obj, group.by.vars = "sample_id", reduction.use = "pca",
                    dims.use = 1:30, reduction.save = "postint_harmony", verbose = TRUE)
```

harmony **1.2.4**. **Every scientific parameter at package defaults** — `theta` 2,
`sigma` 0.1, `lambda` 1, `nclust` 100 (auto cap), `max_iter` 10, `early_stop` TRUE,
`ncores` 1, `project.dim` TRUE, `harmony_options()` untouched. Only technical arguments were
supplied. Converged after **9 of 10** iterations (objective 753.09 → 353.66). A matrix-level
re-run with the same seed reproduced the embedding **bit-identically** (max abs diff 0.000),
so the result is deterministic and the convergence history describes the actual run.

Phase 1 `pca`, `umap_preintegration`, both graphs, both assays and all 93 original metadata
columns proven byte-identical after Harmony and after a disk round-trip.

---

## 5. Integration assessment (M12) — `ACCEPT DEFAULT HARMONY WITH CAVEATS`

All 19,716 cells, no subsampling. PRE = `pca` 1:30; POST = `postint_harmony` 1:30; matched
UMAP parameters and seed. `lisi`/`kBET` are not installed — inverse Simpson was computed
natively on exact kNN and a dominance ratio replaced kBET.

Because sample sizes are unequal, "fully mixed" means matching global composition:
same-sample fraction **0.3065**, entropy **1.2683 nats**, inverse Simpson **3.2627**.

| Metric (k = 15) | PRE | POST | Fully-mixed ref | Gap closed |
| --- | ---: | ---: | ---: | ---: |
| Same-sample NN fraction | 0.9728 | **0.7784** | 0.3065 | 29.2% |
| Neighbourhood entropy (nats) | 0.0490 | **0.3565** | 1.2683 | 25.2% |
| Inverse Simpson | 1.0549 | **1.4408** | 3.2627 | 17.5% |

**The decisive result is that mixing is compartment-specific.** Same-sample fraction
PRE → POST: Panleukocyte 0.933 → **0.527**, Myeloid 0.958 → **0.611**, T/NK 0.953 →
**0.689**, Endothelial 0.962 → **0.709** — versus Schwann/neural-crest 0.990 → **0.900**
and B/plasma 0.986 → **0.918**. Harmony aligned the shared immune and vascular compartments
while leaving the presumptive malignant lineage patient-private, which is the desired
behaviour and explains the modest global figure.

**Biology preserved:** within-sample kNN retention **0.8206**; coherence with the 54
integration-free Phase 1 clusters 0.8644 → 0.8427 (−2.5% relative); canonical-programme
silhouette 0.0627 → **0.0729** (improved), six of ten programmes more cohesive.

**Two caveats carried through every later milestone:**
- **C1** — B/plasma cohesion fell 44% (silhouette 0.437 → 0.243).
- **C2** — fibroblast cohesion fell 44% (0.174 → 0.097); MPNST_4 has the weakest
  within-sample preservation (0.781).

A methodological note worth keeping: the global `sample_id` **silhouette was −0.0104
pre-Harmony** despite 97% same-sample neighbours, because each sample spans malignant,
immune and stromal states. Batch structure here is *local*; only kNN metrics detect it.

Full report: `reports/phase2/HARMONY_ASSESSMENT.md`.

---

## 6. Clustering (M13)

`FindNeighbors(reduction = "postint_harmony", dims = 1:30, k.param = 20)` → graphs
`postint_harmony_nn` / `postint_harmony_snn`;
`RunUMAP(reduction = "postint_harmony", dims = 1:30, seed.use = 42)` →
`postint_umap_harmony`; `FindClusters(algorithm = 1, random.seed = 42)` at resolutions
0.1–1.0, **all ten preserved** as `postint_harmony_clusters_res_*`.

13 → 26 clusters. **No tiny cluster at any resolution** (smallest 141 cells). **Every
adjacent-resolution ARI ≥ 0.89**, so higher resolution splits groups rather than
reorganising them and the choice is not critical downstream.

**Primary = 1.0 (26 clusters)**, chosen by a coded composite rule (0.30 stability + 0.20
compactness + 0.15 compartment coverage + 0.15 programme purity + 0.10 low fragmentation +
0.10 granularity). **UMAP appearance was not a criterion.** *Caveat:* 1.0 is the sweep
boundary, so its stability term is one-sided and inflated. **Alternative = 0.7 (21
clusters)** is interior and two-sided-validated.

11 of 26 clusters are >60% one sample — expected, given M12. Full report:
`reports/phase2/CLUSTERING_ASSESSMENT.md`.

---

## 7. Marker discovery (M14)

Assay/layer structure was **inspected, not assumed**. `SCT` carries four models →
**`PrepSCTFindMarkers()` was mandatory and was applied**. `RNA` was not used (split layers).
**Harmony coordinates were never used for testing.**

`FindAllMarkers`, `SCT`/`data`, `wilcox` (presto), `min.pct = 0.25`,
`logfc.threshold = 0.25`, `only.pos = TRUE`, seed 42 — matching Phase 1
`config/config.yaml`. **35,437 marker rows; 31,774 significant; every cluster returned
markers** (212–2,901). The alternative resolution was run in parallel for robustness.

> These are **cluster-characterisation** markers for annotation. They are **not**
> condition-level differential expression. Cells are not independent biological replicates.

Full report: `reports/phase2/MARKER_REPORT.md`.

---

## 8. Cell-type annotation (M15/M16)

Hierarchical and **evidence-driven**. The annotation is *data, not code*: every label with
its positive/negative markers, supporting pathways, conflicting evidence, confidence and
citation lives in **`config/phase2/annotation_map_M15.tsv`**; the script applies it and
infers nothing.

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

17 Level 2 cell types: Fibroblast (5,064) · Macrophage (3,065) · Plasma cell (1,709) ·
MPNST-like malignant SCP-like (1,680) · T cell (1,352) · Candidate malignant (1,231) ·
Schwann-lineage tumour-like (1,011) · Endothelial cell (960) · Uncertain (720) · cDC2 (612) ·
MPNST-like malignant NC-like (440) · Low-quality/mitochondrial-high (438) · Mural (435) ·
Cycling tumour-like (289) · pDC (287) · B cell/plasmablast (233) · Monocyte (190).

Confidence after the M16 review: **High 9,828 (49.8%) · Moderate 6,637 (33.7%) · Low 3,251
(16.5%)**.

**M16 review:** six canonical compartment panels were scored per cluster and compared with
the assigned compartment under four documented rules. 18 confirmed, 5 confidence-downgraded,
3 deliberately unscored, **0 labels changed** — the panel score adjusts confidence but never
overrides a marker- and literature-grounded call. Both versions are preserved (`*_initial`
vs `*_refined`).

Full reports: `reports/phase2/ANNOTATION_REPORT.md`, `reports/phase2/milestones/M16_REPORT.md`.


---

## 8A. CCC-Oriented Annotation (approved amendment, 2026-09-03)

### 8A.1 Why it was created

The detailed M15/M16 annotation resolves the malignant compartment into five states across
eight clusters. That resolution is scientifically correct, but for cell–cell communication
analysis it fragments the tumour compartment into populations too small and too
patient-specific to act as a communication source or receiver. A simplified layer was
therefore added in which the **evidence-supported** malignant states are collapsed into a
single `MPNST-Tumor` identity, while immune, stromal and endothelial populations keep their
biologically meaningful identities so that tumour↔immune, tumour↔stromal and
tumour↔endothelial axes can each be evaluated.

**The detailed annotation remains the source of truth. It was not replaced or destroyed** —
every detailed column was verified byte-identical after the CCC layer was added.

Built as milestone **M15A** from the validated M16 object. **Harmony (M11), the integration
assessment (M12), clustering (M13) and marker discovery (M14) were not rerun.**

### 8A.2 How `MPNST-Tumor` was defined

> **`MPNST-Tumor` is NOT "everything that is not immune."**

The collapse is keyed on the **detailed Level 2 label**, which already encodes the strength
of the marker and literature evidence. The rule lives in
`config/phase2/ccc_annotation_map.tsv` and is applied mechanically:

| Detailed Level 2 label | → `MPNST-Tumor`? | Basis |
| --- | :---: | --- |
| MPNST-like malignant (SCP-like) | **YES** | L1CAM (C8), MPZ + GFRA3 + ABCB5 (C9); panel-confirmed in M16 |
| MPNST-like malignant (NC-like) | **YES** | SHH + GAL3ST1 + KLK6 (C14); panel-confirmed |
| Schwann-lineage tumour-like | **YES** | Coherent Schwann/neural tumour programme, no competing lineage (C7) |
| Cycling tumour-like | **YES** | Proliferating population inside the malignant compartment (C20) |
| **Candidate malignant** | **NO** | Provisional label, **Low** confidence, evidence is patient-private neural expression only (C13, C17, C18) |
| Everything outside Level 1 `Malignant / tumour` | **NO** | Never eligible |

Machine-verified in M17: no fibroblast, endothelial, pericyte, immune, uncertain or
low-quality cell entered `MPNST-Tumor`.

### 8A.3 Detailed tumour states collapsed

`results/phase2/annotation/MPNST_TUMOR_COMPOSITION.tsv`

| Detailed state | Cells | % of MPNST-Tumor | Cluster(s) | Confidence |
| --- | ---: | ---: | --- | --- |
| MPNST-like malignant (SCP-like) | 1,680 | 49.1% | C8, C9 | Moderate |
| Schwann-lineage tumour-like | 1,011 | 29.6% | C7 | Low |
| MPNST-like malignant (NC-like) | 440 | 12.9% | C14 | Moderate |
| Cycling tumour-like | 289 | 8.5% | C20 | Low |
| **MPNST-Tumor total** | **3,420** | **100%** | C7, C8, C9, C14, C20 | — |

All four states are retained in the final object, so tumour-state-specific communication
analysis remains possible later **without re-annotating**.

### 8A.4 `annotation_ccc` — all 17 identities

| Identity | Cells | % | Compartment | CCC-ready |
| --- | ---: | ---: | --- | :---: |
| Fibroblast | 5,064 | 25.68% | Fibroblast-Stromal | ✅ |
| **MPNST-Tumor** | **3,420** | **17.35%** | MPNST-Tumor | ✅ |
| Macrophage | 3,065 | 15.55% | Immune | ✅ |
| Plasma-cell | 1,709 | 8.67% | Immune | ✅ |
| Candidate-Malignant-Unresolved | 1,231 | 6.24% | Uncertain | ❌ |
| Endothelial | 960 | 4.87% | Endothelial | ✅ |
| Uncertain | 720 | 3.65% | Uncertain | ❌ |
| Dendritic (cDC2) | 612 | 3.10% | Immune | ✅ |
| CD4-T | 606 | 3.07% | Immune | ✅ |
| Low-quality-excluded | 438 | 2.22% | Other | ❌ |
| Pericyte-VSMC | 435 | 2.21% | Fibroblast-Stromal | ✅ |
| CD8-T | 382 | 1.94% | Immune | ✅ |
| Plasmacytoid-DC | 287 | 1.46% | Immune | ✅ |
| T-cell-other | 258 | 1.31% | Immune | ✅ |
| B-cell | 233 | 1.18% | Immune | ✅ |
| Monocyte | 190 | 0.96% | Immune | ✅ |
| NK | 106 | 0.54% | Immune | ✅ |

The immune compartment was **not** collapsed into a single "Immune" label — nine distinct
identities are retained. The single T/NK cluster (C4) was sub-resolved into CD8-T (382),
CD4-T (606), NK (106) and T-cell-other (258) by canonical marker gating restricted to that
cluster: NK = CD3(D/E/G)==0 & (NKG7|GNLY|KLRD1|NCAM1|KLRF1)>0; CD8-T = CD3>0 & (CD8A|CD8B)>0;
CD4-T = CD3>0 & CD8==0 & (CD4|IL7R|CCR7)>0; otherwise T-cell-other. No unsupported subtype
was invented.

### 8A.5 Populations kept separate, and uncertain populations

Kept separate rather than absorbed: `Fibroblast`, `Pericyte-VSMC`, `Endothelial`,
and every immune identity. Fibroblasts were **not** collapsed into tumour despite the
Mes-NC-like hypothesis (§9) — that hypothesis requires CNV evidence and is not established
malignancy; absorbing them would be exactly the prohibited inference.

Retained as explicitly not CCC-ready (12.1% of cells): `Candidate-Malignant-Unresolved`
(1,231), `Uncertain` (720, the hypoxia-vs-perineurial C12 and the low-complexity C23), and
`Low-quality-excluded` (438, the mitochondrial artefact C15).

### 8A.6 Population sizes and sample coverage

**All 14 CCC-ready populations are present in all four samples.** `MPNST-Tumor` has 97–2,428
cells per sample (median 448); `Macrophage` 361–1,133. The thinnest ready population is
`NK` (106 cells, 11–43 per sample), flagged rather than merged, since the amendment forbids
merging biologically distinct immune types merely because they are small. Full audit:
`results/phase2/annotation/CCC_POPULATION_SIZE_AUDIT.tsv`.

### 8A.7 Limitations

1. **The malignant fraction is uncertain.** `MPNST-Tumor` is 17.35% under this conservative
   collapse. If the four fibroblast clusters prove Mes-NC-like malignant, the true compartment
   could reach ~49% — which would materially change any tumour↔stromal result. **CNV
   inference should precede substantive CCC conclusions.**
2. **Two Low-confidence states contribute 38% of `MPNST-Tumor`** (C7, C20). Repeat key
   results with them excluded as a sensitivity check.
3. `Candidate-Malignant-Unresolved` may be promoted later on CNV evidence, raising the
   tumour fraction.
4. `CD4-T` versus `T-cell-other` is a soft boundary — CD4 detection is sparse in droplet data.
5. `NK` is thin; interpret NK results cautiously.
6. No condition variable exists, so cross-sample comparison is descriptive.
7. Dataset = sample = patient = batch, so no between-tumour claim is supportable.

### 8A.8 How Phase 3 should use it

> **Phase 3 must not pool every cell across every patient and treat cells as biological
> replicates.**

`sample_id` is preserved verbatim so communication is computed **per sample** and only then
compared:

```text
Per sample:  MPNST-Tumor -> Macrophage
Per sample:  MPNST-Tumor -> CD8-T
Per sample:  Macrophage  -> MPNST-Tumor
Per sample:  MPNST-Tumor -> Fibroblast / Endothelial
        then: patient/sample-aware comparison across the four tumours
```

Exclude the three not-ready populations. Run CNV inference first. Then sample-aware LIANA
consensus (prioritising `MPNST-Tumor ↔ Macrophage`, the largest and best-replicated pair)
and NicheNet receiver-response modelling using the M14 marker tables as receiver gene sets.
A second pass can substitute the four detailed malignant states for `MPNST-Tumor` without
re-annotating. Full assessment: **`reports/phase2/CCC_READINESS.md`** — verdict
**`READY WITH CAVEATS`**.

---

## 9. MPNST annotation evidence

The rule "non-immune ⇒ tumour" was **not** used. Malignant calls integrate Schwann/
neural-crest lineage markers, MPNST-specific literature, absence of convincing
immune/stromal/endothelial identity, and provenance — never provenance alone. **No CNV
inference** (outside Phase 2 scope). Labels are graded:

| Label | Cells | Clusters | Key evidence |
| --- | ---: | --- | --- |
| MPNST-like malignant (SCP-like) | 1,680 | C8, C9 | **C9 retains MPZ** with GFRA3, CRYAB, NOV, ABCB5; **C8 expresses L1CAM**, the reported SCP-like marker |
| MPNST-like malignant (NC-like) | 440 | C14 | **SHH** with GAL3ST1 and KLK6 — the reported SHH-activated NC-like MPNST-G1 subgroup |
| Schwann-lineage tumour-like | 1,011 | C7 | Coherent neural adhesion programme, no classical Schwann marker — deliberately weaker wording |
| Candidate malignant | 1,231 | C13, C17, C18 | Neural programmes in 93–99% single-patient clusters; **Low confidence** |
| Cycling tumour-like | 289 | C20 | G2/M programme masks lineage |

### The largest open question

The four fibroblast clusters (C0, C3, C5, C25; 5,064 cells) carry canonical fibroblast
markers, but MPNST also contains **Mes-NC-like** malignant cells that express mesenchymal
and ECM programmes.

> **The malignant fraction is 23.6% under the conservative annotation, and could be as high
> as ~49% if those clusters are Mes-NC-like malignant.** C25's CDKN2A expression argues
> against malignancy there (CDKN2A is characteristically deleted in MPNST); C3's PI16+
> epineurial signature argues for genuine fibroblast. **Resolving this requires CNV
> inference — the top Phase 3 priority.**

Key literature (full list in `ANNOTATION_REPORT.md` §Sources):
[Neuro-Oncology doi:10.1093/neuonc/noaf016](https://academic.oup.com/neuro-oncology/advance-article-abstract/doi/10.1093/neuonc/noaf016/7976860) ·
[Sci Adv doi:10.1126/sciadv.abo5442](https://www.science.org/doi/10.1126/sciadv.abo5442) ·
[PMC12204358](https://pmc.ncbi.nlm.nih.gov/articles/PMC12204358/) ·
[PMID 39321200](https://pubmed.ncbi.nlm.nih.gov/39321200/) ·
[PNAS doi:10.1073/pnas.1913444117](https://www.pnas.org/doi/10.1073/pnas.1913444117) ·
[PMID 24719203](https://pubmed.ncbi.nlm.nih.gov/24719203/).

---

## 10. Composition

**Descriptive only** — no inferential condition-level or differential-abundance test was
performed; every table carries that note. `sample_id` is simultaneously dataset and patient,
so those groupings are one table. **No condition table exists** — no such field is present
in the data, and none was fabricated.

**4 of 17 cell types are driven >80% by a single patient** and must not be read as general
MPNST biology:

| Cell type | Cells | Dominant | % |
| --- | ---: | --- | ---: |
| Candidate malignant | 1,231 | MPNST_1 | 97.2% |
| Low-quality / mitochondrial-high | 438 | MPNST_4 | 96.1% |
| MPNST-like malignant (NC-like) | 440 | MPNST_1 | 94.3% |
| MPNST-like malignant (SCP-like) | 1,680 | MPNST_1 | 86.3% |

---

## 11. Final object

| Item | Value |
| --- | --- |
| Path | **`results/phase2/phase2_final_object.rds`** |
| Size | 6,055,929,223 bytes (5.64 GiB) |
| **md5** | **`153d5f6acc70f9c05aa48cabc4f4ac2d`** |
| **sha256** | **`62e97524836309170c38be9335a10a4c76baa48a246b2963d2d440cebead0cf3`** |
| Cells | 19,716 |
| Assays | `RNA` (31,764, split layers), `SCT` (29,113, 4 models) — default `SCT` |
| Reductions | `pca`, `umap_preintegration`, `postint_harmony` (30), `postint_umap_harmony` |
| Graphs | `SCT_nn`, `SCT_snn`, `postint_harmony_nn`, `postint_harmony_snn` |

Metadata: all 93 Phase 1 columns preserved verbatim (including legacy `seurat_clusters`,
`orig.anno` and every `preint_*`), plus
`postint_harmony_clusters_res_0.1`…`_res_1.0`, `postint_harmony_primary_cluster`,
`postint_harmony_alternative_cluster`, `postint_primary_resolution`,
`postint_alternative_resolution`, ten `m12_progscore_*`, `m12_program_argmax`,
`postint_celltype_level1/2/3`, `postint_annotation_confidence`,
`postint_celltype_level1/2/3_initial`, `postint_annotation_confidence_initial`,
`postint_celltype_level1/2/3_refined`, `postint_annotation_confidence_refined`,
`postint_annotation_initial`, `postint_annotation_refined`,
`postint_annotation_review_rule`, `postint_annotation_review_outcome`,
`postint_annotation_source_cluster`, and the CCC layer **`annotation_ccc`**,
`annotation_ccc_compartment`, `annotation_ccc_is_tumor`, `annotation_ccc_ccc_ready`.

Validated by **37 pre-save checks** and **18 post-reload checks** in M17, including
scientific-scope guards confirming no CNV, pseudobulk-DE or trajectory columns exist, and
eleven CCC-specific guards confirming that `annotation_ccc` is complete and that no
fibroblast, endothelial, immune, uncertain or low-quality cell entered `MPNST-Tumor`.

---

## 12. Final figures — `results/phase2/figures/final/`

`01_pre_vs_post_harmony_dataset_{PRE,POST}` · `02_pre_vs_post_harmony_sample` ·
`03_harmony_primary_clusters` · `04_harmony_broad_celltypes` ·
`05_harmony_detailed_celltypes` · `06_canonical_marker_dotplot` ·
`07_cluster_marker_heatmap` · `08_celltype_composition_by_sample` ·
`09_celltype_composition_by_patient` · `11_annotation_confidence` ·
`12_resolution_sweep_panel` · `13_harmony_knn_retention` ·
`14_sample_contribution_per_celltype` (all PDF + PNG) ·
`10_celltype_composition_by_condition_NOT_APPLICABLE.txt`.

Per-milestone suites: `results/phase2/figures/{M12,M13,M14,M15,M16}/` and
`reports/phase2/figures/m12/`. All 178 Phase 2 figure files are registered in
`reports/FIGURE_INDEX.tsv` (438 rows total, including 260 back-filled Phase 1 rows).

---

## 13. Final tables — `results/phase2/tables/final/` (17 files)

`cluster_sizes.tsv` · `cluster_composition.tsv` · `top_markers_per_cluster.tsv` ·
`all_cluster_markers.tsv` · `marker_summary_by_cluster.tsv` · `annotation_evidence.tsv` ·
`annotation_refinement_review.tsv` · `cell_counts_by_annotation.tsv` ·
`cell_counts_by_annotation_level1.tsv` · `cell_proportions_by_sample.tsv` ·
`cell_proportions_by_patient.tsv` · `celltype_sample_dependence.tsv` ·
`clustering_comparison_table.tsv` · `harmony_parameters.tsv` ·
`harmony_pre_post_mixing_summary.tsv` · `harmony_biological_preservation_summary.tsv` ·
`cell_proportions_by_condition_NOT_APPLICABLE.txt`.

---

## 14. Reproduction commands

```bash
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

sbatch scripts/shell/phase2/run_m10_inspect.sh       # validate the Phase 1 handoff
sbatch scripts/shell/phase2/run_m11_harmony.sh       # default Harmony
sbatch scripts/shell/phase2/run_m12_evaluate.sh      # pre/post assessment
sbatch scripts/shell/phase2/run_m12_supplement.sh    # named pre/post panels
sbatch scripts/shell/phase2/run_m13_clustering.sh    # neighbours, UMAP, resolution sweep
sbatch scripts/shell/phase2/run_m14_markers.sh       # cluster markers
sbatch scripts/shell/phase2/run_m15_annotation.sh    # annotation (uses config/phase2/annotation_map_M15.tsv)
sbatch scripts/shell/phase2/run_m16_composition.sh   # refinement + composition
sbatch scripts/shell/phase2/run_m17_freeze.sh        # validation, freeze, manifest
```

Strictly sequential — each consumes the previous milestone's object. Every script accepts
`--validation-mode` for synthetic smoke testing without real data. **No Snakemake is used in
Phase 2.** All seeds are 42.

---

## 15. Environment

R 4.4.3 (`/local/projects-t3/lilab/vmenon/anaconda3/envs/R_env/bin/R`), conda env `R_env`.
Seurat 5.4.0 · SeuratObject 5.3.0 · **harmony 1.2.4** · Matrix 1.7.4 · sctransform 0.4.3 ·
glmGamPoi 1.18.0 · presto 1.0.0 · cluster 2.1.8.1 · RANN 2.6.2 · uwot 0.2.4 · aricode 1.0.3 ·
ggplot2 4.0.1 · patchwork 1.3.2 · reshape2 1.4.5.

**Nothing was installed, upgraded or removed during Phase 2.** `lisi`, `kBET` and `clustree`
are **not** installed; native substitutes were used and are documented in
`results/phase2/harmony/evaluation/harmony_evaluation_parameters.tsv`.

Records: `reports/phase2/PHASE2_ENVIRONMENT.tsv`, `reports/phase2/PHASE2_SESSIONINFO.txt`,
`workflow/envs/R_env_portable.yaml`.

---

## 16. SLURM resources

All jobs: account `ihc`, partition `ihc`, node `ihc-grid-1-1-1`.

| Milestone | JobID | State | Elapsed | CPUs | Mem req | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| M10 inspect | 19886411 | COMPLETED | 00:02:17 | 4 | 96G | 4.73 GiB |
| M11 Harmony | 19886486 | COMPLETED | 00:08:24 | 8 | 64G | 8.38 GiB |
| M12 evaluate | 19886628 | COMPLETED | 00:04:09 | 8 | 96G | 8.14 GiB |
| M12 supplement | 19886682 | COMPLETED | 00:01:32 | 4 | 64G | 7.18 GiB |
| M13 clustering | 19886690 | COMPLETED | 00:09:37 | 8 | 128G | 12.09 GiB |
| M14 markers | 19886699 | COMPLETED | 00:04:26 | 12 | 250G | 16.40 GiB |
| M15 annotation | 19886713 | COMPLETED | 00:08:35 | 8 | 128G | 8.55 GiB |
| M16 composition | 19886728 | COMPLETED | 00:07:34 | 8 | 128G | 8.48 GiB |
| M17 freeze | 19886743 | COMPLETED | 00:07:39 | 8 | 128G | 8.13 GiB |

**Total successful compute: 54 minutes.** Peak memory anywhere: 16.40 GiB (M14,
`PrepSCTFindMarkers`) — 3.6% of the 450G envelope.

**Four documented failures**, all code defects fixed at source; **no resource was ever
increased in response to a failure**:
19886675 (`dpi = NA` rejected by ggplot2 4.0.1) · 19886683 (preservation guard fired; the
whole-frame metadata digest was replaced with a stricter per-column check) · 19886685
(unnamed-vector indexing) · 19886687 (`sprintf("%d", median())`).

---

## 17. Provenance

Per-milestone provenance JSONs (input paths + md5, outputs + md5, parameters, seeds, git
commit, git status, `sessionInfo()`, SLURM JobID) written through the Phase 1
`scripts/R/provenance_utils.R` helpers: `results/phase2/harmony/prov_m11_harmony.json` ·
`results/phase2/harmony/evaluation/prov_m12_evaluation.json` ·
`results/phase2/figures/M12/prov_m12_supplement.json` ·
`results/phase2/clustering/prov_m13_clustering.json` ·
`results/phase2/markers/prov_m14_markers.json` ·
`results/phase2/annotation/prov_m15_annotation.json` ·
`results/phase2/composition/prov_m16_composition.json` ·
`results/phase2/prov_m17_freeze.json`.

Consolidated machine-readable record: **`results/phase2/phase2_manifest.json`**.
Milestone reports: `reports/phase2/milestones/M10_REPORT.md` … `M17_REPORT.md`.

---

## 18. Known limitations

1. **Dataset = sample = patient = batch is one variable.** Residual structure cannot be
   apportioned between uncorrected batch effect and real between-patient tumour biology.
   **No between-tumour, between-patient, between-condition or differential-abundance claim
   can be supported from the integrated embedding.**
2. **No biological-condition, clinical or technical covariate exists.** Several requested
   condition figures and tables were deliberately not produced; each absence is recorded in
   a `*_NOT_APPLICABLE.txt` note rather than fabricated.
3. **No CNV inference.** All malignant calls are expression- and provenance-based.
4. **The malignant fraction is 23.6% conservatively and could be ~49%** if the fibroblast
   clusters are Mes-NC-like malignant. Largest quantitative uncertainty in Phase 2.
5. **M12 caveats C1 and C2** apply to clusters C2/C21/C22 and C0/C3/C5/C7/C15/C25.
6. **Primary resolution 1.0 sits at the sweep boundary**; its stability term is one-sided.
   Alternative 0.7 is interior and two-sided-validated.
7. **Four SCT models persist.** `PrepSCTFindMarkers()` was applied for markers, but residual
   per-sample normalisation differences remain in the SCT residuals.
8. **C15 (438 cells) is a technical mitochondrial artefact** and should be excluded from
   biological interpretation.
9. **Annotation is per cluster**, so C19 (venous + lymphatic) and C4 (T and NK) are
   under-resolved. 16.5% of cells carry Low confidence; 3.7% are Uncertain.
10. **Phase 1 documentation defects D1–D9** (`M10_REPORT.md` §7) remain uncorrected pending
    researcher instruction. D8 — `config/config.yaml: input_rds` points at the deleted
    `/local/projects-t3/lilab/vmenon/Pilot_tumor/` root — would break a Phase 1 re-run from
    the original RDS but does not affect Phase 2.
11. Only Harmony was evaluated. CCA, RPCA and MNN were not benchmarked.

---

## 19. Recommended Phase 3 analyses

**Not executed. None is authorised by Phase 2.**

1. **CNV inference (inferCNV / CopyKAT) — highest priority.** The one analysis that would
   resolve the 23.6%-vs-49% malignant-fraction question by testing whether the fibroblast
   clusters are Mes-NC-like malignant, and would upgrade "candidate malignant" (C13, C17,
   C18) from Low confidence.
2. **Malignant-state refinement** — re-cluster the malignant compartment alone to resolve
   SCP-like / NC-like / Mes-NC-like states and test the MPNST-G1 (SHH) versus G2 (WNT)
   subgrouping suggested by C14.
3. **Tumour–immune interaction analysis** — the immune compartment (T/NK 1,352; myeloid
   4,154 including FOLR2+ and inflammatory macrophages, cDC2, pDC, monocytes) is well
   integrated and well annotated, so this is well supported.
4. **Cell–cell communication (CellChat / CellPhoneDB / LIANA / NicheNet)** — feasible, but
   any result crossing the malignant boundary inherits limitation 4.
5. **Pathway / gene-programme analysis** within compartments — descriptive use is supported.
6. **Differential abundance and pseudobulk DE — NOT recommended on this dataset.** With one
   sample per patient, four patients and no condition variable, there is no replication
   structure and no contrast to test. Both would require additional samples.
7. **Sub-clustering C19** (venous vs lymphatic endothelium) and **C4** (T vs NK).
8. **Exclude C15** (mitochondrial-high) and treat C12/C23 as unresolved in any downstream
   analysis.
9. Optionally, a **Harmony `theta` sensitivity analysis** ({1, 2, 4}) re-scored through the
   M12 metric suite — designed and documented in `HARMONY_ASSESSMENT.md` §J, deliberately
   not run.

### Suitability summary

| Phase 3 analysis | Suitable? |
| --- | --- |
| Tumour–immune interaction | **Yes** |
| Cell–cell communication | **Yes**, with the malignant-identity caveat |
| Pathway / state analysis | **Yes**, descriptively, within compartments |
| CNV inference | **Yes — and it is the top priority** |
| LochNESS | Possible, but interpret against limitation 1 |
| Differential abundance | **No** — no replication structure, no condition variable |
| Pseudobulk differential expression | **No** — same reason |
