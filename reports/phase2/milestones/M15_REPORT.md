# Milestone M15 — Literature-Grounded Cell-Type Annotation

**Phase 2 · MPNST single-cell integration** · *Generated 2026-09-02*
*Status: **COMPLETE** — proceeded automatically to M16.*

> Full scientific argument and citations: [`reports/phase2/ANNOTATION_REPORT.md`](../ANNOTATION_REPORT.md)
> Machine-readable evidence: `results/phase2/annotation/ANNOTATION_EVIDENCE.tsv`

## 1. Scripts and configuration

Created: `scripts/R/phase2/annotate_celltypes.R`, `scripts/shell/phase2/run_m15_annotation.sh`,
**`config/phase2/annotation_map_M15.tsv`** (the auditable annotation map),
`reports/phase2/ANNOTATION_REPORT.md`, `reports/phase2/milestones/M15_REPORT.md`.

**Design decision:** the annotation is *data*, not code. Every label, its positive and
negative markers, supporting pathways, conflicting evidence, confidence and citation live
in one reviewable TSV. The script applies that map and does not infer labels, so the
scientific judgement is visible and reviewable rather than buried in R.

## 2. SLURM

| JobID | State | ExitCode | Elapsed | CPUs | Mem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| **19886713** | **COMPLETED** | **0:0** | **00:08:35** | 8 | 128G | **8,963,852 K (8.55 GiB)** — 6.7% |

No failures or retries. **R warnings: 0.** stderr contains only cosmetic ggplot
"Scale for colour is already present" messages from overriding the FeaturePlot colour
scale — no data-level warning. All 7 round-trip validation checks returned TRUE.

## 3. Literature research performed

Searches covered MPNST single-cell and spatial transcriptomics, MPNST malignant-state
markers, Schwann/neural-crest lineage biology, the MPNST/NF1 immune microenvironment, and
perineurial versus fibroblast identity in peripheral nerve. Ten primary sources with DOIs
or PMIDs are cited in `ANNOTATION_REPORT.md` §Sources and recorded per cluster in the
`source` and `doi_or_url` columns of `ANNOTATION_EVIDENCE.tsv`. No single paper was copied;
evidence was triangulated.

Key findings applied: SCP-like malignant cells are **L1CAM**-high; NC-like are
**APOD**/**NGFR**-high; **SHH** activation defines the NC-like MPNST-G1 subgroup; MPNST
down-regulates SOX10/CNP/PMP22/NGFR relative to neurofibroma, so absent classical Schwann
markers do not exclude a nerve-sheath lineage; CD163+ TAMs dominate the NF1-PNST myeloid
compartment; perineurial cells are ITGB4+/SLC2A1(GLUT1)+ and distinct from PI16+
fibroblasts.

## 4. Result

**Level 1** — Fibroblast/Stromal 5,499 (27.9%) · **Malignant/tumour 4,651 (23.6%)** ·
Myeloid 4,154 (21.1%) · B/Plasma 1,942 (9.8%) · T/NK 1,352 (6.9%) · Endothelial 960 (4.9%) ·
Uncertain 720 (3.7%) · Other/technical 438 (2.2%). Total 19,716.

**Level 2** — 17 cell types (see ANNOTATION_REPORT §2).
**Confidence** — High 10,115 (51.3%), Moderate 7,650 (38.8%), Low 1,951 (9.9%).

New metadata: `postint_celltype_level1/2/3`, `postint_annotation_confidence`,
`postint_annotation_initial`, `postint_annotation_source_cluster`.
Object: `results/phase2/annotation/phase2_harmony_annotated.rds`.

## 5. Malignant annotation — conservative, CNV-free

**No CNV inference** (outside Phase 2 scope). Malignant calls integrate Schwann/neural-crest
markers, MPNST literature, absence of immune/stromal/endothelial identity, and provenance —
never provenance alone, and never "non-immune ⇒ tumour". Labels are graded:

- **MPNST-like malignant (SCP-like)** — 1,680 (C8, C9). C9 retains **MPZ** with GFRA3,
  CRYAB, NOV and ABCB5; C8 expresses **L1CAM**.
- **MPNST-like malignant (NC-like)** — 440 (C14). **SHH** with GAL3ST1 and KLK6.
- **Schwann-lineage tumour-like** — 1,011 (C7). Coherent neural programme, no classical
  Schwann marker; deliberately weaker wording.
- **Candidate malignant** — 1,231 (C13, C17, C18), **Low confidence**. Neural programmes in
  93–99% single-patient clusters.
- **Cycling tumour-like** — 289 (C20). Proliferation programme masks lineage.

**Largest open question, stated explicitly:** the four fibroblast clusters (5,064 cells)
could contain MPNST **Mes-NC-like** malignant cells, which express mesenchymal/ECM
programmes. **The malignant fraction is 23.6% conservatively, and could be up to ~49% if
those clusters are malignant.** C25's CDKN2A expression argues against malignancy there
(CDKN2A is characteristically deleted in MPNST); C3's PI16+ epineurial signature argues for
genuine fibroblast. Resolving this needs CNV inference — the top Phase 3 priority.

## 6. Deliberate non-assignments

C12 (500 cells) — carries both a hypoxia programme and the reported perineurial
GLUT1/ITGB4 signature; both readings defensible, left **Uncertain**.
C23 (220 cells) — ribosomal-pseudogene/lncRNA dominated, left **Uncertain**.
C15 (438 cells) — mitochondrial-dominated, labelled `Other → Low-quality /
mitochondrial-high`; technical, matching the Phase 1 MPNST_4 flag.

A notable correction: the coarse M12 programme score had called **C21** B/plasma (pDCs
share IGJ and MZB1); the M14 markers LILRA4/CLEC4C/SPIB/GZMB are decisive for
**plasmacytoid dendritic cells**, and the annotation follows the markers.

## 7. Outputs

Tables: `ANNOTATION_EVIDENCE.tsv` (26 rows × 21 columns including `positive_markers`,
`negative_markers`, `supporting_pathways`, `conflicting_evidence`, `confidence`, `source`,
`doi_or_url`, `notes`, plus data-derived top-15 markers and per-sample counts),
`cluster_to_annotation_map.tsv`, `m15_annotation_record.json`, `prov_m15_annotation.json`,
`figure_index_m15.tsv`.

Figures (`results/phase2/figures/M15/`, PDF + PNG): broad-compartment UMAP, detailed
cell-type UMAP, cell-state UMAP, annotation-confidence UMAP, labelled cluster UMAP,
canonical marker dot plot, canonical marker heatmap, eight per-compartment FeaturePlot
panels, cluster→annotation mapping, confidence-by-compartment barplot.

## 8. Next

Proceeded automatically to **M16 — Annotation Refinement and Composition**.

---

# AMENDMENT ADDENDUM — M15A: CCC-Oriented Annotation Layer

*Added 2026-09-03 · SLURM JobID 19893067 · script `scripts/R/phase2/build_ccc_annotation.R`*

The M15 detailed annotation above is **unchanged and remains the source of truth**. This
addendum records the approved amendment that added a derived CCC-oriented layer.

## Detailed annotation method (recap)

Per cluster, from the M14 markers plus ten cited primary sources, recorded in the auditable
map `config/phase2/annotation_map_M15.tsv` and `ANNOTATION_EVIDENCE.tsv`. Hierarchical
Level 1/2/3 plus a confidence grade. No CNV inference.

## Malignant-cluster evidence (recap)

C9 retains **MPZ** with GFRA3, CRYAB, NOV, ABCB5 · C8 expresses **L1CAM** (reported SCP-like
marker) · C14 expresses **SHH** with GAL3ST1 and KLK6 (SHH-activated NC-like MPNST-G1
subgroup) · C7 shows a coherent Schwann/neural tumour programme with no competing lineage ·
C20 is a G2/M proliferation cluster inside the malignant compartment · C13, C17, C18 show
only patient-private neural programmes (Low confidence).

## CCC mapping strategy

Keyed on the **detailed Level 2 label**, applied mechanically from
`config/phase2/ccc_annotation_map.tsv`. The prohibited inference
`if (!immune) annotation_ccc <- "MPNST-Tumor"` was never used, and M17 machine-verifies that
no fibroblast, endothelial, immune, uncertain or low-quality cell entered `MPNST-Tumor`.

### Clusters mapped to `MPNST-Tumor` — 5 clusters, 3,420 cells (17.35%)

| Cluster | Detailed state | Cells | Evidence |
| --- | --- | ---: | --- |
| C8 | MPNST-like malignant (SCP-like) | 967 | L1CAM+, TGFA, PTHLH, MAL, OLFM1 |
| C9 | MPNST-like malignant (SCP-like) | 713 | MPZ+, GFRA3, CRYAB, NOV, ABCB5 |
| C14 | MPNST-like malignant (NC-like) | 440 | SHH, GAL3ST1, KLK6 |
| C7 | Schwann-lineage tumour-like | 1,011 | NLGN1, TENM2, NDST4, EPHB1, COL25A1, ROBO2 |
| C20 | Cycling tumour-like | 289 | TOP2A, UBE2C, NEK2, KIF20A, DLGAP5, ASPM |

### Clusters explicitly NOT mapped to tumour — 21 clusters

- **C13, C17, C18 (1,231 cells) — `Candidate-Malignant-Unresolved`.** The detailed label is
  explicitly provisional, confidence is **Low**, and the only evidence is neural-programme
  expression in clusters that are 93–99% one patient. Patient-private expression alone is
  weak evidence for malignancy, and absorbing these cells would inflate `MPNST-Tumor` with
  cells whose malignant identity is unestablished. Retained as a separate, visible,
  excludable identity that can be promoted later if CNV evidence supports it.
- **C0, C3, C5, C25 (5,064 cells) — `Fibroblast`.** Canonical fibroblast markers. The MPNST
  Mes-NC-like hypothesis is an unresolved possibility requiring CNV evidence, not established
  malignancy; absorbing them would be exactly the prohibited inference.
- **C16 — `Pericyte-VSMC`** · **C11, C19 — `Endothelial`** · **C1, C6 — `Macrophage`** ·
  **C24 — `Monocyte`** · **C10 — `Dendritic`** · **C21 — `Plasmacytoid-DC`** ·
  **C2 — `Plasma-cell`** · **C22 — `B-cell`** · **C4 — sub-resolved into CD8-T / CD4-T / NK /
  T-cell-other**.
- **C12, C23 (720 cells) — `Uncertain`**, retained, never absorbed.
- **C15 (438 cells) — `Low-quality-excluded`**, the mitochondrial artefact.

## Uncertain clusters and confidence

12.1% of cells (2,389) sit in the three not-CCC-ready identities. Annotation confidence is
carried unchanged from the detailed annotation and travels with the CCC label; two
Low-confidence states (C7, C20) contribute 38% of `MPNST-Tumor` and are quantified in
`MPNST_TUMOR_COMPOSITION.tsv` so their influence can be tested by exclusion.

## Figures (`results/phase2/figures/CCC_annotation/`, 14 figures × PDF + PNG)

`01_harmony_clusters` · `02_detailed_annotation` · **`03_ccc_annotation`** ·
**`04_mpnst_tumor_highlight`** · `05_immune_celltypes` · `06_tme_compartments` ·
**`07_ccc_annotation_marker_dotplot`** · `08_ccc_annotation_marker_heatmap` ·
`09_ccc_population_sizes` · `10_ccc_population_proportions` ·
`11_ccc_composition_by_sample` · `12_ccc_composition_by_patient` ·
**`14_mpnst_tumor_internal_composition`** · `15_annotation_confidence`.
Figure 13 (by condition) was **not** generated — no condition field exists;
`13_ccc_composition_by_condition_NOT_APPLICABLE.txt` records why.

## Tables

`CCC_ANNOTATION_MAPPING.tsv` (one row per cluster, every decision and every exclusion with
rationale, evidence and citation) · `CCC_ANNOTATION_SUMMARY.tsv` ·
`CCC_POPULATION_SIZE_AUDIT.tsv` · `MPNST_TUMOR_COMPOSITION.tsv` · `CCC_MARKER_SUMMARY.tsv` ·
`ccc_counts_by_{sample,patient}.tsv` · `ccc_proportions_by_{sample,patient}.tsv` ·
`ANNOTATION_EVIDENCE.tsv` augmented with `annotation_ccc`, `collapsed_to_mpnst_tumor`,
`ccc_compartment` and `ccc_mapping_rationale` (nothing deleted).

## Marker validation

**No marker rediscovery.** `CCC_MARKER_SUMMARY.tsv` reuses the M14 `filtered_cluster_markers.tsv`
aggregated to CCC identities, and the dot plot and heatmap use mean **SCT** expression.
Harmony coordinates were never reinterpreted as expression.

## SLURM

| JobID | State | ExitCode | Elapsed | CPUs | Mem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| **19893067** | **COMPLETED** | **0:0** | **00:07:33** | 8 | 128G | **8,853,320 K (8.44 GiB)** — 6.6% |

No failures or retries. All 16 round-trip validation checks TRUE, including proof that every
detailed annotation column is byte-identical and that no prohibited population entered
`MPNST-Tumor`.

## Warnings and limitations

- **1 recorded warning**: CD4 transcript detection is sparse in droplet scRNA-seq, so the
  `CD4-T` gate relies partly on IL7R/CCR7 and its boundary with `T-cell-other` is soft.
- The malignant fraction under this conservative collapse is 17.35%; if the fibroblast
  clusters prove Mes-NC-like malignant it could reach ~49%. CNV inference is required.
- `NK` (106 cells) is the thinnest CCC-ready population; flagged, not merged.
