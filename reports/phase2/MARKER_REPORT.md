# Cluster Marker Discovery Report

**Phase 2 · Milestone M14 · MPNST single-cell integration**
*Generated: 2026-09-02 · SLURM JobID 19886699*

> **STATISTICAL SCOPE.** These are exploratory **cluster-characterisation** markers, used
> to support the M15 annotation. They are **not** condition-level differential expression.
> Cells are not independent biological replicates, so no between-sample,
> between-patient or between-condition claim may be made from these tables. Every marker
> file carries this statement.

---

## 1. Assay and layer — inspected, not assumed

Phase 2 §14 requires that Seurat v5 layer structure be inspected explicitly. It was
(`results/phase2/markers/assay_layer_inspection.tsv`):

| Assay | Class | Default | Features | Layers | SCT models |
| --- | --- | :---: | ---: | --- | ---: |
| `RNA` | `Assay5` | no | 31,764 | `counts.MPNST_1.1; data.MPNST_1.1; scale.data.1; … ×4` (split four ways) | — |
| `SCT` | `SCTAssay` | **yes** | 29,113 | `counts; data; scale.data` (joined) | **4** |

**Consequence, established at M10 and confirmed here: the `SCT` assay carries four
SCTransform models, so `PrepSCTFindMarkers()` is mandatory.** It was run before any test,
recorrecting counts to a common sequencing depth so that cells normalised under different
models are comparable. Without it, differential expression across the four models would be
invalid.

The `RNA` assay was **not** used: its counts and data layers are split four ways and would
require `JoinLayers()` plus renormalisation, discarding the Phase 1 SCT model. Testing on
Harmony coordinates was never considered — Harmony is an embedding-level method and
carries no expression values.

## 2. Parameters

| Parameter | Value |
| --- | --- |
| Assay / layer | `SCT` / `data` |
| Preparation | `PrepSCTFindMarkers()` (4 models) |
| Test | `wilcox` (accelerated by `presto` 1.0.0) |
| `min.pct` | 0.25 |
| `logfc.threshold` | 0.25 |
| `only.pos` | TRUE |
| `recorrect_umi` | FALSE (already done by `PrepSCTFindMarkers`) |
| `random.seed` | 42 |
| Identities | `postint_harmony_primary_cluster` (resolution **1.0**, 26 clusters) |
| Robustness | repeated on `postint_harmony_alternative_cluster` (resolution 0.7, 21 clusters) |

`min.pct`, `logfc.threshold`, `only.pos` and the test match Phase 1
`config/config.yaml: markers`, keeping M14 comparable to the Phase 1 per-sample marker
runs.

## 3. Results

- **35,437 marker rows** across the 26 primary clusters; **31,774** pass adjusted
  p < 0.05 and log2FC ≥ 0.25.
- **Every cluster returned markers** — the minimum is 212 (C15) and the maximum 2,901
  (C18). No cluster is marker-less.
- A parallel run on the alternative resolution (0.7) is stored separately for robustness.

### Per-cluster summary (top markers by avg_log2FC)

| Cluster | n | Sig. markers | Dominant sample | Top markers |
| --- | ---: | ---: | --- | --- |
| C0 | 2,192 | 1,484 | MPNST_4 56% | SFRP2, CXCL14, C7, SFRP4, GAS1, MGST1, FBLN1, FGF7 |
| C1 | 1,907 | 769 | MPNST_1 36% | RNASE1, CD163, FOLR2, F13A1, SLCO2B1, C1QA, SIGLEC1, STAB1 |
| C2 | 1,709 | 338 | MPNST_4 43% | IGJ, MZB1, FKBP11, XBP1, SSR4, SEC11C, PRDX4 |
| C3 | 1,418 | 1,167 | MPNST_4 50% | ANGPTL7, SHISA3, APOD, PTGDS, ADH1B, CLDN1, PI16 |
| C4 | 1,352 | 494 | MPNST_3 49% | CD3D, KLRB1, CD3E, GZMA, CD3G, ICOS, CD2, GZMK, LCK, CCL5 |
| C5 | 1,313 | 1,271 | MPNST_4 78% | COMP, C1QTNF3, CILP2, DPT, TNMD, THBS4, CILP, COL1A1, FMOD |
| C6 | 1,158 | 1,275 | MPNST_1 39% | IL1A, CCL3L3, IL1B, CXCL3, CCL4, CCL3, BCL2A1, CXCL2 |
| C7 | 1,011 | 1,651 | MPNST_1 50% | NDST4, TENM2, NLGN1, KCNIP1, SULT1E1, RGS6, EPHB1, COL25A1 |
| C8 | 967 | 2,055 | MPNST_1 91% | ESPN, TGFA, MT1L, PTHLH, AZGP1, MAL, OLFM1, **L1CAM** |
| C9 | 713 | 2,263 | MPNST_1 80% | GFRA3, XKR4, SEMA3B, ABCB5, NOV, CRYAB, **MPZ**, RSPO3 |
| C10 | 612 | 1,129 | MPNST_4 33% | CD1C, CD1E, FCER1A, LGALS2, CLEC10A, FLT3 |
| C11 | 609 | 1,422 | MPNST_1 45% | ESM1, ANGPT2, APLN, DLL4, FLT1 |
| C12 | 500 | 600 | MPNST_1 57% | SLC2A1, PTHLH, IGFBP6, ADIRF, TGFA, VEGFA, NDRG1, ITGB4 |
| C13 | 485 | 2,500 | MPNST_1 99% | SULT1E1, COL26A1, ROBO2, CHRNA1, EDN3, SEMA3A, GAD1 |
| C14 | 440 | 2,692 | MPNST_1 94% | CLDN10-AS1, GAL3ST1, GRIA2, ANO3, ACPP, KLK6, **SHH**, LPL |
| C15 | 438 | 212 | MPNST_4 96% | MT-CYB, MT-CO3, MT-ATP6, MT-ND5, MT-CO2, MT-ND4 |
| C16 | 435 | 1,222 | MPNST_1 59% | FHL5, MYOCD, ACTG2, HIGD1B, COX4I2, PLN, MYH11, RGS5 |
| C17 | 384 | 763 | MPNST_1 93% | KIRREL3, TMEM178B, PRDM16, HS3ST4, NRXN1, NCAM2, SORCS3 |
| C18 | 362 | 2,901 | MPNST_1 99% | KCNJ6, GPC5, PPARGC1A, UNC5D, CHRM3, SLC8A3, KCNH1 |
| C19 | 351 | 964 | MPNST_4 33% | DARC(ACKR1), SELP, CCL14, MMRN1, SLCO2A1, RAMP3 |
| C20 | 289 | 1,198 | MPNST_4 74% | HIST1H3G, NEK2, KIF20A, DLGAP5, HJURP, UBE2C, TOP2A, ASPM |
| C21 | 287 | 997 | MPNST_3 57% | LILRA4, LRRC26, CLEC4C, GZMB, SCT, CLIC3, SPIB, PLAC8 |
| C22 | 233 | 230 | MPNST_3 40% | POU2AF1, CD79A, MZB1, FCRL5, XBP1, ZBP1 |
| C23 | 220 | 242 | MPNST_1 99% | RPL13AP5, SCG2, SNHG16, RPS7P10, RPS7P1 (ribosomal pseudogenes) |
| C24 | 190 | 986 | MPNST_3 38% | S100A12, S100A8, APOBEC3A, FCN1, S100A9, CD300E, AQP9 |
| C25 | 141 | 949 | MPNST_1 74% | CST2, ADH1C, COL10A1, CDKN2A, CST1, ITGBL1, MMP11 |

### Observations carried into M15

- **C9 retains MPZ**, the canonical Schwann myelin protein, alongside GFRA3 and the
  neural-crest stem marker ABCB5 — the strongest Schwann-lineage evidence in the dataset.
- **C8 expresses L1CAM**, reported as a marker of SCP-like malignant cells in peripheral
  nerve sheath tumours.
- **C14 expresses SHH** together with the Schwann-associated GAL3ST1 and KLK6.
- **C21's markers (LILRA4, CLEC4C, SPIB, GZMB) are decisive for plasmacytoid dendritic
  cells**, overriding the coarse M12 canonical-programme score, which had labelled it
  B/plasma because pDCs share IGJ and MZB1. This is an example of marker evidence
  correcting a programme-score prior.
- **C15's top markers are almost entirely mitochondrial** — a technical signature, not a
  cell type. This is the MPNST_4 stress population Phase 1 flagged (cluster C09,
  `percent.mt` R² = 0.47).
- **C23 is dominated by ribosomal pseudogenes and lncRNAs**, a low-complexity signature
  rather than a lineage.
- **C12 carries both a hypoxia programme (SLC2A1, VEGFA, NDRG1, ADM) and the reported
  perineurial signature (SLC2A1/GLUT1, ITGB4, ITGA6, CAV1)** — genuinely ambiguous, and
  left `Uncertain` in M15.

## 4. Tables (`results/phase2/markers/`)

`all_cluster_markers.tsv` (35,437 rows) · `filtered_cluster_markers.tsv` (31,774) ·
`top_markers_per_cluster.tsv` (= top 10) · `top10/top20/top50_markers_per_cluster.tsv` ·
`marker_summary_by_cluster.tsv` · `all_cluster_markers_alternative_resolution.tsv` ·
`marker_parameters.tsv` · `assay_layer_inspection.tsv` · `m14_marker_record.json` ·
`prov_m14_markers.json` · `figure_index_m14.tsv`

Columns: `cluster`, `gene`, `avg_log2FC`, `pct.1`, `pct.2`, `p_val`, `p_val_adj`,
`pct_diff`, `resolution`.

## 5. Figures (`results/phase2/figures/M14/`, PDF + PNG)

`M14_01_primary_cluster_umap` · `M14_02_marker_heatmap_top5` ·
`M14_03_marker_dotplot_top3` · `M14_04_featureplots_top_marker_per_cluster` ·
`M14_05_violin_top_markers` · `M14_06_markers_per_cluster`

## 6. Limitations

1. Markers are computed on the integrated clustering. Cluster boundaries inherit the M12
   caveats: B/plasma (C1) and fibroblast (C2) cohesion each fell 44% under Harmony, so
   markers for C2/C21/C22 and C0/C3/C5/C7/C15/C25 should be cross-checked against the
   non-integrated baseline before being relied on.
2. `only.pos = TRUE` means negative markers are not tabulated. Negative evidence quoted in
   the M15 annotation comes from inspecting canonical marker expression, not from these
   tables.
3. Wilcoxon p-values on 19,716 cells are not a meaningful measure of biological
   significance; effect size (`avg_log2FC`) and detection rates (`pct.1`, `pct.2`) carry
   the interpretation.
4. Clusters dominated by a technical programme (C15, mitochondrial) or a proliferation
   programme (C20) have lineage identity masked by that programme.
