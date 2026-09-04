# Phase 4 · M28 — SCEVAN Feasibility Report

*SLURM 19896711 · COMPLETED · 00:05:07 · 8 CPUs · ReqMem 96 G · MaxRSS 26.90 GiB (28% of request)*

**Verdict: SCEVAN is feasible on this dataset, and every §67 stop condition was checked and
cleared.** All four patients can be analysed independently at modest cost.

---

## 1. The required declarations (§11)

| Item | Value |
| --- | --- |
| **input object** | `results/phase2/phase2_final_object.rds` |
| **checksum** | md5 `153d5f6acc70f9c05aa48cabc4f4ac2d` — **verified equal** to the Phase 2 manifest |
| **input assay** | `RNA` (the Phase 2 default assay is `SCT`, which is *not* used — see §3) |
| **input layer** | the four per-sample `counts` layers: `counts.MPNST_1.1` … `counts.MPNST_4.4` |
| **gene identifier type** | **gene symbols** (0 of 31,764 features match Ensembl gene IDs) |
| **number of cells** | **19,716** |
| **number of genes** | 31,764 RNA features; 28,340 (89.2%) map to SCEVAN's chr1–22 annotation |
| **samples** | MPNST_1, MPNST_2, MPNST_3, MPNST_4 (`sample_id` = patient = dataset) |
| **normal-reference strategy** | primary: SCEVAN automatic (`norm_cell = NULL`). Sensitivity: 7,448 high-confidence immune cells, `FIXED_NORMAL_CELLS = FALSE`. See §5. |
| **planned SCEVAN mode** | per-patient `pipelineCNA(SUBCLONES = TRUE, ClonalCN = TRUE, plotTree = TRUE, organism = "human", beta_vega = 0.5)`, `par_cores` = SLURM CPUs |
| **resources** | M29: 8 CPUs, 64 G, 12 h, one array task per patient |

---

## 2. Object structure

```text
cells        19,716          duplicated cell IDs: 0
assays       RNA (31,764 features), SCT (29,113 features)   default: SCT
reductions   pca · umap_preintegration · postint_harmony · postint_umap_harmony
metadata     149 columns
RNA layers   counts.MPNST_1.1 … counts.MPNST_4.4 (+ data / scale.data per sample)
```

The four RNA count layers were mapped to samples **by cell membership, not by parsing the
layer name** — each layer was confirmed to contain cells from exactly one sample, and the four
layers to cover all four samples exactly once.

---

## 3. Raw counts, confirmed rather than assumed (§12)

SCEVAN must not receive Harmony coordinates, PCA, SCT residuals, scaled data or
log-normalised expression. The RNA `counts` layers were audited directly:

| Sample | genes expressed | cells | class | integer? | max | sparsity | median counts/cell | median genes/cell |
| --- | ---: | ---: | --- | :--: | ---: | ---: | ---: | ---: |
| MPNST_1 | 30,121 | 7,615 | dgCMatrix | **yes** | 6,215 | 0.8888 | 7,503 | 3,139 |
| MPNST_2 | 19,922 | 2,284 | dgCMatrix | **yes** | 3,397 | 0.8641 | 7,944 | 2,490 |
| MPNST_3 | 20,150 | 2,940 | dgCMatrix | **yes** | 10,839 | 0.9033 | 4,265 | 1,594 |
| MPNST_4 | 21,867 | 6,877 | dgCMatrix | **yes** | 5,007 | 0.9043 | 5,008 | 2,114 |

Every value satisfies `x == floor(x)`, so these are raw integer counts. The assertion is
re-executed at the top of each M29 job, not just here.

**The SCT assay is deliberately not used.** It carries four separate SCTransform models
(Phase 1 ran SCTransform layer-wise), so SCT values are not on a common footing across
samples — the same reason Phase 3 chose RNA + LogNormalize as its expression basis. SCEVAN
wants raw counts in any case.

### Memory-conscious extraction (§13)

`as.matrix()` was never applied to the full dataset. Each sample's counts were extracted as a
sparse `dgCMatrix`, restricted to genes expressed in that sample, and written to
`results/phase4/scevan/by_sample/<sample>/counts_raw.rds`. The four M29 jobs read those files
and **never reload the 6.0 GB Phase 2 object**.

SCEVAN does densify internally — `annotateGenes()` ends in `cbind(edb, as.matrix(mtx))` — so
the dominant single allocation was estimated per sample:

| Sample | dense annotated matrix |
| --- | ---: |
| MPNST_1 | 1.71 GiB |
| MPNST_4 | 1.12 GiB |
| MPNST_3 | 0.44 GiB |
| MPNST_2 | 0.34 GiB |

Adding SCEVAN's smoothed, relative and segmented copies plus a `parallelDist` distance matrix
(7,615² × 8 B = 0.46 GiB for the largest sample), **64 G is several times the working set**.
That is why M29 requests 64 G and not the 450 G maximum, and why `--cpus-per-task=8` is matched
exactly by `par_cores = 8` — no oversubscription.

---

## 4. Gene-identifier compatibility

`SCEVAN::annotateGenes()` matches row names against its internal `EnsDB_Hsapiens_v86` table,
restricted to chromosomes 1–22, preferring `gene_name` and falling back to `gene_id`.

```text
SCEVAN EnsDB_Hsapiens_v86 rows (chr1-22)  55,691
overlap with gene_name                    28,340   (89.2% of 31,764 RNA features)
overlap with gene_id                           0
=> SCEVAN matches on gene_name (symbols)
```

Per-sample annotatable genes: MPNST_1 26,970 · MPNST_2 18,153 · MPNST_3 18,371 · MPNST_4 19,848.

> **Limitation to carry forward**: SCEVAN's annotation covers **chromosomes 1–22 only**.
> Chromosome X and Y events are not assessable in this analysis. SCEVAN additionally removes
> cell-cycle genes and all `HLA-*` genes before inference — worth noting because HLA-E and
> HLA-F are Phase 3 findings, so the CNV analysis is structurally blind to their loci.

---

## 5. Normal-reference strategy, and the anti-circularity constraint (§15, §16)

| Group | Populations | Cells | Role |
| --- | --- | ---: | --- |
| High-confidence immune | CD4-T, CD8-T, NK, T-cell-other, B-cell, Plasma-cell, Macrophage, Monocyte, Dendritic, Plasmacytoid-DC | **7,448** | sensitivity-run reference only |
| Endothelial | Endothelial | 960 | **excluded** from the reference set |
| **Disputed** | MPNST-Tumor, Candidate-Malignant-Unresolved, Fibroblast, Pericyte-VSMC, Uncertain, Low-quality-excluded | **11,308** | **never a fixed reference — this is the question** |

Immune reference cells per patient: MPNST_1 2,036 · MPNST_2 762 · MPNST_3 2,323 · MPNST_4 2,327.

* **Primary run — `norm_cell = NULL`.** SCEVAN's `getConfidentNormalCells()` operates with **no
  input from Phase 2 labels whatsoever**. This is the run that carries inferential weight.
* **Sensitivity run** — supplies the 7,448 immune cells, `FIXED_NORMAL_CELLS = FALSE`.
* Endothelium is kept out of the reference set on purpose: it is a Phase 3 receiver of
  interest (VEGFA→KDR/FLT1/NRP1, the JAG/NOTCH circuit), so its malignancy call must remain
  independent rather than assumed.

### `FIXED_NORMAL_CELLS = TRUE` is prohibited in this project

Reading the installed SCEVAN source:

```r
if (FIXED_NORMAL_CELLS) {
    cellType_pred[!cellType_pred %in% norm_cell_names] <- "malignant"
    ...
}
```

Every cell outside the reference set is **forced** to malignant. That is exactly the
`if (!immune) MPNST-Tumor` inference the Phase 2 amendment banned, arriving through a function
argument. With `FIXED_NORMAL_CELLS = FALSE`, SCEVAN instead calls `classifyCluster()` on the
CNA-based hierarchical clustering, so classification retains genuine inferential content.

---

## 6. Cell-ID integrity (§67)

```text
cells in object                     19,716
duplicated cell IDs                      0
cells across extracted matrices     19,716   duplicated: 0
setequal(extracted, object cells)     TRUE
```

No cell was lost or duplicated in extraction. Example IDs: `MPNST_1_AAACCCAAGGACAGTC-1`.
IDs already carry the sample prefix, so cross-patient collision is structurally impossible.

Snapshots written once so downstream milestones are cheap:
`results/phase4/malignancy/phase4_cell_metadata.tsv.gz` (19,716 × 155, includes both UMAP
embeddings) and `results/phase4/malignancy/phase2_harmony_embedding.rds`
(`postint_harmony`, 19,716 × 30 — a **frozen copy**, never written back).

---

## 7. Technical smoke test — **not a scientific result**

600 randomly chosen MPNST_2 cells, seed 42, full pipeline (`SUBCLONES = TRUE`,
`ClonalCN = TRUE`):

```text
status                SUCCESS            elapsed        3.86 min
raw                   18,887 x 600       past filtering 8,069 -> 7,252 genes
confident normals     30                 after preprocessing  577 cells
returned              600 rows           columns: class, confidentNormal, subclone
class table           tumor 363 · normal 214 · filtered 23
subclones             5
files written         74
```

**Scientific use: none.** 600 random cells cannot support a malignancy call, and this output is
never propagated. Its only purpose was to prove the API, the gene annotation, the
classification path, the clonal-CN path, the subclone path and the plotting path all execute
end to end on real data of this shape. The scientific classification is the full per-patient
M29 run.

Two defects were found and fixed getting here — both documented in
`reports/phase4/environment/PHASE4_INSTALL_LOG.md` §3:

1. **SLURM 19896576** — SCEVAN 1.0.3 ignores `output_dir` in `getScevanCNV`,
   `getScevanCNVfinal`, `plotAllClonalCN`, `plotAllSubclonalCN`, `plotConsensusCNA` and
   `analyzeSegm2` (all hardcode `path = "./output"`), and `plotCNclonal()` does not forward it.
   Classification succeeded, then plotting died with `cannot open the connection`. Fixed by
   running each SCEVAN call with the process working directory set to its own folder and
   letting `output_dir` take SCEVAN's default. **The package was not patched.**
2. **SLURM 19896623** — `subcloneAnalysisPipeline` → `plotTSNE` requires the **Python**
   `umap-learn` via reticulate, and `plotTSNE` runs *before* the line that writes subclone
   labels into `classDf`, so its failure loses the clone assignments. Fixed by installing
   umap-learn 0.5.12 into an **isolated** `p4_umap_env` and pointing `RETICULATE_PYTHON` at it;
   `R_env` was verified unchanged before and after.

Neither fix involved raising a resource request.

---

## 8. Environment

```text
R              4.4.3 (2025-02-28)
SCEVAN         1.0.3      (Date 2025-02-12, built for R 4.4.3)
yaGST          2017.8.25
Seurat         5.4.0      unchanged from Phase 2/3
SeuratObject   5.3.0      unchanged
Matrix         1.7.4      unchanged
harmony        1.2.4      unchanged
sctransform    0.4.3      unchanged
```

SCEVAN and yaGST were **already installed**, so no R dependency was resolved for the CNV work.
The only addition anywhere was Python umap-learn in a separate environment.

---

## 9. Stop conditions checked (§67)

| Condition | Result |
| --- | --- |
| raw counts unavailable | **cleared** — four integer `counts` layers found and verified |
| cell IDs cannot be mapped | **cleared** — 1:1, zero duplicates, round-trip `setequal` TRUE |
| Phase 2 object corrupted | **cleared** — md5 matches the Phase 2 manifest exactly |
| SCEVAN fundamentally incompatible | **cleared** — 89.2% gene-symbol overlap; full pipeline ran |
| major unexplained cell loss | **cleared** — 19,716 in, 19,716 out |
| storage/environment failure | **cleared** |

**Proceeding to M29 automatically.**
