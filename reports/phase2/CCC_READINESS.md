# Cell–Cell Communication Readiness Assessment

**Phase 2 · Milestone M17 · MPNST single-cell integration**
*Generated: 2026-09-03 · CCC layer built in M15A, SLURM JobID 19893067*

**Verdict: `READY WITH CAVEATS`** — see §7.

> **No cell–cell communication method was executed.** No CellChat, CellPhoneDB, LIANA,
> NicheNet or LochNESS was run. This is a readiness audit only.

---

## 1. Annotation

| Requirement | Status |
| --- | --- |
| `annotation_ccc` exists in the final object | ✅ 17 identities, no missing labels |
| Malignant population defined | ✅ `MPNST-Tumor`, 3,420 cells (17.35%) |
| Malignant definition is evidence-based, not "non-immune" | ✅ keyed on the detailed Level 2 label; machine-verified that no fibroblast, endothelial, immune, uncertain or low-quality cell entered `MPNST-Tumor` |
| Immune populations defined and resolved | ✅ 9 identities: CD8-T, CD4-T, NK, T-cell-other, Macrophage, Monocyte, Dendritic, Plasmacytoid-DC, B-cell, Plasma-cell |
| Stromal / endothelial retained separately | ✅ `Fibroblast`, `Pericyte-VSMC`, `Endothelial` |
| Uncertain cells documented, not absorbed | ✅ `Uncertain` (720), `Candidate-Malignant-Unresolved` (1,231), `Low-quality-excluded` (438) |
| Detailed annotation preserved as source of truth | ✅ every detailed column byte-identical after the CCC layer was added |
| Tumour heterogeneity recoverable | ✅ four detailed malignant states retained inside `MPNST-Tumor` |

Full per-cluster rationale, including every *exclusion*:
`results/phase2/annotation/CCC_ANNOTATION_MAPPING.tsv`.

---

## 2. Population size and replication

Thresholds applied: ≥100 cells total; present in ≥2 samples with ≥10 cells each.
Source: `results/phase2/annotation/CCC_POPULATION_SIZE_AUDIT.tsv`.

| `annotation_ccc` | Cells | Samples | Min/sample | Median/sample | Readiness |
| --- | ---: | ---: | ---: | ---: | --- |
| Fibroblast | 5,064 | 4 | 115 | 995 | **READY** |
| **MPNST-Tumor** | **3,420** | **4** | **97** | **448** | **READY** |
| Macrophage | 3,065 | 4 | 361 | 786 | **READY** |
| Plasma-cell | 1,709 | 4 | 72 | 450 | **READY** |
| Endothelial | 960 | 4 | 112 | 244 | **READY** |
| Dendritic | 612 | 4 | 100 | 156 | **READY** |
| CD4-T | 606 | 4 | 16 | 113 | **READY** |
| Pericyte-VSMC | 435 | 4 | 29 | 76 | **READY** |
| CD8-T | 382 | 4 | 26 | 105 | **READY** |
| Plasmacytoid-DC | 287 | 4 | 30 | 46 | **READY** |
| T-cell-other | 258 | 4 | 18 | 63 | **READY** |
| B-cell | 233 | 4 | 15 | 63 | **READY** |
| Monocyte | 190 | 4 | 27 | 46 | **READY** |
| NK | 106 | 4 | 11 | 26 | **READY** (thinnest population) |
| Candidate-Malignant-Unresolved | 1,231 | 4 | 6 | 15 | **EXCLUDE** — malignant identity unestablished |
| Uncertain | 720 | 4 | 21 | 98 | **EXCLUDE** — identity unresolved |
| Low-quality-excluded | 438 | 3 | 6 | 11 | **EXCLUDE** — technical artefact |

**Every CCC-ready population is present in all four samples.** That is the key result for
sample-aware communication analysis: each of the 14 ready identities can be evaluated
per sample and then compared across samples, rather than pooled.

**Replication caveat.** "Present in all four samples" is not the same as statistical
replication. Four tumours, one sample per patient, no condition variable — so
cross-sample comparison is descriptive. NK (106 cells, 11–43 per sample) is the thinnest
population and communication estimates involving it will be noisy in the smaller samples.
It was **not** merged into another type: the amendment forbids merging biologically distinct
immune types merely because they are small, so it is flagged instead.

---

## 3. `MPNST-Tumor` composition

`results/phase2/annotation/MPNST_TUMOR_COMPOSITION.tsv`

| Detailed malignant state | Cells | % of MPNST-Tumor | Cluster(s) | Confidence |
| --- | ---: | ---: | --- | --- |
| MPNST-like malignant (SCP-like) | 1,680 | 49.1% | C8, C9 | Moderate |
| Schwann-lineage tumour-like | 1,011 | 29.6% | C7 | Low |
| MPNST-like malignant (NC-like) | 440 | 12.9% | C14 | Moderate |
| Cycling tumour-like | 289 | 8.5% | C20 | Low |

The collapsed population is **not a black box**: about half of it is the best-supported
SCP-like compartment (L1CAM, MPZ), and the two Low-confidence contributors (C7, C20) are
identified and quantified so their influence on any communication result can be tested by
exclusion.

---

## 4. Expression data

| Requirement | Status |
| --- | --- |
| Appropriate expression layer retained | ✅ `SCT` (`counts`, `data`, `scale.data`; 29,113 features) and `RNA` (31,764 features, split per sample) |
| Expression **not** replaced by Harmony | ✅ Harmony wrote only the `postint_harmony` reduction; no expression value was modified anywhere in Phase 2 |
| Ligand/receptor genes available | ✅ spot-checked across the canonical panels used for annotation (e.g. `IL1B`, `CCL3`, `CCL4`, `CXCL2`, `CXCL3`, `CXCL14`, `ANGPT2`, `DLL4`, `FLT1`, `PDGFRA`, `PDGFRB`, `NOTCH3`, `CSF1R`, `TGFA`, `SHH`, `VEGFA`, `IL7R`, `CD4`, `CD8A`) |
| Assay for downstream testing documented | ✅ `SCT` with `PrepSCTFindMarkers()` applied, per M14 |

**Caveat for Phase 3 tooling.** The `RNA` assay's `counts` and `data` layers are **split
four ways** (`counts.MPNST_1.1` … `scale.data.4`). Any Phase 3 method requiring a single
joined RNA matrix must call `JoinLayers()` first. The `SCT` assay is already joined but
carries **four SCTransform models**, so cross-sample expression comparisons on `SCT`
require `PrepSCTFindMarkers()` as in M14.

---

## 5. Metadata

| Field | Present | Note |
| --- | :---: | --- |
| `sample_id` | ✅ | 4 levels; **simultaneously the dataset and the patient** |
| `orig.ident` | ✅ | cell-for-cell identical to `sample_id` |
| dataset | ✅ | = `sample_id` |
| patient | ✅ | = `sample_id`, one sample per patient |
| condition | ❌ | **does not exist** — no treatment, stage, site, NF1 status or grade anywhere in the data |
| QC covariates | ✅ | `nCount_RNA`, `nFeature_RNA`, `percent.mt`, `percent.ribo` |
| Phase 1 provenance | ✅ | all 93 original columns including every `preint_*` |

---

## 6. Statistical rule that Phase 3 must follow

> **Phase 3 must not pool every cell across every patient and treat cells as biological
> replicates.**

`sample_id` is preserved verbatim precisely so communication can be computed **per sample**
and only then compared:

```text
Per sample:  MPNST-Tumor -> Macrophage
Per sample:  MPNST-Tumor -> CD8-T
Per sample:  Macrophage  -> MPNST-Tumor
Per sample:  MPNST-Tumor -> Fibroblast / Endothelial
        then: patient/sample-aware comparison across the four tumours
```

With four tumours, one sample per patient and no condition variable, that comparison is
**descriptive**. There is no replication structure for condition-level inference, so
differential abundance and pseudobulk differential expression remain unsupportable on this
dataset.

A second constraint follows from M12: the integrated embedding was produced by correcting
`sample_id`, which is inseparable from patient identity. Communication *within* a sample is
unaffected by that correction, which is another reason to compute per sample rather than on
the pooled integrated space.

---

## 7. Verdict

# `READY WITH CAVEATS`

**Ready because:** `annotation_ccc` exists and is validated; `MPNST-Tumor` is defined from
positive malignant evidence rather than by exclusion; 14 populations pass the size audit and
**all 14 are present in all four samples**, enabling sample-aware analysis; the malignant,
immune, stromal and endothelial compartments are all separately addressable, so
tumour↔immune, tumour↔stromal and tumour↔endothelial axes can each be evaluated; expression
data is intact and was never replaced by Harmony; and the detailed annotation and tumour
heterogeneity are fully recoverable.

**Caveats that must be carried into Phase 3:**

1. **Exclude the three flagged populations** — `Candidate-Malignant-Unresolved` (1,231),
   `Uncertain` (720) and `Low-quality-excluded` (438), together 12.1% of cells.
2. **The malignant fraction is uncertain.** `MPNST-Tumor` is 17.35% under this conservative
   collapse. If Phase 3 CNV analysis shows the four fibroblast clusters are Mes-NC-like
   malignant, the true tumour compartment could reach ~49%, which would materially change
   any tumour↔stromal result. **CNV inference should precede substantive CCC conclusions.**
3. **Two Low-confidence states contribute 38% of `MPNST-Tumor`** (C7 Schwann-lineage
   tumour-like, C20 cycling tumour-like). Repeat key results with them excluded as a
   sensitivity check.
4. **NK is thin** (106 cells, 11–43 per sample). Interpret NK-involving results cautiously;
   do not merge it away.
5. **`CD4-T` versus `T-cell-other` is a soft boundary** — CD4 transcript detection is sparse,
   so the gate relies partly on IL7R/CCR7.
6. **No condition variable exists.** Cross-sample differences are descriptive only.
7. **Layer handling**: `JoinLayers()` for RNA-based tools; `PrepSCTFindMarkers()` for
   cross-sample SCT comparisons.
8. **Dataset = sample = patient = batch**, so no between-tumour claim is supportable from
   the integrated embedding.

---

## 8. Recommended Phase 3 architecture (not executed)

1. **CNV inference first** (inferCNV / CopyKAT) — resolves caveat 2 and would let
   `Candidate-Malignant-Unresolved` be promoted or discarded on evidence.
2. **Sample-aware LIANA consensus** — run per `sample_id`, aggregate consensus ranks across
   the four tumours, and report a communication event only where it reproduces in multiple
   samples. Prioritise `MPNST-Tumor ↔ Macrophage` (the largest, best-replicated pair, 361+
   macrophages in every sample) and `MPNST-Tumor ↔ CD8-T`.
3. **NicheNet receiver-response modelling** — use the M14 marker tables as the receiver gene
   sets, with `MPNST-Tumor` as sender and the immune identities as receivers, then reverse.
4. **Tumour-state-specific CCC** as a second pass, substituting the four detailed malignant
   states for `MPNST-Tumor` — possible without re-annotating, since those labels are retained.
5. **Sensitivity analyses**: exclude C7 and C20; exclude the flagged populations; repeat on
   the preserved non-integrated Phase 1 baseline for the M12-caveated B/plasma and fibroblast
   compartments.

**Not recommended on this dataset:** differential abundance, pseudobulk differential
expression, and any condition-level inference — no replication structure and no condition
variable exist.
