# Harmony Grouping-Variable Decision

**Phase 2 · Milestone M10 · MPNST single-cell integration**

*Created: 2026-09-02*
*Status: **PROPOSED — requires explicit researcher approval before M11***

---

## 0. Why this document exists

Phase 2 §8 forbids integrating on `sample_id`, `orig.ident`, patient or any other field
merely because it exists. Phase 2 §8 further requires that, if the intended Harmony
variable is strongly or perfectly confounded with biology, execution **stops** and the
researcher reviews the decision.

**That condition is met here.** In this dataset the batch variable and the biological
donor are the same variable. This document sets out the evidence, the proposed decision,
and the risks the researcher is being asked to accept.

---

## 1. Evidence base

All statements below are derived from the frozen Phase 1 records and from a read-only
structural inspection of the Phase 1 handoff object executed inside SLURM job
`19886411` (script `scripts/R/phase2/inspect_phase1_handoff.R`; outputs in
`results/phase2/handoff/`). No Harmony was run.

| Source | Path |
| --- | --- |
| Phase 1 manifest | `results/phase1_manifest.json` |
| Phase 1 metadata inventory | `reports/combined/pre_integration/metadata_inventory.tsv` |
| Phase 1 confounding cross-tab | `reports/combined/pre_integration/confounding_summary.tsv` |
| Phase 1 pre-integration assessment | `reports/PRE_INTEGRATION_ASSESSMENT.md` |
| Phase 1 integration preparation | `reports/INTEGRATION_PREPARATION.md` |
| M10 object structure | `results/phase2/handoff/handoff_structure.json` |
| M10 batch cross-tab | `results/phase2/handoff/candidate_batch_crosstab.tsv` |
| M10 per-sample QC | `results/phase2/handoff/qc_summary_by_sample.tsv` |

---

## 2. Candidate variables considered

Every one of the 93 metadata columns in the handoff object was screened. Only the
following were viable candidates for a Harmony grouping variable.

| Candidate | Cardinality | Verdict | Reason |
| --- | --- | --- | --- |
| `sample_id` | 4 | **SELECTED** | Canonical dataset splitter established and validated in Phase 1 M1; the field Phase 1 used to construct every constituent object. |
| `orig.ident` | 4 | Rejected as redundant | Verified cell-for-cell **identical** to `sample_id` (M10 inspection: `sample_id_equals_orig_ident = true`; off-diagonal counts all zero in `candidate_batch_crosstab.tsv`). Integrating on it would be the same operation under a different name. |
| `orig.anno` | 16 | Rejected — prohibited and unsound | Legacy *biological* cell-type labels from the original analysis. Phase 2 §6 forbids using legacy analysis results to choose Harmony variables, and using a cell-type label as a batch variable would deliberately destroy the biology Harmony must preserve. |
| `unintegrated_clusters`, `seurat_clusters`, `pca100_*_clusters`, `pca100.sct_*_clusters` | 10–21 | Rejected — prohibited | Legacy clustering and legacy *integration* output (Harmony/CCA/RPCA/MNN). Phase 2 §6 explicitly bars reuse. |
| `doublet_class`, `retained`, `fail_*` | 1 | Rejected | Zero variance after Phase 1 filtering (all cells `singlet` / `TRUE` / `FALSE`). Nothing to correct. |
| `percent.mt`, `percent.ribo`, `nCount_RNA`, `nFeature_RNA` | continuous | Not applicable | Harmony's `vars_use` requires categorical grouping variables. These are continuous QC covariates and are handled at the normalization stage, not by Harmony. They will instead be used as **diagnostics** in M12. |
| `preint_*` | various | Rejected | Phase 1 provenance columns (per-dataset cluster sweeps, recommended resolutions). Immutable provenance, not batch structure. |
| Library prep batch, capture channel, sequencing run, flow cell, operator, date | — | **Do not exist** | Confirmed absent in Phase 1 M8 (`INTEGRATION_PREPARATION.md` Q3/Q4) and re-confirmed in the M10 metadata inventory. |
| Patient ID, sex, age, anatomical site, primary/metastatic, treatment status, NF1 status, tumour grade | — | **Do not exist** | No clinical covariate of any kind is recorded in the object. |

---

## 3. Cardinality and cross-tabulation

### 3.1 `sample_id` × `orig.ident` (from `candidate_batch_crosstab.tsv`)

| | MPNST_1 | MPNST_2 | MPNST_3 | MPNST_4 |
| --- | ---: | ---: | ---: | ---: |
| **MPNST_1** | **7,615** | 0 | 0 | 0 |
| **MPNST_2** | 0 | **2,284** | 0 | 0 |
| **MPNST_3** | 0 | 0 | **2,940** | 0 |
| **MPNST_4** | 0 | 0 | 0 | **6,877** |

Perfectly diagonal. The two fields are interchangeable.

### 3.2 Group sizes

| `sample_id` | Cells | Share |
| --- | ---: | ---: |
| MPNST_1 | 7,615 | 38.62% |
| MPNST_2 | 2,284 | 11.58% |
| MPNST_3 | 2,940 | 14.91% |
| MPNST_4 | 6,877 | 34.88% |
| **Total** | **19,716** | 100% |

Group sizes are unbalanced by a factor of 3.3× (MPNST_1 vs MPNST_2). Harmony's default
`theta = 2` is applied uniformly to all levels of the grouping variable and does not
compensate for imbalance; the smaller samples therefore carry more leverage per cell
in the diversity penalty. This is recorded as a diagnostic to watch in M12, not as a
reason to deviate from defaults in M11.

---

## 4. Relationship to patient, dataset and biological condition

| Relationship | Status |
| --- | --- |
| `sample_id` ↔ dataset | **Identical**. Each Phase 1 constituent dataset is exactly one `sample_id`. |
| `sample_id` ↔ patient | **Identical (1:1)**. Phase 1 established that each dataset is one clinical sample from one donor. There is no donor contributing two samples, and no sample split across batches. |
| `sample_id` ↔ technical batch | **Unknown but assumed identical**. No technical batch variable was ever recorded. Library chemistry, capture channel, sequencing run and operator are unmeasured (Phase 1 M8 Q4). |
| `sample_id` ↔ biological condition | **No biological condition variable exists.** There is no treatment arm, no disease stage, no site, no NF1 status. There is therefore no condition-level contrast to protect and none to destroy. |

### 4.1 The core confounding statement

> Dataset = sample = patient = (presumed) technical batch. These four concepts are a
> single 4-level variable. Correcting for technical batch is mathematically
> indistinguishable from erasing between-patient tumour biology. No design, no
> covariate and no statistical method available in this dataset can separate them.

This is not a defect introduced by Phase 2. It is an inherent property of a four-sample,
one-sample-per-patient study without recorded technical metadata, and it was already
documented in `reports/PRE_INTEGRATION_ASSESSMENT.md` §5 and `PHASE1_HANDOFF.md` §7.

---

## 5. Is integration justified at all?

Two facts argue that some correction is needed for the intended downstream goal:

1. **Pre-integration structure is essentially total.** In the frozen Phase 1 shared PCA
   space, the mean same-dataset fraction among the 15 nearest neighbours is **0.9613**
   and the mean neighbourhood dataset entropy is **0.0966**. Cells sit almost
   exclusively next to cells from their own sample. Joint clustering of the
   non-integrated space would largely reproduce four per-sample clusterings, which
   Phase 1 has already produced independently and better.
2. **The samples do appear to share major compartments.** *(Legacy evidence, provenance
   only — see the caveat below.)* The legacy `orig.anno` labels place malignant,
   macrophage, fibroblast, T-cell and endothelial populations in all four samples:

   | Legacy label *(provenance only)* | Total | M1 | M2 | M3 | M4 |
   | --- | ---: | ---: | ---: | ---: | ---: |
   | Malignant NC-MES or MES-like | 7,607 | 2,670 | 707 | 176 | 4,054 |
   | Macrophages | 4,068 | 1,374 | 573 | 921 | 1,200 |
   | Malignant SCP-like | 1,618 | 1,571 | 23 | 14 | 10 |
   | Fibroblast | 1,355 | 286 | 605 | 24 | 440 |
   | T cells | 1,075 | 229 | 52 | 561 | 233 |
   | Immature SC-like | 997 | 761 | 57 | 112 | 67 |
   | Endothelial | 861 | 313 | 112 | 151 | 285 |
   | B cells | 813 | 19 | 4 | **776** | 14 |
   | VSMCs/pericytes | 460 | 265 | 42 | 71 | 82 |
   | Cycling tumor cells | 367 | 51 | 8 | 7 | 301 |

   > **Caveat.** These are *legacy* labels from the original analysis. Under Phase 2 §6
   > they are inventoried for provenance and risk assessment only. They are **not** used
   > to choose the Harmony variable, dimensions, resolution, markers or annotations, and
   > Phase 2 will re-derive every compartment independently in M13–M15.

The shared immune and stromal compartments are the populations for which cross-sample
alignment is scientifically meaningful. The malignant compartment is the population for
which alignment is scientifically dangerous (see §6.1).

---

## 6. Risk assessment

### 6.1 Risk of over-correction — **HIGH, and unavoidable in principle**

- **Patient-private malignant programs will be pushed together.** Harmony's objective is
  to make each cell's neighbourhood batch-diverse. Malignant cells are, by construction,
  patient-private (distinct driver mutations, distinct copy-number landscapes, distinct
  subclonal programs). Forcing MPNST_1 tumour cells to neighbour MPNST_4 tumour cells is
  exactly what the algorithm is asked to do, and exactly what is biologically wrong.
- **Sample-restricted populations are the most fragile.** The legacy labels indicate at
  least three strongly sample-restricted populations: B cells (95% in MPNST_3),
  Malignant SCP-like (97% in MPNST_1), and Cycling tumour cells (82% in MPNST_4). A
  population present in one batch has no counterpart to be mixed with; Harmony can
  distort or dissolve such populations into the nearest abundant cluster.
- **The `MPNST_4` stress/high-mitochondrial state** (flagged in Phase 1 as cluster C09
  with percent.mt R² = 0.47) is a documented candidate for erasure.
- **No independent ground truth is available** to prove biology was preserved, because
  the only cell-type labels in the object are legacy ones that Phase 2 must not treat as
  truth.

**Mitigation:** the non-integrated Phase 1 baseline (`pca`, `umap_preintegration`) is
preserved untouched, M12 computes explicit over-correction diagnostics (pre→post
nearest-neighbour retention, per-cluster batch composition, silhouette of biological
structure), and the malignant compartment is annotated conservatively in M15/M16 with
explicit uncertainty labels.

### 6.2 Risk of under-correction — **MODERATE**

Starting from 96.1% same-dataset neighbours, Harmony with default `theta = 2` may leave
substantial residual sample structure. If it does, the correct response is **not** to
raise `theta` until the UMAP looks mixed. It is to report the residual structure, and to
propose a documented sensitivity analysis for researcher approval (Phase 2 §12).

### 6.3 Unresolved confounding — **TOTAL and permanent**

Anatomical site, primary vs. metastatic status, treatment history, NF1 germline status,
tumour grade, donor age and sex, library chemistry, capture channel, sequencing run and
operator are all unrecorded. Any biological difference attributable to any of these is
irreversibly mixed into the `sample_id` term that Harmony will remove.

**Consequence to be carried into every later Phase 2 document:** after integration, no
statement of the form "population X differs between tumours" can be supported, because
the between-tumour axis is the axis that was deliberately removed.

---

## 7. Decision

**Proposed Harmony grouping variable: `sample_id` (4 levels), used as the sole entry in
`group.by.vars`.**

Rationale:

1. It is the only variable in the object that carries the technical batch structure at all.
2. It is cell-for-cell identical to `orig.ident`, so the choice between them is cosmetic;
   `sample_id` is the canonical Phase 1 splitter (`config/config.yaml: dataset_id_column`)
   and keeps Phase 1 and Phase 2 provenance aligned.
3. Every alternative is either prohibited (legacy integration/annotation output),
   degenerate (zero-variance QC flags), continuous (not a Harmony grouping variable), or
   non-existent (all clinical and technical covariates).
4. No second grouping variable is available, so no multi-variable Harmony design is
   possible and none is proposed.

**This decision is proposed, not taken.** Because §4.1 establishes perfect
batch–biology confounding, Phase 2 §8 requires the researcher to review and approve it
before any Harmony execution.

### 7.1 The three options actually open to the researcher

| Option | What it means | Consequence |
| --- | --- | --- |
| **A. Approve (recommended)** | Run default Harmony on `sample_id`, keep the Phase 1 baseline as the permanent reference, and interpret the integrated space as a *joint annotation space* rather than as evidence about between-tumour differences. | Enables unified cell-type annotation across all four tumours. Between-patient malignant heterogeneity is knowingly suppressed in the integrated embedding and must be studied in the preserved non-integrated baseline instead. |
| **B. Approve with a restricted scope** | Same run, but the researcher additionally directs that the malignant compartment be interpreted only in the non-integrated space. | Safer for tumour biology; adds a Phase 3 workstream. |
| **C. Decline integration** | Continue Phase 2 on the non-integrated baseline only. | Preserves all between-tumour biology but leaves joint clustering dominated by sample identity, and makes cross-tumour cell-type annotation substantially harder. |

The recommendation is **Option A**, with the M12 over-correction diagnostics treated as a
genuine gate rather than a formality: if M12 shows that shared immune/stromal
compartments failed to align *or* that sample-restricted populations were dissolved,
default Harmony will be reported as inadequate and a sensitivity analysis proposed,
rather than silently accepted.

---

## 8. Exact Harmony configuration proposed for M11

Defaults of the installed package (`harmony 1.2.4`) are used for every scientific
parameter. Only technical arguments — the object, which reduction to correct, which
dimensions, where to store the result, and the grouping variable — are specified.

```r
set.seed(42)

obj <- harmony::RunHarmony(
  object         = obj,                 # Phase 1 handoff object, loaded read-only
  group.by.vars  = "sample_id",         # the decision made in this document
  reduction.use  = "pca",               # frozen Phase 1 pre-Harmony PCA (assay = SCT)
  dims.use       = 1:30,                # matches the frozen Phase 1 baseline exactly
  reduction.save = "postint_harmony",   # Phase 1 -> Phase 2 namespace contract
  verbose        = TRUE                 # convergence messages captured in the SLURM log
)
```

### 8.1 Parameters explicitly set vs. left at package defaults

| Parameter | Value | Explicit or default | Note |
| --- | --- | --- | --- |
| `group.by.vars` | `"sample_id"` | explicit | This decision. |
| `reduction.use` | `"pca"` | explicit | Required technical argument. |
| `dims.use` | `1:30` | explicit | Required technical argument; see §8.2. |
| `reduction.save` | `"postint_harmony"` | explicit | Required by the Phase 1 namespace contract. |
| `verbose` | `TRUE` | explicit | Diagnostics only; no effect on the result. |
| `project.dim` | `TRUE` | **default** | |
| `theta` | `NULL` → 2 per variable | **default** | Diversity clustering penalty. Not tuned. |
| `sigma` | `0.1` | **default** | Soft-clustering width. |
| `lambda` | `1` | **default** | Ridge regression penalty. |
| `nclust` | `NULL` (auto) | **default** | Number of soft clusters. |
| `max_iter` | `10` | **default** | Harmony outer iterations. |
| `early_stop` | `TRUE` | **default** | |
| `ncores` | `1` | **default** | |
| `plot_convergence` | `FALSE` | **default** | Objective history captured numerically instead — see §8.3. |
| `.options` | `harmony_options()` | **default** | `alpha = 0.2`, `tau = 0`, `block.size = 0.05`, `max.iter.cluster = 20`, `epsilon.cluster = 1e-3`, `epsilon.harmony = 1e-2`. |
| Random seed | `42` | explicit | Inherited from Phase 1 `config/config.yaml`. |

### 8.2 Why `dims.use = 1:30`

`dims.use` has no package default (`NULL` means "all dimensions of the reduction",
i.e. all 50 computed PCs), so a value must be chosen. `1:30` is proposed because:

- Phase 1 M8 built the frozen baseline neighbour graph and UMAP on **exactly** PCs 1–30.
  Using the same 30 dimensions makes the pre/post comparison in M12 a comparison of
  *integration*, not a comparison of dimensionality.
- PCs 1–30 capture **91.1%** of the variance in the shared PCA
  (`reports/combined/pre_integration/pca_variance_explained.tsv`); PC 30 explains 0.61%
  and PCs 31–50 each explain <0.61%, i.e. the noise floor.
- The Phase 1 *per-dataset* PC recommendations (5–9 PCs) are not applicable: they were
  derived from four separate single-sample PCAs, whereas this is one shared 19,716-cell
  PCA that must span four tumours' worth of cell states.

### 8.3 Convergence diagnostics

`RunHarmony.Seurat` returns a Seurat object and discards Harmony's internal objective
history. To satisfy Phase 2 M11 item 7 without altering the primary result, M11 will
additionally call the matrix-level entry point on the *same* input with the *same*
defaults and `return_object = TRUE`, extract `objective_harmony` / `objective_kmeans`
per iteration into `results/phase2/harmony/harmony_convergence.tsv`, and then **assert
that the resulting corrected embedding is numerically identical** to the one stored in
`postint_harmony`. If the two disagree, M11 fails loudly rather than reporting a
mismatched diagnostic.

---

## 9. What this decision does NOT authorise

- It does not authorise tuning `theta`, `lambda`, `sigma`, `nclust` or `max_iter`.
- It does not authorise a second grouping variable.
- It does not authorise using the integrated embedding for any between-tumour,
  between-condition or differential-abundance claim.
- It does not authorise overwriting, replacing or regenerating any Phase 1 object,
  reduction, graph or metadata column.

---

## 10. Researcher sign-off

| Field | Value |
| --- | --- |
| Proposed grouping variable | `sample_id` |
| Proposed dimensions | `1:30` of reduction `pca` |
| Proposed output reduction | `postint_harmony` |
| Harmony version | 1.2.4 |
| Scientific parameters | all package defaults |
| Decision status | ☐ Approved ☐ Approved with scope restriction (Option B) ☐ Declined (Option C) |
| Approved by | _________________ |
| Date | _________________ |
