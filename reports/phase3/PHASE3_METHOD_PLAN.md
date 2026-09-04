# Phase 3 Method Plan

**Phase 3 · Milestone M18 · MPNST tumour–immune cell–cell communication**
*Generated: 2026-09-03*

---

## 1. Biological objective

Determine how MPNST tumour cells interact with immune and other tumour-microenvironment
populations, identify recurrent and biologically plausible signalling circuits, establish
which findings are concordant across independent CCC frameworks, and connect candidate
tumour-derived signals to receiver-cell transcriptional states.

**This is not a methods-benchmarking study.** Multiple tools are used because CCC inference
is method-dependent; an interaction supported by several independent frameworks earns higher
confidence. Method *disagreement* is preserved and reported, never averaged away.

---

## 2. Input — verified from repository evidence

| Property | Value |
| --- | --- |
| Object | `results/phase2/phase2_final_object.rds` |
| md5 | `153d5f6acc70f9c05aa48cabc4f4ac2d` (verified at M19 runtime) |
| Cells | 19,716 · Features `RNA` 31,764 / `SCT` 29,113 |
| CCC grouping field | **`annotation_ccc`** (17 identities) |
| Detailed annotation | `postint_celltype_level1/2/3_refined` (retained; source of truth) |
| Sample field | **`sample_id`** — simultaneously the **dataset** and the **patient** |
| Condition field | **does not exist** (audited across all 93 metadata columns in Phase 2 M10) |
| Tumour population | `MPNST-Tumor`, 3,420 cells (17.35%) |

---

## 3. Populations analysed

The three populations flagged not CCC-ready in Phase 2 are **excluded outright**:
`Candidate-Malignant-Unresolved` (1,231), `Uncertain` (720),
`Low-quality-excluded` (438) — 2,389 cells, 12.1%.

That leaves **14 populations / 17,327 cells**, and the M19 audit found **every one of them
evaluable in all four samples**:

| Compartment | Populations |
| --- | --- |
| Tumour | `MPNST-Tumor` |
| Immune — myeloid | `Macrophage`, `Monocyte`, `Dendritic`, `Plasmacytoid-DC` |
| Immune — lymphoid | `CD8-T`, `CD4-T`, `NK`, `T-cell-other`, `B-cell`, `Plasma-cell` |
| Stromal | `Fibroblast`, `Pericyte-VSMC` |
| Vascular | `Endothelial` |

No Mast population exists at Harmony resolution 1.0 (only ~40 cells by programme score,
never resolving into a cluster), so no Mast analysis is possible and none is fabricated.

---

## 4. Directional priorities

| Priority | Direction | Pairs |
| --- | --- | ---: |
| **P1** | `MPNST-Tumor → TME` | 13 |
| **P1** | `TME → MPNST-Tumor` | 13 |
| P2 | `TME ↔ TME` | 156 |

All 182 directed pairs are evaluable in all four samples. P2 pairs are retained but reported
separately so they cannot overwhelm the tumour-centric analysis.

---

## 5. Expression basis — a documented decision

CCC uses the **`RNA` assay, `JoinLayers()` + `LogNormalize` (scale.factor 1e4)**, *not* SCT.

*Rationale:* the `SCT` assay carries **four SCTransform models** (Phase 1 M8 ran SCTransform
layer-wise), so SCT values are not normalised on a common footing across samples. Every CCC
framework expects one uniformly log-normalised matrix. The RNA layers were split four ways
and were joined for this purpose.

**Harmony coordinates are never used as expression.** Harmony is embedding-level; it informs
the neighbourhood structure used by LochNESS, never the ligand–receptor arithmetic.

---

## 6. Minimum-cell policy — stated, not silently chosen

1. Populations flagged not CCC-ready in Phase 2 are excluded outright.
2. A population is **EVALUABLE in a sample** only with **≥ 10 cells** in that sample.
   *Basis:* 10 is CellChat's `min.cells` default and the conventional floor for per-group
   mean-expression estimates in LR inference.
3. A sample is **ANALYSABLE** only if ≥ 3 populations are evaluable in it. All four qualify.
4. **A population below the minimum is recorded NOT EVALUABLE, never as "no signalling".**
   Absence of evidence is not evidence of absence.

Recorded in `results/phase3/tables/CCC_SAMPLE_CELL_COUNTS.tsv`,
`CCC_POPULATION_EVALUABILITY.tsv`, `CCC_PAIR_EVALUABILITY.tsv`.

---

## 7. Methods and their distinct roles

| Framework | Role in Phase 3 | Output used |
| --- | --- | --- |
| **LIANA** | Primary multi-resource / multi-score **consensus LR** framework | Aggregate rank per `sender·receiver·ligand·receptor` |
| **CellChat** | Independent **pathway- and network-oriented** framework | Pathway-level communication, outgoing/incoming signalling roles, information flow |
| **CellPhoneDB** | Independent **permutation-based LR** framework | Significant means and permutation p-values |
| **NicheNet** | **Not** a fourth LR list generator — it models **sender ligand → receiver transcriptional response** | Ligand activity (AUPR) against receiver gene sets, ligand→target regulatory potential |

**Raw scores are never averaged across tools.** A CellChat probability, a CellPhoneDB
p-value, a LIANA aggregate rank and a NicheNet ligand activity are different quantities on
different scales. Each contributes a **method-specific boolean support flag**, and those
flags are combined into a transparent concordance class.

---

## 8. Concordance model

Canonical key: `sender | receiver | ligand | receptor` (plus `pathway` where a framework
provides one). Gene symbols and multi-subunit complexes are normalised before joining.

Per interaction: `supported_LIANA`, `supported_CellChat`, `supported_CellPhoneDB`,
`n_LR_methods_supported`, `concordance_class`.

With three independent LR frameworks:

| Class | Rule |
| --- | --- |
| **High concordance** | 3/3 LR frameworks |
| **Moderate concordance** | 2/3 |
| **Single-method** | 1/3 |
| **Discordant/ambiguous** | supported by a framework but contradicted by another where the pair was testable |

**NicheNet is represented as orthogonal receiver-response support**
(`NicheNet_ligand_support`), *not* counted as a fourth equivalent LR framework, because it
supports a ligand→receiver-programme claim rather than the identical receptor pair.

If a framework fails to install or run, the denominator is adjusted and the change is stated
explicitly in the concordance report rather than silently reducing the thresholds.

---

## 9. Biological prioritisation

Priority combines, without collapsing into an unjustified composite score:

```
method concordance + sample recurrence + patient recurrence
+ sender ligand expression + receiver receptor expression
+ receiver-response evidence + literature support
```

An **evidence matrix** is preferred over a single number. Where an ordering is needed it is
lexicographic on (concordance class, sample recurrence, receiver-response support), and the
rule is printed in the report.

---

## 10. LochNESS — role and adaptation

**LochNESS is not a ligand–receptor method.** Phase 3 will never write or imply
"high LochNESS = communication". It is used as an **orthogonal transcriptional-neighbourhood
phenotype** characterising receiver-state structure *after* CCC nominates candidate
signalling:

```
CCC -> candidate sender->receiver signalling -> receiver-response analysis
    -> LochNESS receiver-state structure -> additional evidence
```

Implementation follows the **official MMCA R formulation** (`k = round(0.5·√N)`,
L2-normalised PCA, exact `FNN::get.knnx`, disjoint query/reference sets so same-sample
neighbours are excluded, global fraction over the reference set, label-permutation null with
the exclusion structure preserved). The perturb-seq Python implementation is retained as an
independent cross-check. See `LOCHNESS_IMPLEMENTATION_AUDIT.md`; the MPNST-specific target /
reference / null design is written in `LOCHNESS_MPNST_DESIGN.md` **before** any LochNESS run.

---

## 11. Statistics

No cell-level pseudoreplication. Emphasis on **sample recurrence, patient recurrence, effect
size, concordance and receiver response** rather than p-values computed over 17,327 cells.
Where a stratified count comparison is genuinely warranted, sample is used as the stratum.
Permutation nulls are **sample-aware**.

---

## 12. Milestone plan

| Milestone | Content |
| --- | --- |
| M18 | Reconstruction, perturb-seq and MMCA audits, environment, method plan |
| M19 | Sample-aware CCC input preparation and minimum-cell policy ✅ |
| M20 | LIANA, CellChat, CellPhoneDB execution (independent SLURM jobs) |
| M21 | Concordance: canonical keys, support flags, classes, method overlap |
| M22 | Tumour ↔ immune/TME biological prioritisation + focused literature |
| M23 | NicheNet receiver-response analysis |
| M24 | LochNESS design, R-vs-Python comparison, receiver-state runs |
| M25 | CCC + LochNESS integration as an evidence matrix |
| M26 | Robustness: sample/patient consistency, leave-one-sample-out, method disagreement |
| M27 | Final freeze, manifest, handoff, Phase 4 spatial recommendations |

---

## 13. Prohibitions carried through Phase 3

No survival analysis · no treatment-response modelling · no deep learning or Transformers ·
no trajectory inference · no RNA velocity · **no CNV inference** · **no spatial inference** ·
no pseudobulk condition DE unrelated to receiver-response questions · no large unrelated
pathway screens · no arbitrary composite scores without stated justification.

---

## 14. Interpretation language

Conclusions use **predicted / inferred / supported / consistent with / concordant /
associated with**. Every report states:

> scRNA-seq ligand–receptor analysis infers communication *potential* from expression and
> associated transcriptional programmes. It does not establish physical cell adjacency or
> direct signalling.
