# Receiver-Response Report (NicheNet)

**Phase 3 · Milestone M23 · MPNST** · *Generated 2026-09-03 · SLURM JobID 19896063*

> NicheNet models **regulatory potential** from a prior network. It does not prove physical
> signalling, and its AUPR is never averaged with LR scores.

---

## 1. What was asked

Not "which ligand–receptor pairs exist" — the LR frameworks answered that. NicheNet answers
**Question 5**: *does a tumour-derived candidate ligand explain the receiver's own
transcriptional programme?*

```
MPNST-Tumor -> candidate ligand -> receiver receptor -> receiver target genes
```

| Design element | Choice |
| --- | --- |
| Prior model | Zenodo record 7074291; ligand–target matrix 33,354 × 1,226; lr_network 4,986 rows |
| Sender | `MPNST-Tumor`, ligands detected in ≥10% of tumour cells |
| Gene set of interest | **receiver population marker genes** (receiver vs all other CCC populations, RNA LogNormalize, Wilcoxon, adj p < 0.05, only.pos) — cluster characterisation, **not** condition-level DE |
| Background | genes detected in ≥10% of receiver cells and present in the prior matrix |
| Potential ligands | expressed by tumour **and** with a receptor expressed in the receiver |
| Score | `aupr_corrected` |

---

## 2. Results — a clean split by lineage

| Receiver | Cells | Gene set | Ligands tested | Top ligand | AUPR | Top-20 also LR-prioritised |
| --- | ---: | ---: | ---: | --- | ---: | ---: |
| Macrophage | 3,065 | 1,514 | 229 | **CSF1** | 0.166 | 14/20 |
| Monocyte | 190 | 969 | 224 | **CSF1** | 0.123 | 15/20 |
| Dendritic | 612 | 1,342 | 227 | **CSF1** | 0.107 | 15/20 |
| CD8-T | 382 | 539 | 141 | **IL15** | 0.116 | 15/20 |
| NK | 106 | 278 | 138 | **IL15** | 0.089 | 16/20 |
| CD4-T | 606 | 604 | 141 | **IL15** | 0.084 | 16/20 |
| Endothelial | 960 | 2,941 | 264 | **TGFB1** | 0.067 | 15/20 |
| Fibroblast | 5,064 | 2,941 | 291 | TGFB1 | **0.022** | 16/20 |

### 2.1 Myeloid receivers converge on CSF1

CSF1 is the top-ranked ligand for **all three** myeloid receivers, independently of the LR
analysis. This is a genuinely orthogonal convergence: NicheNet asks whether CSF1 explains the
macrophage/monocyte/dendritic *gene programme*, while the LR frameworks ask whether
CSF1 and CSF1R are co-expressed in the right cells. Both say yes. It is also the axis with
existing clinical traction in neurofibroma (pexidartinib).

Runners-up for macrophage — CCL2 (0.160), IL15 (0.138), CCL3 (0.131), TGFB1 (0.118) — are a
coherent monocyte-recruitment and polarisation set.

### 2.2 Lymphoid receivers converge on IL15

IL15 is the top ligand for CD8-T, NK **and** CD4-T. IL15 is a survival and homeostatic
cytokine for CD8 and NK cells, so a tumour-derived IL15 signal would be expected to *support*
lymphocyte persistence. Read against §2.3 of the interaction report — where the tumour also
engages NKG2A and LILRB1/2 — the lymphoid picture is genuinely mixed rather than uniformly
immunosuppressive.

### 2.3 Endothelial

TGFB1 (0.067), HMGB1 (0.059), VEGFA (0.059), ANGPT1 (0.044) — a coherent angiogenic set,
converging with the LR finding of VEGFA → KDR/FLT1/NRP1 at 4 evidence streams.

### 2.4 Fibroblast — a clear negative

The best ligand reaches **AUPR 0.022**, and the remainder are ≈0 or negative. **Tumour-derived
ligands do not explain the fibroblast transcriptional programme.** This is reported as a
negative result, not omitted. It is also the most interesting negative in Phase 3: the
tumour↔fibroblast LR signal is abundant, but there is no accompanying evidence that tumour
ligands drive the fibroblast state. Two readings are possible and cannot be separated here —
either the LR co-expression is not functional, or the fibroblast programme is dominated by
inputs other than tumour ligands.

---

## 3. Divergence from the LR analysis, reported not reconciled

**APP does not appear among the top NicheNet ligands for macrophages**, despite APP → CD74
being the single strongest LR finding (3 frameworks, 4/4 patients, 5 receivers).

The two methods are asking different questions — LR co-expression versus explanation of a
downstream gene programme — so this is not a contradiction, but it is a real limit on the
APP-CD74 claim: **we have strong co-expression evidence and no receiver-programme evidence.**
The correct statement is "predicted APP–CD74 engagement", not "APP drives the macrophage
state".

Conversely CSF1 is top by NicheNet for all myeloid receivers but only reaches 2-framework
support in the LR analysis. The two evidence streams nominate **different** leading
candidates, and both are reported.

---

## 4. Limitations

1. NicheNet's prior network is human-generic, not MPNST-specific.
2. The gene set of interest is receiver *identity* markers, so the analysis asks what
   explains the receiver's steady-state programme, not its response to a perturbation.
3. AUPR values are low in absolute terms (0.02–0.17), as is usual for this framework; only
   the ranking is interpreted.
4. NK (106 cells) and Monocyte (190 cells) give small gene sets and unstable rankings.
5. No sample-stratified NicheNet was run — the gene sets are cohort-level, so this analysis is
   not sample-aware in the way the LR analysis is.

## 5. Files

`results/phase3/receiver_response/NicheNet_ligand_activity.tsv` ·
`NicheNet_receiver_summary.tsv` · `receiver_target_programs.tsv` ·
`NicheNet_receiver_response_support.tsv` · figures `results/phase3/figures/M23/`.
