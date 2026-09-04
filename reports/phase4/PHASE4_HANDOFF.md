# MPNST Phase 4 Handoff — SCEVAN Malignancy Refinement, Tumour-State Resolution, Targeted CCC Reassessment

*Milestones M28–M35 · executed and frozen 2026-09-03 · readable without conversation history*

> **SCEVAN infers copy number from gene-expression patterns. This is not DNA sequencing.**
> Inference is reliable for broad chromosomal, arm-level and large-segment events only.
> **Single-gene CNV calls are not asserted anywhere in Phase 4.** A SCEVAN non-malignant call is
> **not proof** of non-malignancy: some MPNST malignant cells may be copy-number quiet.
> **n = 4 patients.** No population-level or epidemiological claim follows from any result here.

---

## 1. Objective

Phase 2 defined the malignant compartment **conservatively and deliberately**: `MPNST-Tumor`
= 3,420 cells, **17.35%** of 19,716, collapsed only where positive MPNST/Schwann/neural-crest
marker evidence existed. That was correct on marker evidence alone — the alternative was the
prohibited `if (!immune) MPNST-Tumor` inference — but it left 5,064 `Fibroblast`, 1,231
`Candidate-Malignant-Unresolved`, 720 `Uncertain` and 435 `Pericyte-VSMC` cells whose malignant
status **marker expression cannot settle**, and that uncertainty propagated into every
tumour-centric Phase 3 interaction.

Phase 4 brought **orthogonal evidence** — inferred copy number — to settle malignant identity,
then resolved tumour states and clones, then asked whether the Phase 3 architecture survives.

**Phase 4 does not invalidate Phase 2 or Phase 3.** Phase 2 annotation is historical biological
annotation; Phase 4 malignancy is additional orthogonal evidence. Phase 4 **only adds** metadata
fields, and that is asserted, not assumed (§23 below).

## 2. Phase 2 / Phase 3 input

| Item | Value |
| --- | --- |
| Phase 2 object | `results/phase2/phase2_final_object.rds`, md5 `153d5f6acc70f9c05aa48cabc4f4ac2d` — **verified equal** to the Phase 2 manifest at every stage |
| cells / features | 19,716 / 31,764 RNA, 29,113 SCT |
| Phase 3 baseline | `results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv` (36,486 supported) and `results/phase3/tables/final/CCC_CONCORDANCE.tsv` (502,846 keys) |
| Phase 3 CCC object | `results/phase3/ccc/ccc_input_object.rds`, md5 `aeb02c59df8a42a24074df509c662f7f` |
| input assay / layer | **RNA**, the four per-sample `counts` layers — integer-verified; the Phase 2 default `SCT` assay carries four SCTransform models and was deliberately unused |
| gene identifiers | **symbols**; 28,340 of 31,764 (89.2%) map to SCEVAN's chr1–22 annotation |

## 3. Why SCEVAN

De Falco, Caruso, Sudmant & Ceccarelli, *A variational algorithm to detect the clonal copy
number substructure of tumors from scRNA-seq data*, **Nature Communications 14, 1074 (2023)**,
doi:10.1038/s41467-023-36790-9 · https://github.com/AntonioDeFalco/SCEVAN

Chosen because it starts from raw counts, identifies confident normal cells **automatically**
(which is what makes a non-circular primary run possible at all), separates malignant from
non-malignant cells, infers large-scale CNA profiles, characterises clones and subclones,
supports multi-sample clonal comparison, and is substantially cheaper than inferCNV in the
datasets the original study evaluated. It is also native R.

**Not claimed to be universally superior to inferCNV.** No second CNV method was added for
benchmarking, and none was needed: SCEVAN did not fail in a way requiring independent CNV
validation (§66).

## 4. SCEVAN implementation

Version **1.0.3** (package date 2025-02-12) with **yaGST 2017.8.25**. Both **already installed**,
so no R dependency was resolved, upgraded or downgraded for the CNV work.

Per patient, twice, with identical parameters:

```r
SCEVAN::pipelineCNA(count_mtx = <sparse integer RNA counts for that patient>,
  sample = "<PATIENT>_<primary|sensitivity>", par_cores = 8,
  norm_cell = NULL | <immune cells>, SUBCLONES = TRUE, beta_vega = 0.5,
  ClonalCN = TRUE, plotTree = TRUE, organism = "human",
  ngenes_chr = 5, perc_genes = 10, FIXED_NORMAL_CELLS = FALSE, output_dir = "./output")
```

Seed 42 throughout (SCEVAN, clustering, UMAP, module scores, subsampling). `par_cores` matched
`--cpus-per-task` exactly — no oversubscription.

### Three SCEVAN 1.0.3 defects, worked around without patching the package

1. **`output_dir` not honoured by the read-back/plot helpers.** `getScevanCNV`,
   `getScevanCNVfinal`, `plotAllClonalCN`, `plotAllSubclonalCN`, `plotConsensusCNA` and
   `analyzeSegm2` hardcode `path = "./output"`, and `plotCNclonal()` does not forward it.
   Classification completed, then plotting died with `cannot open the connection`
   (SLURM 19896576). **Fix**: each run sets its own working directory and lets `output_dir` take
   SCEVAN's default, so both paths coincide.
2. **`plotTSNE` requires Python `umap-learn`** via reticulate, and it runs **before** the line
   writing subclone labels into `classDf` — so it is not optional and skipping it would lose the
   clone assignments (SLURM 19896623). **Fix**: umap-learn 0.5.12 in an isolated `p4_umap_env`
   reached via `RETICULATE_PYTHON`; `R_env` verified unchanged before and after.
3. **`plotCloneTree` fails under ggplot2 4.x** — ggtree calls the removed
   `ggplot2:::is.waive()`. SCEVAN catches this internally, so all eight runs returned normally
   and **every clone assignment, CN profile, CNA matrix and onco-heatmap is intact**; only the
   clone **phylogeny plot** is absent. **ggplot2 was deliberately not downgraded** — that would
   destabilise the Phase 3 figure suite for a dendrogram that is not a deliverable. The
   clone/state relationship is reported as an explicit cross-tabulation instead.

## 5. Normal-reference design, and the anti-circularity constraint

Two circular arguments were available and both are prohibited (§16): declaring the disputed
fibroblast-like cells normal and then reporting that SCEVAN "confirms" it; or declaring all
`MPNST-Tumor` cells malignant and then reporting that SCEVAN "independently confirms" them.

| Group | Populations | Cells | Role |
| --- | --- | ---: | --- |
| High-confidence immune | CD4-T, CD8-T, NK, T-cell-other, B-cell, Plasma-cell, Macrophage, Monocyte, Dendritic, Plasmacytoid-DC | 7,448 | sensitivity run only |
| Endothelial | Endothelial | 960 | **excluded** from the reference set |
| **Disputed** | MPNST-Tumor, Candidate-Malignant-Unresolved, Fibroblast, Pericyte-VSMC, Uncertain, Low-quality-excluded | **11,308** | **never a fixed reference — the question** |

* **Primary run: `norm_cell = NULL`** — SCEVAN's confident-normal detection ran with **no input
  from Phase 2 labels at all**. This is the run that carries inferential weight.
* **Sensitivity run**: the immune cells as `norm_cell`, `FIXED_NORMAL_CELLS = FALSE`.
* Endothelium kept out on purpose: it is a Phase 3 receiver of interest (VEGFA→KDR/FLT1/NRP1,
  the JAG/NOTCH circuit), so its malignancy call had to stay independent.

**`FIXED_NORMAL_CELLS = TRUE` is prohibited in this project.** The installed source executes
`cellType_pred[!cellType_pred %in% norm_cell_names] <- "malignant"` — forcing every
non-reference cell to malignant, which is precisely the banned inference arriving through a
function argument.

**Decision rules were written and committed before any SCEVAN result was inspected**
(`MALIGNANCY_DECISION_RULES.md`, 20 ordered rules, every threshold fixed a priori). Two
amendments were added afterwards and **both are dated and labelled post-hoc** — see §10.

## 6. Per-patient malignant classification

| Patient | cells | assessed | **SCEVAN malignant** | non-malignant | filtered | subclones | confident normals | primary↔sensitivity agreement |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| MPNST_1 | 7,615 | 7,043 | **1,952** | 5,091 | 572 | 7 | 30 | **0.9951** |
| MPNST_2 | 2,284 | 2,205 | **694** | 1,511 | 79 | 4 | 30 | **0.9663** |
| MPNST_3 | 2,940 | 2,749 | 1,414 | 1,335 | 191 | 3 | 25 | **0.0864** |
| MPNST_4 | 6,877 | 6,452 | **4,734** | 1,718 | 425 | 8 | 30 | **0.9974** |
| **total** | **19,716** | **18,449** | **8,794** | 9,655 | 1,267 | 22 | 115 | — |

`filtered` cells are recorded as **not assessed**, never as non-malignant.

**MPNST_3 failed a basic sanity check.** Both its runs assign large blocks of unambiguous immune
cells to the malignant class in **opposite directions** — the primary run calls 89/89 B cells,
327/328 CD4-T, 135/136 CD8-T, 32/32 NK and 165/165 pDC malignant while calling macrophages
normal; the sensitivity run inverts it. Six of nine testable immune populations individually
exceed 25% malignant, only 25 confident normal cells were found, MPNST_3 has the lowest depth of
the four (median 1,594 genes/cell), and its three "clones" are pure immune lineages. Malignant
**promotions** from MPNST_3 are disabled. **Cost stated plainly: any Mes-NC-like malignant
fibroblast population in MPNST_3 cannot be detected by this analysis.**

**MPNST_4's high overall immune rate (0.293) is a different phenomenon** — 98.2% of its plasma
cells, an immunoglobulin-locus artefact (§10) — with 0 of 9 populations failing and the highest
agreement of the four.

## 7. CNA results

59 broad clonal segments (≥10 Mb, non-neutral). Recurrence, reported **twice** because MPNST_3
is unreliable:

| Chromosome | direction | reliable 3 | all 4 |
| --- | --- | :--: | :--: |
| **18** | **loss** | **3/3** | 4/4 |
| **2** | **gain** | **3/3** | 3/4 |
| **7** | **gain** | **3/3** | 3/4 |
| 3, 5, 6, 11, 15, 19, 22 | mixed | 2/3 | 2/4 |

Pairwise Jaccard of broad-event sets among the reliable patients is only **0.20–0.23**, and only
**1 of 59** events appears in all four — so CNV architecture is **predominantly patient-specific
with a small shared core**. These are **shared or recurrent CNA patterns**, never "the same
clone" across patients (§50).

**NF1 and NF2 loci, reported with deliberate caution.** MPNST_4 carries a **38.1 Mb CN = 1
segment on chr17 (0.44–38.51 Mb) whose interval contains the NF1 locus (17q11.2)**, with gain
across the rest of chromosome 17. MPNST_1 and MPNST_4 each carry a **34.2 Mb CN = 1 segment on
22q containing NF2 (22q12.2)**. Consistent with a nerve-sheath malignancy — but this analysis
**cannot** distinguish a focal deletion from large-segment loss that happens to include the
locus, and does not attempt to. Confirming NF1/NF2 status requires DNA sequencing.

## 8. Phase 2 versus SCEVAN

Malignant fraction **of assessed cells**, primary run:

| Phase 2 `annotation_ccc` | MPNST_1 | MPNST_2 | MPNST_3 † | MPNST_4 |
| --- | ---: | ---: | ---: | ---: |
| **Fibroblast** | **0.628** | **0.559** | 0.078 | **1.000** |
| **Candidate-Malignant-Unresolved** | **0.733** | 0.167 | 0.789 | 0.429 |
| MPNST-Tumor | 0.225 | 0.089 | 0.429 | 0.912 |
| Pericyte-VSMC | 0.016 | 0.071 | 0.081 | 0.942 |
| Endothelial | 0.000 | 0.009 | 0.046 | 0.527 |
| Uncertain | 0.008 | 0.056 | 0.438 | 0.933 |
| CD4-T / CD8-T / NK / T-other | 0.000 | 0.000 | ~0.99 | ≤0.031 |
| Macrophage / Monocyte / Dendritic | 0.000 | 0.000 | ≤0.055 | 0.000 |
| Plasma-cell | 0.197 | 0.635 | 0.993 | 0.982 |

† sanity-gate failure; shown for completeness, carries no weight.

## 9. The ambiguous populations

**`Fibroblast` — mixed, and predominantly malignant.** Of 5,064: **4,036 refined Malignant
(79.7%)**, 908 Non-malignant (17.9%), 120 Ambiguous. Four independent legs: majority-malignant
in **3 of 4** patients; **threshold-independent** (4,036 at every `pop_frac_low` from 0.15 to
0.40); reference-independent (agreement 0.966–0.997); and **clone composition** — all four
MPNST_2 subclones and seven of eight MPNST_4 subclones are fibroblast-dominated.

**`Candidate-Malignant-Unresolved` — resolved as malignant.** 836 Malignant, 395 Ambiguous,
**0 Non-malignant**. Highest CNA burden of any population (median 0.328) and highest
CNV-elevated fraction (65.8%) — higher than `MPNST-Tumor`. But 97% of these cells are MPNST_1,
so the malignancy call is resolved while its generality is not.

**`MPNST-Tumor` — the uncomfortable result.** SCEVAN corroborates only **36.7%**; 1,405 retained
Malignant, **2,015 → Ambiguous** via rule R5. Its marker evidence is `malignant-consistent`,
so lineage and copy-number evidence genuinely disagree, and `Ambiguous` records that rather than
picking a winner (§65C).

**`Pericyte-VSMC` — largely genuine.** 348/435 Non-malignant. Malignant only in MPNST_4 (0.890),
one patient.

**`Uncertain` — mostly non-malignant.** 476/720 Non-malignant.

**`Endothelial` — flagged, not promoted.** MPNST_4's endothelium is 52.7% malignant and forms
its own subclone (`MPNST_4_clone8`, 146/152 endothelial). Rule R11 sent those 168 cells to
`Ambiguous/Low`: a majority-malignant canonical endothelial population in **one** patient is
more likely an artefact than a discovery — flagged, never silently overridden.

## 10. Refined malignant definition

`malignancy_refined` ∈ {Malignant, Non-malignant, Ambiguous, Excluded-low-quality} with
`malignancy_confidence` ∈ {High, Moderate, Low} and a printed `malignancy_reason` naming the
rule. **`Ambiguous` is a terminal outcome and never becomes tumour.** **R99 fail-safe fired for
0 cells**, so rule coverage is complete.

**Amendment A1 (post-hoc, labelled)**: a per-sample immune sanity gate; a run calling >25% of
canonical immune cells malignant cannot drive malignant **promotions**. Deliberately asymmetric
— an unreliable CNV run may still support a `Non-malignant` call, because those always require
concordant lineage evidence too.

**Amendment A2 (post-hoc, labelled, and prompted by the gate excluding the strongest
evidence)**: the flat gate failed MPNST_4, sending 2,967 cells — including its 2,887 malignant
fibroblasts — to `Ambiguous`. Investigation showed MPNST_4's failure was **98.2% of its plasma
cells**, an immunoglobulin-locus (IGH 14q32 / IGK 2p11 / IGL 22q11) artefact documented for this
method family, with 0 of 9 populations failing and the highest run agreement of the four — a
completely different phenomenon from MPNST_3's global inversion. A2 therefore excludes plasma
cells from the sanity denominator, adds a **breadth** criterion so a global inversion cannot
dilute past the gate (which independently condemns MPNST_3), and adds rule **R7p** so 1,125
plasma cells are called `Non-malignant` with the artefact named rather than dropped.

**The rules document states this ordering explicitly**, publishes both gates for all four
samples with a `flat_gate_would_have_failed` column, and records the two independent checks a
sceptical reader should apply. **The fibroblast conclusion's direction survives rejecting A2**
(still majority-malignant in MPNST_1 and MPNST_2, which pass under both gates); **its magnitude
does not**.

## 11. Refined malignant fraction

| Stage | Malignant | Fraction of 19,716 |
| --- | ---: | ---: |
| Phase 2 conservative | 3,420 | 17.35% |
| SCEVAN primary | 8,794 / 18,449 assessed | 47.7% of assessed |
| **Phase 4 refined** | **6,434** | **32.63%** |
| — High confidence | 3,261 | 16.54% |
| — Moderate confidence | 3,173 | 16.09% |

Plus 9,078 Non-malignant, **3,766 Ambiguous (19.10%)**, 438 Excluded.

**Was 17.35% too conservative? Yes — but the correction is two-directional, and Phase 2 did not
so much undercount as look in the wrong place.** Phase 4 adds 4,036 fibroblasts, 836
candidate-malignant, 92 uncertain and 65 pericytes, while withdrawing confidence from 2,015 of
the 3,420 cells Phase 2 called malignant.

**Critical caveat**: dropping MPNST_4 gives a refined fraction of **0.214** against a Phase 2
fraction of **0.215** — fold change 1.00. **The 32.63% headline is not patient-robust**, because
outside MPNST_4 the promotions and demotions roughly cancel. The *fibroblast* conclusion is a
separate claim and **is** robust.

**Threshold sensitivity, published rather than discovered later**: at `pop_frac_low` 0.15–0.20
the refined count would be 8,295 (42.07%) with 3,266 `MPNST-Tumor` retained; at 0.25–0.40 it is
6,434 (32.63%) with 1,405 retained. **Fibroblast malignant is 4,036 at every grid point.** The
`MPNST-Tumor` retention is the most fragile number in Phase 4, because MPNST_1's pop_frac is
**0.222**, just under the a priori 0.25.

## 12. SCEVAN clone structure

22 subclones, **patient-scoped by construction**: MPNST_1 7, MPNST_2 4, MPNST_3 3 (immune
artefacts), MPNST_4 8. No patient is dominated by one subclone — the largest holds 24–34% of
the malignant compartment in the reliable patients, so subclonal structure is genuinely
polyclonal.

**Clone composition is the most direct evidence in Phase 4**, because SCEVAN builds clones from
CNA profiles without reference to any Phase 2 label. Every MPNST_1 clone mixes
`Candidate-Malignant-Unresolved`, `Fibroblast` and `MPNST-Tumor` together: marker-based
annotation split into three identities what copy-number evidence sees as **one malignant
compartment with seven subclones**.

## 13. Tumour transcriptional states

8 states from a malignant-only reduction (`malignant_pca` → `malignant_harmony` →
`malignant_umap`, dims 1:30, 93.0% of variance). Phase 2's reductions were never modified —
asserted in §23. The new reduction was justified by measurement: kNN(k=20) overlap with the
Phase 2 Harmony embedding is only **0.1733**.

| state | cells | frac | patients | dominant share | clones | defining programme |
| --- | ---: | ---: | :--: | ---: | ---: | --- |
| Mesenchymal_ECM-1 | 1,851 | 0.287 | 1 | 1.00 | 7 | ECM 0.619 |
| Mesenchymal_ECM-2 | 1,598 | 0.248 | 4 | 0.99 | 13 | ECM 0.956 |
| Mesenchymal_ECM-3 | 1,291 | 0.200 | 1 | 1.00 | 8 | ECM 1.095 |
| Mesenchymal_ECM-4 | 648 | 0.100 | 3 | 1.00 | 5 | ECM 1.255 |
| Mesenchymal_ECM-5 | 505 | 0.078 | 1 | 1.00 | 7 | ECM 0.887 |
| Cycling | 290 | 0.045 | 4 | 0.83 | 19 | Cycling 0.960 |
| **Schwann_like** | 179 | 0.027 | 2 | 0.67 | 3 | Schwann 1.222 |
| Interferon | 72 | 0.011 | 1 | 1.00 | 2 | Interferon 1.140 |

**Headline negative result (§36): 0 of 8 states are recurrent** (≥3 patients each contributing
≥5% and ≥10 cells); 1 is shared between 2 patients; **97.2% of malignant cells sit in
patient-private states**. Expected with n = 4 and `sample_id` = patient = dataset, but reported
as a negative result about tumour-state generality.

**The coherent finding underneath**: Phase 2's `MPNST-Tumor` cells map to Cycling (248),
Schwann_like (175) and Interferon (72) — the marker-legible states — while the promoted
fibroblasts map to the five Mesenchymal_ECM states (1,245 + 1,128 + 637 + 512 + 481).
**Phase 2 found the Schwann-like and cycling malignant cells and missed the Mes-NC-like ECM
ones, because an ECM programme is what a fibroblast looks like.** That is why copy-number
evidence, not a better marker panel, was the right instrument.

## 14. Clone versus tumour-state relationship

Kept as distinct concepts (§33). `tumor_clone_phase4` is CNV architecture;
`tumor_state_phase4` is transcriptional phenotype. Clone counts per state range from **2 to
19**, so a single state routinely spans many clones — one-to-one correspondence is absent.
`CLONE_VS_TUMOR_STATE.tsv` cross-tabulates them and the descriptive Cramér's V is reported as an
**association, not an equivalence**.

## 15. Phase 3 CCC sensitivity analysis

Phase 3's LIANA, CellChat, CellPhoneDB and concordance scripts were reused **completely
unmodified** at identical versions, resources, thresholds, seed and expression basis. Only
`ccc_label` differs — which is what makes this a label-sensitivity analysis rather than a
tool-drift confound. Refined input: `ccc_input_object_refined.rds`, **15,036 cells, 14
populations**, md5 `99fd641faf0cc1f7e38448406c883f32`.

| | Phase 3 | Phase 4 | change |
| --- | ---: | ---: | ---: |
| distinct keys | 502,846 | 494,973 | −1.6% |
| **supported** | **36,486** | **34,893** | **−4.4%** |
| **High concordance** | **5,848** | **5,698** | **−2.6%** |
| Moderate | 499 | 513 | +2.8% |
| Single-method | 24,654 | 23,302 | −5.5% |
| Discordant | 5,485 | 5,380 | −1.9% |

**4,036 cells entered the tumour compartment and 2,015 left, and high-concordance interactions
fell by only 2.6%.** Change classes over 506,605 union keys: Stable 18,276 · Weakened 8,446 ·
Newly-supported 5,421 · Lost 3,963 · **Sender-reassigned 3,051** · Strengthened 2,746.

**Median change across the 23 named axes: −3.1%. No axis was lost.**

* **effectively unchanged (≤5%)**: APP→CD74 −3.1% · ANXA1→FPR1 0.0% · HLA-F→LILRB1 −1.4% ·
  HLA-F→LILRB2 −3.1% · HLA-E→KLRC1 −2.0% · CD58→CD2 −1.4% · JAG1→NOTCH3 0.0% · JAG2→NOTCH2 0.0%
  · FGF2→FGFR1 −4.3%
* **strengthened**: FN1→ITGAV_ITGB8 **+33.3%** · CD99→PILRA +14.0% · VEGFA→KDR +11.1% ·
  JAG1→NOTCH2 +7.5%
* **genuinely weakened (>20%)**: COL1A1→ITGAV_ITGB8 **−44.4%** · SLIT2→ROBO1 −38.9% ·
  DLL4→NOTCH2 −31.6% · COL6A2→ITGAV_ITGB8 −30.0% · VEGFA→FLT1 −28.0%
* **sender reassigned**: all four ECM→integrin axes, plus VEGFA→KDR, JAG1→NOTCH3, DLL4→NOTCH2

**The fibroblast question (§43)**: of 9,483 Phase 3 supported interactions involving
`Fibroblast`, **6,547 (69.0%) retained**, **1,565 (16.5%) reassigned to the refined tumour
compartment**, 1,371 (14.5%) lost or no longer testable. All four collagen/FN1 → ITGAV_ITGB8
axes changed sender: the ECM→integrin interaction Phase 3 read as stroma-to-tumour is
substantially **tumour-autocrine**.

## 16. Tumour-state-specific communication

**Cannot be established in this cohort, and the §68 model structure was NOT imposed.** Of 176
evidence rows, **12 reach ≥3 patients, 34 reach ≥2, 121 rest on one patient** — and **all 12 of
the ≥3-patient rows belong to `Cycling`**, a 290-cell state that is itself 83% one patient.
`Cycling` leads all six programmes purely because it is the only state present in enough
patients to be evaluated: an **evaluability artefact, not a biological preference**. The five
`Mesenchymal_ECM` states hold **91.3%** of the malignant compartment and have **exactly one
evaluable patient each**.

Deliberately targeted (§45): LR frameworks were **not** run per state. Evidence is per-state,
per-patient expression of the prioritized ligands (40/40 present) and receptors (39 present) at
the Phase 3 detection floor, joint only when both sides pass **in the same patient**.
**No composite score was invented** — every column is an independent stream.

## 17. Receiver-response integration

Phase 3's NicheNet and LochNESS results are **reused verbatim, never recomputed**, and NicheNet
is **not** forced to support APP–CD74. Phase 3 established that predicted APP–CD74 engagement
and the CSF1-explained macrophage state are **distinct findings**, and Phase 4 preserves that:
the myeloid ligand programme (APP, CD99, ANXA1, HLA-F, THBS1) and the receiver-state cytokine
programme (CSF1, IL15, TGFB1) are carried as **separate** streams and never merged. LochNESS
enters only as receiver-lineage context and remains a Phase 3 negative (best descriptive
p = 0.333 against a 0.167 floor at n = 4, with the score itself unstable).

## 18. Robustness

| claim | robust? | basis |
| --- | :--: | --- |
| Phase 2's `Fibroblast` compartment is predominantly malignant | **yes** | 3/4 patients; threshold-independent 0.15–0.40; reference-independent 0.966–0.997; clone composition |
| The Phase 3 CCC architecture survives refinement | **yes** | High concordance −2.6%; median axis change −3.1%; no axis lost; LOSO 0.727–0.955 |
| ECM→integrin interactions are substantially tumour-autocrine | **yes** | all four axes sender-reassigned; 1,565 reassignments |
| A small shared CNA core (chr18 loss, chr2/chr7 gain) | qualified | 3/3 reliable patients, but n = 3 and Jaccard only 0.20–0.23 |
| `Candidate-Malignant-Unresolved` is malignant | qualified | 0 non-malignant, highest CNA burden — but 97% one patient |
| **Refined malignant fraction is 32.63%** | **no** | equals the Phase 2 fraction once MPNST_4 is dropped |
| Malignant transcriptional states are recurrent | **no** | 0 of 8; 97.2% patient-private |
| A tumour-state → TME model | **no** | 12 of 176 rows at ≥3 patients, all in one 290-cell state |
| Newly-supported interactions after refinement | **no** | **91.8% single-patient**; only 27 of 5,421 at ≥3 patients |

Leave-one-patient-out retention of the 34,893 refined supported interactions: MPNST_1 **0.727**,
MPNST_2 0.825, MPNST_3 0.916, MPNST_4 0.955 — no patient dominates.

## 19. Biological model

Stated only as far as the evidence reaches.

**The malignant compartment of these MPNSTs is roughly twice as large as marker-based annotation
indicated, and it is dominated by a mesenchymal/ECM-programme population that marker-based
annotation read as fibroblasts.** 4,036 of 5,064 `Fibroblast` cells carry inferred CNAs,
majority-malignant in three of four patients, and SCEVAN's CNA-defined subclones are built from
those cells. The five mesenchymal/ECM states hold 91.3% of the refined malignant compartment,
while the marker-legible Schwann-like, cycling and interferon states account for under 9%.

**The Phase 3 tumour–TME architecture is not an artefact of the Phase 2 tumour definition.**
Myeloid-directed signalling (APP→CD74, CD99→PILRA, ANXA1→FPR1, HLA-F→LILRB1/2), the NK
inhibitory axis (HLA-E→KLRC1), the angiogenic axis (VEGFA→KDR/FLT1/NRP1) and the reciprocal
perivascular Notch circuit all survive refinement essentially unchanged, with a median axis
change of −3.1% and no axis lost.

**What refinement changed is attribution, not existence.** 3,051 interactions were
sender-reassigned, concentrated in the ECM–integrin axis: the collagen and fibronectin signalling
Phase 3 attributed to fibroblasts is substantially **tumour-autocrine**, from malignant
Mes-NC-like cells to malignant cells.

**Copy number is patient-individual with a small shared core** — chr18 loss and chr2/chr7 gain
in all three reliable patients, pairwise Jaccard only 0.20–0.23 — and one patient carries a
38.1 Mb single-copy chr17 segment containing the NF1 locus.

**What cannot be claimed**: that any tumour state recurs across patients (0 of 8), that a
specific state drives a specific TME programme (unassessable — 97.2% of malignant cells are in
patient-private states), or that the 32.63% malignant fraction is a cohort property (it is not
patient-robust).

## 20. Limitations

**A.** SCEVAN infers copy number from expression. **Not DNA sequencing.**
**B.** Reliable for broad chromosomal, arm-level and large-segment events; **single-gene CNV
calls are not asserted**.
**C.** Some MPNST malignant cells may be copy-number quiet, so **a SCEVAN normal call is not
proof of non-malignancy** — this is why 2,015 `MPNST-Tumor` cells became `Ambiguous` rather than
`Non-malignant`.
**D.** Tumour state and CNV clone are not equivalent and are kept separate.
**E.** **n = 4.** No population-level or epidemiological claim follows.
**F.** `sample_id` = patient = dataset, so biological and technical effects cannot be fully
separated.
**G.** SCEVAN's annotation covers **chr1–22 only**; X and Y events are not assessable.
**H.** SCEVAN removes cell-cycle genes **and all `HLA-*` genes** before inference, so the CNV
analysis is **structurally blind to the HLA-E / HLA-F loci** Phase 3 highlighted.
**I.** Per-sample reliability differs. MPNST_3 failed the sanity gate and its promotions were
disabled, costing sensitivity in that patient.
**J.** `Ambiguous` (3,766 cells) is a terminal outcome, excluded from the refined CCC analysis
as **NOT EVALUABLE** — which is not the same as "no signalling".
**K.** The clone phylogeny plot is absent (ggtree/ggplot2 4.x); ggplot2 was not downgraded and
no clone data is affected.
**L.** Amendment A2 was written after observing that the flat gate excluded the sample carrying
the strongest evidence. Both gates are published; the fibroblast conclusion's direction survives
rejecting A2, its magnitude does not.
**M.** The 32.63% refined fraction is **not patient-robust** (§11).
**N.** `Fibroblast` retains 1 cell in MPNST_4 and `Pericyte-VSMC` 4, so those populations are
NOT EVALUABLE there, which is why the change classes separate `Ambiguous` from `Lost`.

## 21. Final figures

`results/phase4/figures/final/` — `01_scevan_malignancy_umap` · `02_phase2_vs_scevan_malignancy`
· `03_refined_malignancy_umap` · `04_scevan_cna_heatmap` · `05_scevan_clones_by_patient` ·
`06_ambiguous_population_malignancy` · `07_refined_malignant_fraction` · `08_tumor_state_umap` ·
`09_tumor_state_marker_heatmap` · `10_clone_vs_tumor_state` · `11_phase3_vs_phase4_ccc` ·
`12_refined_tumor_to_immune` · `13_refined_tumor_to_vascular` · `14_tumor_state_specific_ccc` ·
`15_final_tumor_state_tme_model` · `16_cnv_profiles_by_annotation`, plus every supporting
figure. PDF primary, PNG alongside.

## 22. Final tables

`results/phase4/tables/final/` — `PHASE4_MALIGNANCY_CALLS.tsv` (19,716 × 37, per-cell rule and
printed reason) · `SCEVAN_CELL_CLASSIFICATION.tsv` · `SCEVAN_CLONES.tsv` ·
`SCEVAN_CNV_SUMMARY.tsv` · `SCEVAN_VS_PHASE2_ANNOTATION.tsv` ·
`AMBIGUOUS_POPULATION_MALIGNANCY_AUDIT.tsv` · `SCEVAN_SAMPLE_RELIABILITY.tsv` ·
`MALIGNANCY_THRESHOLD_SENSITIVITY.tsv` · `MALIGNANCY_SUMMARY_BY_{CLUSTER,SAMPLE}.tsv` ·
`TUMOR_STATE_{ASSIGNMENTS,MARKERS,PATIENT_DISTRIBUTION}.tsv` · `CLONE_VS_TUMOR_STATE.tsv` ·
`PHASE3_VS_PHASE4_CCC.tsv` · `PHASE3_AXIS_SENSITIVITY.tsv` · `REFINED_TUMOR_CCC.tsv` ·
`TUMOR_STATE_TME_EVIDENCE.tsv` · the M34 robustness tables.

## 23. Final object

```text
path      results/phase4/phase4_final_object.rds
size      5.64 GB
md5       e85ba8486e456917e2483f2773bdbaf3
sha256    a658447744e5ee8621fc5a3127492ceb11bce37bac4d6f3f0f7255bdcb647ae4
cells     19,716          features  31,764 (RNA) / 29,113 (SCT)
assays    RNA, SCT        reductions  pca, umap_preintegration, postint_harmony, postint_umap_harmony
metadata  149 Phase 1-3 columns + 32 Phase 4 columns
```

**Validation: 23/23 checks passed** after saving and reloading under SLURM. Preservation guards,
all passed:

```text
all 149 pre-existing metadata columns are byte-identical
reductions unchanged; EVERY Phase 2 embedding numerically identical (PCA and Harmony included)
assays and cell order unchanged
annotation_ccc_phase3 is a verbatim copy of the frozen Phase 3 layer
no non-Malignant cell carries the refined MPNST-Tumor label
tumour states exist only on Malignant cells
save/load md5 stable
0 non-finite values across all four embeddings
```

**Phase 4 fields added** (nothing overwritten): `scevan_call` · `scevan_confident_normal` ·
`scevan_clone` · `scevan_subclone_raw` · `scevan_sample` · `scevan_sample_reliable` ·
`population_class` · `cnv_burden` · `cnv_frac_gain` · `cnv_frac_loss` · `cnv_mean_abs` ·
`cnv_ref_threshold` · `cnv_elevated` · `cnv_summary` · `malignancy_phase2` ·
`malignancy_scevan` · `malignancy_scevan_sens` · `malignancy_refined` ·
`malignancy_confidence` · `malignancy_rule` · `malignancy_reason` · `annotation_ccc_phase3` ·
`annotation_ccc_refined` · `marker_evidence` · `score_malignant_schwann_nc` ·
`score_nonmalignant_max` · `pop_frac` · `pop_n_assessed` · `runs_agree` ·
`tumor_state_phase4` · `tumor_clone_phase4` · `tumor_cluster_phase4`.

## 24. Reproduction

```bash
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
sbatch scripts/shell/phase4/install_umap_learn.sh          # one-off, isolated env
J28=$(sbatch --parsable scripts/shell/phase4/run_m28_feasibility.sh)
J29=$(sbatch --parsable --dependency=afterok:$J28 scripts/shell/phase4/run_m29_scevan.sh)
J30=$(sbatch --parsable --dependency=afterok:$J29 scripts/shell/phase4/run_m30.sh)
bash scripts/shell/phase4/submit_m32_to_m35.sh             # after M30 completes
```

Seeds fixed at 42 throughout. Every script exports `RETICULATE_PYTHON` to the isolated
interpreter and matches `par_cores` to `--cpus-per-task`.

## 25. Environment

**The frozen Phase 2/3 R stack did not move.** Verified before and after every operation:

```text
R 4.4.3 · Seurat 5.4.0 · SeuratObject 5.3.0 · Matrix 1.7.4 · harmony 1.2.4 · sctransform 0.4.3
SCEVAN 1.0.3 · yaGST 2017.8.25          (both ALREADY INSTALLED - no installation attempted)
liana 0.1.14 · CellChat 2.2.0.9001 · nichenetr 2.2.1.1 · CellPhoneDB 5.0.1 (isolated cpdb_env)
```

**Dependency changes to `R_env` in Phase 4: none.** One addition anywhere: Python
**umap-learn 0.5.12** in an isolated `p4_umap_env`, required because SCEVAN's subclone stage
calls `umap::umap(method = "umap-learn")` before writing clone labels. Snapshots in
`reports/phase4/environment/`; full account in `PHASE4_INSTALL_LOG.md`.

## 26. SLURM accounting

| JobID | Milestone / stage | State | Elapsed | ReqMem | MaxRSS |
| --- | --- | --- | --- | --- | --- |
| 19896576 | M28 attempt 1 | **FAILED** | 00:02:15 | 96 G | 27.60 GiB |
| 19896623 | M28 attempt 2 | **FAILED** | 00:03:14 | 96 G | 27.60 GiB |
| 19896654 | umap env attempt 1 | **FAILED** | 00:00:11 | 16 G | 0 |
| 19896672 | umap env | COMPLETED | 00:05:39 | 16 G | 0.48 GiB |
| 19896711 | **M28** | COMPLETED | 00:05:07 | 96 G | 26.90 GiB |
| 19896712_1 | M29 MPNST_1 | COMPLETED | 00:34:20 | 64 G | **60.87 GiB** |
| 19896712_2 | M29 MPNST_2 | COMPLETED | 00:08:13 | 64 G | 6.24 GiB |
| 19896712_3 | M29 MPNST_3 | COMPLETED | 00:11:20 | 64 G | 18.10 GiB |
| 19896712_4 | M29 MPNST_4 | COMPLETED | 00:42:07 | 64 G | 50.49 GiB |
| 19897096 | M30 (pre-A2) | COMPLETED | 00:01:06 | 64 G | 0.57 GiB |
| 19897148 | **M30** | COMPLETED | 00:01:00 | 64 G | 1.38 GiB |
| 19897175 | M32 refined input | COMPLETED | 00:07:04 | 96 G | 27.78 GiB |
| 19897176 | M32 LIANA | COMPLETED | 00:34:36 | 64 G | 38.87 GiB |
| 19897177_1..4 | M32 CellChat | COMPLETED | 00:19:40–00:39:48 | 48 G | 14.49–18.60 GiB |
| 19897178 | M32 CellPhoneDB | COMPLETED | 00:02:31 | 48 G | 1.64 GiB |
| 19897179 | M32 concordance+sens | **FAILED** | 00:03:39 | 48 G | 1.34 GiB |
| 19897616 | M32 sensitivity rerun | COMPLETED | ~00:01 | 48 G | — |
| 19897617 | M31/33/34 (pre-fix) | COMPLETED | 00:02:24 | 64 G | 5.16 GiB |
| 19897639 | **M31/M33/M34** | COMPLETED | 00:02:22 | 64 G | 5.41 GiB |
| 19897640 | **M35** | COMPLETED | — | 96 G | — |

Cancelled without scientific output: 19897097–19897111 (two duplicate dependency chains fired
simultaneously by a waiter and a monitor) and 19897180/19897181 (dependency never satisfied
after 19897179 failed).

**Peak memory 60.87 GiB — 13.5% of the 450 G envelope. No failure was ever addressed by
increasing RAM or walltime, and no Phase 2/3 package was downgraded to make a tool work.**
MPNST_1's M29 task peaked at 95% of its 64 G request; that is recorded because the estimate was
closer than intended and a larger sample would need more.

## 27. Recommended next phase

**Phase 5 is NOT initiated and requires separate authorization.**

Two things Phase 4 makes worth doing, and one it rules out:

1. **Targeted spatial validation of the tumour-autocrine ECM axis.** The firmest and most
   actionable Phase 4 result is that collagen/FN1–ITGAV_ITGB8 signalling Phase 3 read as
   stroma-to-tumour is substantially tumour-to-tumour. That is directly testable *in situ*:
   dual staining for a malignant marker and COL1A1/FN1 would distinguish malignant Mes-NC-like
   cells from genuine CAFs in tissue, which no amount of further scRNA-seq analysis can do.
2. **Experimental confirmation of malignant identity in the fibroblast compartment.** The
   claim rests on inferred copy number. NF1/NF2 status by DNA sequencing, or FISH on the
   recurrent chr18 loss / chr2 / chr7 gains, would convert an inference into a determination.
3. **A tumour-state → TME model is NOT worth attempting on this cohort.** With 0 of 8 states
   recurrent and 97.2% of malignant cells in patient-private states, more analysis of these four
   patients cannot separate state from patient. That question needs more patients, not more
   method.

---

## 28. M35A — SCEVAN Figure Consolidation

*Appended 2026-09-03. **Nothing above this line was rewritten.** M35A is a visualization
milestone: it surfaces and consolidates evidence that already existed. SCEVAN was not re-run, and
no malignancy call, decision rule, amendment A1/A2, clone assignment, tumour state, CCC result or
threshold was changed. `results/phase4/phase4_final_object.rds` was neither loaded nor modified;
its md5 `e85ba8486e456917e2483f2773bdbaf3` is unchanged.*

**Why this milestone was needed.** The SCEVAN CNA evidence existed but was not visible: 267 native
SCEVAN files sat behind `SCEVAN_NATIVE_FIGURE_INDEX.tsv` and no final figure showed a native CNA
heatmap. The claim that carries Phase 4 — that a large fibroblast-like compartment belongs to
CNA-defined malignant populations — was defensible in the tables and only partly visible in the
figures. Full audit: `reports/phase4/M35A_SCEVAN_FIGURE_AUDIT.md`.

### 28.1 Native figures surfaced

Per patient and run, SCEVAN produced a CNA heatmap over all cells with a malignant/normal track
(`heatmap.png`), a tumour-cells-only heatmap (`onlytumorheatmap.png`), a subclone-annotated
heatmap (`heatmap_subclones.png`), a consensus clonal CN profile (`consensus.png`), cytoband
onco-heatmaps (`OncoHeat.png`, `OncoHeat2.png`), a CNA-space UMAP (`umap_CNA.png`) and an
expression-space UMAP (`umap_scRNA.png`), plus 151 per-segment subclone DE panels and 44 pathway
panels. Machine-readable CN profiles exist as `_Clonal_CN.seg` and `_subclone<N>_CN.seg`.

**Surfaced into `results/phase4/figures/final/`:** the subclone-annotated heatmap, the all-cell
heatmap, the tumour-only heatmap and the consensus profile for the three reliable patients, as
`17_scevan_native_cna_MPNST1.pdf`, `18_..._MPNST2.pdf`, `19_..._MPNST4.pdf` (2 pages each) and
`20_scevan_cna_reliable_patients.pdf`.

**The native heatmaps were sufficient as CNA evidence and insufficient as an argument.** They
carry SCEVAN's own subclone track but no Phase 2 identity, no clone size and no refined
malignancy call, and redrawing them would have meant reconstructing SCEVAN's internal cell
ordering — fragile, and it would put the CNA matrix at risk for no scientific gain. So the native
plot is **embedded verbatim as a raster** and a **companion panel built from the same frozen clone
assignments** is placed beside it. **No CNA value was redrawn anywhere in M35A.**

**Not surfaced, with reasons.** The cytoband onco-heatmaps are unreadable at page scale and invite
gene-level over-reading, so broad-event evidence is shown at segment resolution in
`27_broad_cna_recurrence.pdf` instead. The CNA-space and expression-space UMAPs are redundant with
existing figures `05`, `31_03` and `01/02`. **All 8 `CloneTree.png` files are blank** (pixel
standard deviation exactly 0) — the documented ggtree/ggplot2 4.x failure of §4 and limitation K.
ggplot2 was not downgraded and clone phylogeny plotting was not forced.

### 28.2 Custom figures created

| figure | what it establishes |
| --- | --- |
| `21_scevan_clone_composition_phase2.pdf` | clone × Phase 2 identity, per patient, with clone size, refined-malignant fraction, fibroblast dominance and disputed-triple marks |
| `22_fibroblast_cna_burden_by_patient.pdf` | the four cell-level CNA metrics per patient, with per-patient Cliff's delta |
| `23_fibroblast_malignant_vs_nonmalignant_cna_profile.pdf` | genome-wide mean CNA profile per group, from the stored native CNA matrices |
| `24_phase2_to_phase4_malignancy_transition.pdf` | the identity transition, disputed compartments only |
| `25_MPNST3_scevan_failure_qc.pdf` | MPNST_3 shown failing |
| `26_fibroblast_malignancy_threshold_robustness.pdf` | robust fibroblast count vs fragile cohort fraction |
| `27_broad_cna_recurrence.pdf` | broad recurrent architecture, segment resolution, no gene-level claim |
| `28_phase4_scevan_evidence_summary.pdf` | the whole argument on one page, panels A–H |

Supporting tables: `results/phase4/tables/final/SCEVAN_CLONE_COMPOSITION_VISUALIZATION.tsv` (176
rows), `M35A_FIBROBLAST_CNA_EFFECT_SIZES.tsv`, `M35A_FIBROBLAST_CNA_PROFILE_CORRELATION.tsv`.

### 28.3 Why clone composition is the strongest malignancy evidence

SCEVAN builds its subclones from inferred CNA profiles **without reference to any Phase 2 label**,
so asking what a CNA-defined clone is made of is not circular — which is exactly what the
malignant-fraction-per-population evidence cannot claim, because that is a per-population summary
of the same call.

* **MPNST_2 — 4 of 4 clones fibroblast-dominated** (229/234, 160/186, 138/158, 110/116 Fibroblast).
* **MPNST_4 — 7 of 8 clones fibroblast-dominated**; the eighth (`MPNST_4_clone8`, 152 cells) is
  146/152 endothelial and was sent to `Ambiguous/Low` by rule R11 rather than promoted.
* **MPNST_1 — all 7 clones mix `Candidate-Malignant-Unresolved` with `MPNST-Tumor`, and 6 of the 7
  also carry `Fibroblast` cells.** Marker-based annotation split into three identities what copy
  number sees as one malignant compartment with seven subclones. See §28.7 for the correction this
  forced to the wording of §12.

### 28.4 How malignant and non-malignant fibroblast-labelled cells compare

Both comparisons point the same way, and both are reported per patient because cells within a
patient are not independent — a pooled cell-level p-value over 5,064 cells would manufacture
significance from n = 4 patients.

**CNA burden** (Cliff's delta, Fibroblast Malignant vs Non-malignant; +1 = complete separation):

| patient | n Mal / n Non | cnv_burden | cnv_frac_gain | cnv_frac_loss | cnv_mean_abs |
| --- | ---: | ---: | ---: | ---: | ---: |
| MPNST_1 | 512 / 303 | **0.894** | 0.915 | 0.834 | 0.917 |
| MPNST_2 | 637 / 503 | **0.836** | 0.781 | 0.811 | 0.890 |
| MPNST_3 | 0 / 101 | NOT EVALUABLE (promotions disabled) | | | |
| MPNST_4 | 2,887 / 1 | NOT EVALUABLE (1 non-malignant fibroblast) | | | |

Median `cnv_burden` 0.369 vs 0.184 (MPNST_1) and 0.311 vs 0.186 (MPNST_2).

**CNA profile shape** (Pearson correlation of the group mean profiles, from the stored native CNA
matrices):

| patient | r(Fib→Mal, MPNST-Tumor→Mal) | r(Fib→Mal, Fib→Non-mal) | r(Fib→Mal, immune control) |
| --- | ---: | ---: | ---: |
| MPNST_1 | **0.934** | 0.184 | −0.055 |
| MPNST_2 | not evaluable (7 MPNST-Tumor→Malignant cells) | 0.346 | −0.387 |
| MPNST_4 | **0.972** | not evaluable (1 cell) | −0.381 |

**Answer, stated at the strength the data supports.** Where the contrast can be made, malignant-called
fibroblast-labelled cells resemble the malignant compartment far more than they resemble
non-malignant fibroblasts — 0.934 against 0.184 in MPNST_1, 0.972 against an immune baseline of
−0.381 in MPNST_4. **The contrast exists in 2 of 4 patients only**: MPNST_4 retains a single
non-malignant fibroblast and MPNST_3 contributes no malignant fibroblasts. The separation is not
weak; its breadth is limited, and both facts are printed on the figures.

### 28.5 How the MPNST_3 failure is demonstrated

`25_MPNST3_scevan_failure_qc.pdf` shows four things: primary↔sensitivity agreement **0.0864**
against 0.9663–0.9974 elsewhere (read from the frozen `m29_*_summary.json`, which counts filtered
cells as agreeing — 254 of 2,940); MPNST_3's three "clones" resolving to T/NK (clone 3, 674 cells),
plasma/pDC/B (clone 2, 402) and plasma-dominated (clone 1, 338); six of nine canonical immune
populations called ~100% malignant against the a priori 25% gate; and sequencing depth as
**context only** — median 1,594 genes/cell against 3,139 / 2,490 / 2,114 — with the figure stating
explicitly that low depth alone is not claimed to have caused the failure.

> MPNST_3 malignant promotions were disabled because SCEVAN's inferred tumour partitions
> corresponded predominantly to canonical immune lineages and the primary and sensitivity
> classifications were incompatible.

### 28.6 How threshold sensitivity affects the interpretation

`26_fibroblast_malignancy_threshold_robustness.pdf` separates two claims that are easy to conflate.
**Fibroblast → Malignant is 4,036 at every tested `pop_frac_low`** (0.15, 0.20, 0.25, 0.30, 0.40),
as is Candidate-Malignant-Unresolved at 836. **Retained `MPNST-Tumor` moves 3,266 → 1,405** and the
refined fraction moves **42.07% → 32.63%**, because MPNST_1's `pop_frac` is 0.222, just under the a
priori 0.25. The 32.63% headline is additionally not patient-robust (§11). The fibroblast
conclusion and the cohort fraction are separate claims and only the first is robust.

### 28.7 One discrepancy found against this handoff, recorded not smoothed

All frozen counts reproduced exactly: 19,716 cells; Fibroblast 5,064 → 4,036 / 908 / 120;
Candidate-Malignant-Unresolved 836 / 395 / 0; MPNST-Tumor 1,405 / 2,015; refined 6,434 / 9,078 /
3,766 / 438; clones 7 / 4 / 3 / 8; agreement 0.9951 / 0.9663 / 0.0864 / 0.9974.

**§12 above states that "Every MPNST_1 clone mixes `Candidate-Malignant-Unresolved`, `Fibroblast`
and `MPNST-Tumor` together". The frozen clone assignments show this holds for 6 of the 7 clones.**
`MPNST_1_clone6` (176 cells) is `Candidate-Malignant-Unresolved` = 108 and `MPNST-Tumor` = 68 with
**no `Fibroblast` cells**; `MPNST_1_clone7` (54 cells) contains 1. The accurate statement is that
**all 7 clones mix `Candidate-Malignant-Unresolved` with `MPNST-Tumor`, and 6 of 7 additionally
carry `Fibroblast` cells**. No call, count, clone assignment or conclusion changes — the two
exceptions are the smallest mixed clones, 230 of MPNST_1's 1,952 malignant cells, and the
underlying claim stands on all seven. §12 is left as written; this is the correction of record and
figure 21 plots the composition as it actually is.

### 28.8 SLURM accounting

| JobID | stage | State | Elapsed | ReqMem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: |
| 19899301 | M35A figures, first pass | COMPLETED | 00:04:29 | 24 G | 2.11 GiB |
| 19899307 | M35A figures, label/wrapping fixes | COMPLETED | 00:04:25 | 24 G | 1.95 GiB |
| 19899308 | M35A figures, summary-layout fix | COMPLETED | 00:04:27 | 24 G | 1.96 GiB |
| 19899311 | M35A figures, genome-axis and panel-order fix | COMPLETED | 00:04:25 | 24 G | 2.01 GiB |
| 19899312 | **M35A figures, final** | COMPLETED | 00:04:24 | 24 G | 2.97 GiB |
| 19899313 | M35A finalize, first attempt | COMPLETED | 00:00:10 | 8 G | — |
| 19899315 | **M35A finalize** (checksums, figure index, manifest) | COMPLETED | 00:00:09 | 8 G | — |

Every run exited 0; the reruns are figure-quality iterations found by inspecting the rendered
output, not failures. **Peak memory 2.97 GiB against a 24 G request — 0.7% of the 450 G envelope.**
No plotting problem was addressed by increasing memory, and the 6 GB final object was never loaded:
the plotting input is the 5 MB `results/phase4/malignancy/phase4_malignancy_metadata.rds`, the
small frozen TSVs, and the stored native SCEVAN PNG/RData outputs.

### 28.9 Environment

`R_env` unchanged. ggplot2 4.0.3, ggtree and every other package left exactly as frozen; nothing
installed, upgraded or downgraded. No CNV method was added.

### 28.10 A defect caught in the manifest writer, recorded because it nearly did damage

Job **19899313** completed successfully and produced a manifest that was **wrong**. `jsonlite::toJSON`
defaults to `digits = 4`, so re-serialising the manifest silently rounded frozen values already in
it — run agreement `0.99514117 → 0.9951`, refined fraction `0.32633394 → 0.3263`, every
`elapsed_min`, every `pop_frac` and every CNA reference threshold. Nothing errored; the job exited
0.

It was caught by diffing the manifest against a pre-write backup rather than by trusting the exit
code. `m35a_finalize.R` now writes with `digits = NA` and, after writing, **re-reads the manifest and
asserts every pre-existing section is identical to a snapshot taken before the write** — 33 sections
verified unchanged, with only `figures`, `tables` (appended to) and the new
`m35a_figure_consolidation` block differing. The manifest was restored from backup and regenerated
under job **19899315**.
