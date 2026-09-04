# Phase 3 Handoff — MPNST Tumour–Immune Cell–Cell Communication

*Frozen 2026-09-03 · Phase 3 status: **COMPLETE***

Self-contained: a new researcher or agent can reconstruct Phase 3 from this document without
any conversation history.

> **Interpretation rule that governs every statement below.** scRNA-seq ligand–receptor
> analysis infers communication **potential** from expression and associated transcriptional
> programmes. It does **not** establish physical cell adjacency or direct signalling.
> Conclusions use *predicted, inferred, supported, consistent with, concordant, associated
> with* — never *proves* or *demonstrates*.

---

## 1. Objective

Determine how MPNST tumour cells interact with immune and other tumour-microenvironment
populations, identify recurrent and biologically plausible signalling circuits, establish
which findings are concordant across independent CCC frameworks, and connect candidate
tumour-derived signals to receiver-cell transcriptional states.

**Not a methods-benchmarking study.** Multiple tools were used because CCC inference is
method-dependent; agreement across independent frameworks raises confidence, and disagreement
is reported rather than averaged away.

## 2. Phase 2 input

| Item | Value |
| --- | --- |
| Object | `results/phase2/phase2_final_object.rds`, md5 `153d5f6acc70f9c05aa48cabc4f4ac2d` |
| Cells | 19,716 · RNA 31,764 features · SCT 29,113 (4 SCT models) |
| CCC field | `annotation_ccc` (17 identities) |
| Sample field | `sample_id` — simultaneously **dataset and patient** |
| Condition field | **does not exist** |
| Tumour | `MPNST-Tumor`, 3,420 cells (17.35%) |

## 3. CCC populations

Three Phase 2 not-CCC-ready populations were excluded outright
(`Candidate-Malignant-Unresolved` 1,231, `Uncertain` 720, `Low-quality-excluded` 438).
**14 populations · 17,327 cells** remained:

`MPNST-Tumor` | `Macrophage` `Monocyte` `Dendritic` `Plasmacytoid-DC` `CD8-T` `CD4-T` `NK`
`T-cell-other` `B-cell` `Plasma-cell` | `Fibroblast` `Pericyte-VSMC` | `Endothelial`

No Mast population exists at Harmony resolution 1.0, so none was analysed and none fabricated.

## 4. Sample / patient coverage

**Every one of the 14 populations, and all 182 directed pairs, is evaluable in all four
samples.** Minimum-cell policy (stated, not silently chosen): ≥10 cells per population per
sample (CellChat's `min.cells` default); ≥3 populations for a sample to be analysable; and
**a population below the minimum is recorded NOT EVALUABLE, never as absence of signalling**.

## 5. CCC methods and their distinct roles

| Framework | Version | Role |
| --- | --- | --- |
| **LIANA** | 0.1.14 | primary multi-resource **consensus LR** framework |
| **CellChat** | 2.2.0.9001 | independent **pathway/network** framework |
| **CellPhoneDB** | 5.0.1 | independent **permutation-based LR** framework |
| **NicheNet** | 2.2.1.1 | **receiver-response modelling**, not a fourth LR list |
| **LochNESS** | MMCA formulation | **orthogonal neighbourhood phenotype**, not a CCC method |

**Raw scores were never averaged across tools.** Each LR framework contributes a boolean
support flag. **Declared partial non-independence:** LIANA's method set includes a
CellPhoneDB-style score.

Expression basis: **RNA assay, joined, LogNormalize** — the SCT assay's four models are not on
a common footing across samples. **Harmony coordinates were never used as expression.**

## 6. Database / resource versions

LIANA Consensus resource via **OmnipathR 3.14.0** · **CellChatDB.human** (3,233 interactions,
290 pathways) · **CellPhoneDB v5.0.0** (1,000 permutations, threshold 0.10) · **NicheNet prior
model** Zenodo record 7074291 (ligand–target 33,354 × 1,226, md5
`b09606b04b2d4490418d9028c0e58b9f`; lr_network 4,986 rows; weighted networks 73.7 MB) ·
**MMCA** `af629c498f421d6f0bfb325e9fd6b2dab0e22837` (vendored read-only).

## 7–9. Per-framework results

| Framework | Rows pooled | Supported | Tumour-involving |
| --- | ---: | ---: | ---: |
| LIANA | 192,015 | 36,684 | 6,280 |
| CellChat | 17,506 (140 pathways) | 17,506 (significant only) | — |
| CellPhoneDB | 1,460,984 | 25,394 | 5,327 |

## 10. Method concordance

**502,846 distinct keys · 36,486 supported.**

| Class | n |
| --- | ---: |
| High concordance | **5,848** |
| Moderate concordance | 499 |
| Single-method | 24,654 |
| **Discordant/ambiguous** | **5,485** |

**Resource overlap governs interpretation:** testable by 3 frameworks 1,347; by 2 13,064; by 1
**488,435**. 97% of keys exist in only one resource, so "single-method" usually means resource
non-overlap. Among the 1,347 jointly testable interactions, **763 (57%) were supported by all
three** — a high agreement rate for CCC inference.

## 11. Tumour → immune findings

**A myeloid-directed, partly immunosuppressive, antigen-presenting tumour.**

Three frameworks, 4/4 patients: **APP → CD74** (Macrophage, Dendritic, Monocyte, B-cell, pDC —
the most consistent tumour signal in the dataset) · **CD99 → PILRA** · **ANXA1 → FPR1** ·
**THBS1 → CD36** · **HLA-F → LILRB1/LILRB2**.

Lymphoid, and deliberately reported as **mixed**: **HLA-A/HLA-E → CD8A/CD8B** (class I
retained), **HLA-E → KLRC1/NKG2A** and **HLA-F → LILRB1/2** (inhibitory) alongside
**BAG6 → NCR3/NKp30** and **CD58 → CD2** (activating). This dataset cannot resolve which
dominates functionally.

## 12. Immune → tumour findings

**JAG1 / JAG2 / DLL4 → NOTCH2** from pericytes and endothelium — a perivascular Notch niche ·
**CRTAM → CADM1** (CD8-T) · **TNF → TNFRSF1A** (CD4-T) · **NCAM1 → FGFR1** and
**SEMA4D → PLXNB2** (NK) · **COL1A1/COL1A2/COL6A2/FN1 → ITGAV_ITGB8** (fibroblast, pericyte).

The Notch axis runs **both ways** (tumour JAG1 → pericyte NOTCH3; pericyte/endothelial ligands
→ tumour NOTCH2), suggesting a reciprocal perivascular circuit.

## 13. Stromal / endothelial findings

**VEGFA → KDR / FLT1 / NRP1** (4 evidence streams, 4/4 patients) · **JAG1 → NOTCH4** ·
**COL1A2 → ITGA1/2/3/9_ITGB1** · **ANGPTL4 → CDH5** (endothelial); **FGF2 → FGFR1** and
**SLIT2 → ROBO1** (fibroblast); **JAG1 → NOTCH3** (pericyte).

Notably, direct **TGFB1 → TGFBR1/TGFBR2 was only weakly supported** (1 framework, P5_low); the
better-supported route to that biology here is stromal ECM → tumour `ITGAV_ITGB8`, an integrin
that activates latent TGF-β.

## 14. Receiver-response (NicheNet)

**Myeloid receivers all converge on CSF1** (Macrophage 0.166, Monocyte 0.123, Dendritic 0.107).
**Lymphoid receivers all converge on IL15** (CD8-T 0.116, NK 0.089, CD4-T 0.084).
**Endothelial**: TGFB1, HMGB1, VEGFA, ANGPT1. **Fibroblast: a clear negative** — best AUPR
0.022, so tumour ligands do not explain the fibroblast programme.

**Reported divergence:** APP is **not** among the top NicheNet ligands for macrophages. The two
methods ask different questions, but this bounds the APP-CD74 claim to *predicted engagement*,
not *drives the macrophage state*.

## 15. LochNESS design

Official MMCA formulation, R: `k = round(0.5·√N)`, L2-normalised `postint_harmony` per lineage,
exact `FNN::get.knnx`, **disjoint query/reference sets (same-sample exclusion)**, global
fraction over the reference set. Context = sample-level split on tumour **APP** expression
(MPNST_1, MPNST_4 high; MPNST_2, MPNST_3 low) computed from **tumour cells only**, so the label
and the receiver-space score are independent. Two nulls: exact sample-level permutation
(6 labelings, **p-floor 0.167**) and 100 cell-level shuffles. Full rationale and limitations in
`LOCHNESS_MPNST_DESIGN.md`.

## 16. LochNESS results

| Lineage | mean | Null A p | Harmony-vs-PCA ρ | LOSO ρ |
| --- | ---: | ---: | ---: | ---: |
| Macrophage | −0.163 | 0.667 | 0.419 | 0.497 |
| Fibroblast | +0.104 | 0.833 | 0.023 | 0.347 |
| Endothelial | −0.176 | 0.333 | 0.219 | 0.448 |
| CD8-T | −0.311 | 0.333 | −0.265 | 0.310 |

**A negative result: no lineage shows receiver-state structure associated with the APP context
beyond chance**, and the score is unstable across embeddings and reference samples.
**This does not refute the APP-CD74 LR finding** — LochNESS is not a communication method.

**Implementation comparison (spec §36):** R (MMCA) versus Python given mathematically
equivalent inputs → **Pearson 1.0000, max abs diff 0**. Versus the perturb-seq formulation as
written → Pearson 0.3764, because that formulation omits same-sample exclusion. Both codebases
implement the formula correctly; only the surrounding design differs.

## 17. Integrated biological findings

Evidence matrix, not a composite score. Streams: ≥2 LR frameworks · ≥3/4 patients ·
expression support · NicheNet receiver response.

| Streams | Interactions |
| --- | ---: |
| 4 | **325** |
| 3 | 3,406 |
| 2 | 7,365 |

**547 tumour-centric interactions reach ≥3 streams.** The four-stream tumour-centric core:
VEGFA→KDR/FLT1/NRP1 · HLA-E→KLRC1 · HLA-F→LILRB1/LILRB2 · HLA-A/HLA-E→CD8A/CD8B ·
FGF2→FGFR1 · COL1A2→ITGA1/2/3/9_ITGB1 · NCAM1→FGFR1 · ANGPTL4→CDH5 · CD58→CD2.

### The architecture, in one paragraph

The MPNST tumour compartment is predicted to act principally on **myeloid cells and the
vasculature**. It supplies a dense, highly reproducible myeloid-directed signal set
(APP→CD74 across five receivers, CD99→PILRA, ANXA1→FPR1, HLA-F→LILRB1/2) whose published
counterparts in other cancers are immunosuppressive, while independently the myeloid
*transcriptional programme* is best explained by CSF1 — the axis already in clinical
development for neurofibroma. Toward the vasculature it is canonically angiogenic (VEGFA→KDR/
FLT1/NRP1, ANGPTL4→CDH5, COL1A2→endothelial integrins) and engages a **reciprocal perivascular
Notch circuit** with pericytes and endothelium, in a tumour type where Notch is already
implicated in Schwann-cell transformation. Toward lymphocytes the picture is **mixed rather
than uniformly suppressive**: class I presentation is retained and IL15 best explains the
CD8/NK/CD4 programmes, yet inhibitory HLA-E→NKG2A and HLA-F→LILRB engagement is among the
best-supported findings in the dataset. The fibroblast compartment is the weakest link:
abundant LR co-expression but **no receiver-programme evidence at all**.

## 18. Sample / patient robustness

Leave-one-patient-out retention of ≥2-patient support: 71.1% / 76.3% / 79.2% / 79.4% — **no
single patient dominates**. But **53.6% of supported interactions rest on one patient** and are
flagged, not promoted. Rare populations (Monocyte, B-cell, pDC, T-cell-other, NK) have a median
of 1–2 supporting patients and are systematically less reproducible; none was merged away.

## 19. Literature interpretation

`results/phase3/tables/CCC_LITERATURE_EVIDENCE.tsv` — 12 prioritised axes with DOIs.

**Established in MPNST/NF1 and recovered here:** NOTCH (JAG1/NOTCH2/NOTCH3 —
[doi:10.1093/noajnl/vdae188](https://academic.oup.com/noa/article/6/1/vdae188/7888919),
[Oncogene 2004](https://www.nature.com/articles/1207068)); CSF1→CSF1R (pexidartinib in
neurofibroma); SPP1→CD44 ([doi:10.1038/s41698-025-01078-2](https://www.nature.com/articles/s41698-025-01078-2)).
**Established elsewhere, apparently not characterised in MPNST:** APP→CD74
([doi:10.1002/path.6343](https://pathsocjournals.onlinelibrary.wiley.com/doi/full/10.1002/path.6343)),
ANXA1→FPR1, HLA-E→KLRC1, BAG6→NCR3. **Apparently less characterised anywhere here:** CD99→PILRA.

**No novelty is claimed on the basis of a search that failed to find a paper.**

## 20. Limitations

1. Communication **potential** only — no adjacency, no directionality proof.
2. `sample_id` = dataset = patient = batch; no between-tumour or condition claim is supportable.
3. **No condition variable exists**; no condition-stratified CCC was possible and none was fabricated.
4. 53.6% of supported interactions rest on one patient.
5. 97% of interaction keys exist in only one framework's resource.
6. LIANA and CellPhoneDB are partially non-independent.
7. LochNESS found no receiver-state structure and is itself unstable at n = 4.
8. With four samples the exact permutation null has a p-floor of 0.167.
9. **The Phase 2 malignant fraction is uncertain (17.35% conservatively, up to ~49%).** Every
   tumour-centric result inherits this; the tumour↔fibroblast axis is the most exposed.
10. NicheNet did not corroborate APP; the two evidence streams nominate different leaders.
11. NicheNet was run cohort-level, not sample-stratified.

## 21. Final figures — `results/phase3/figures/final/` (36 files)

`01_global_ccc_network` · `02_mpnst_outgoing_signaling` · `03_mpnst_incoming_signaling` ·
`04_tumor_to_immune_LR_heatmap` · `05_immune_to_tumor_LR_heatmap` · `06_method_concordance` ·
`07_method_overlap` · `08_sample_recurrence` · `09_receiver_ligand_activity` ·
`10_lochness_receiver_state` · `11_ccc_receiver_response_integration` ·
`12_final_mpnst_tme_network` · `13_lochness_implementation_comparison` ·
`14_lochness_null_permutation` · `15_leave_one_sample_out_robustness` ·
`16_method_disagreement` · `17_top_prioritised_interactions` · `18_ccc_input_cell_counts`
(all PDF + PNG). Per-milestone suites in `results/phase3/figures/M19..M26/` (108 files total),
all registered in `reports/FIGURE_INDEX.tsv` (578 rows).

## 22. Final tables — `results/phase3/tables/final/` (29 files)

`MPNST_CCC_MASTER_TABLE.tsv` (the primary deliverable) · `CCC_CONCORDANCE.tsv` ·
`method_overlap_counts.tsv` · `tumor_to_immune_interactions.tsv` ·
`immune_to_tumor_interactions.tsv` · `tumor_to_stromal_interactions.tsv` ·
`sample_recurrence.tsv` · `patient_recurrence.tsv` · `CCC_LITERATURE_EVIDENCE.tsv` ·
`CCC_INTEGRATED_EVIDENCE_TOP.tsv` · `NicheNet_ligand_activity.tsv` ·
`receiver_target_programs.tsv` · `lochness_scores.tsv` · `lochness_summary.tsv` ·
`lochness_permutation_summary.tsv` · `LOCHNESS_IMPLEMENTATION_COMPARISON.tsv` ·
`CCC_LEAVE_ONE_SAMPLE_OUT.tsv` · `CCC_METHOD_DISAGREEMENT_SUMMARY.tsv` and others.

The master table lets a researcher trace: interaction → methods supporting it → patients
supporting it → receiver-response evidence → LochNESS context → literature → priority tier.

## 23. Reproduction

```bash
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

sbatch scripts/shell/phase3/install_phase3_deps.sh          # then _retry.sh, _retry2.sh
sbatch scripts/shell/phase3/run_m19_ccc_inputs.sh           # M19
sbatch scripts/shell/phase3/run_m20_liana.sh                # M20 (parallel)
sbatch scripts/shell/phase3/run_m20_cellphonedb.sh          #     (cpdb_env)
for s in MPNST_1 MPNST_2 MPNST_3 MPNST_4; do
  sbatch --job-name=p3cc_$s scripts/shell/phase3/run_m20_cellchat_one.sh $s; done
sbatch scripts/shell/phase3/prefetch_nichenet3.sh
sbatch scripts/shell/phase3/run_m21_m22.sh                  # M21 + M22
sbatch scripts/shell/phase3/run_m23.sh                      # M23
sbatch scripts/shell/phase3/run_m24.sh                      # M24
sbatch scripts/shell/phase3/run_m25_m26.sh                  # M25 + M26
```

All seeds are 42. No Snakemake.

## 24. Environment

`R_env`: R 4.4.3, **Seurat 5.4.0, SeuratObject 5.3.0, harmony 1.2.4, Matrix 1.7.4 — all
unchanged from Phase 2 and verified after every install stage.** Added: liana 0.1.14,
CellChat 2.2.0.9001, nichenetr 2.2.1.1, OmnipathR 3.14.0, decoupleR 2.12.0, NMF 0.28,
FNN 1.1.4.1, systemfonts 1.3.2, ggraph 2.2.2, ggpubr 1.0.0.
**One dependency change recorded honestly: igraph 2.2.1 → 2.1.4**, downgraded by the conda
solver; Phase 2 outputs are frozen artefacts and unaffected.
`cpdb_env` (isolated, deliberate): Python 3.10, cellphonedb 5.0.1, anndata 0.11.4, scanpy 1.11.5.
Exports and full install log in `reports/phase3/environment/`.

## 25. SLURM accounting

| Stage | JobID | State | Elapsed |
| --- | --- | --- | --- |
| Network probe | 19894719 | COMPLETED | 00:00:02 |
| Install 1 / retry 1 / retry 2 | 19894761 / 19895124 / 19895459 | COMPLETED | 00:43:55 / 00:04:45 / 00:11:02 |
| M19 inputs | 19894939 | COMPLETED | 00:08:24 (27.43 GiB) |
| M20 LIANA | 19895133 FAILED → **19895159** | COMPLETED | 00:38:22 |
| M20 CellPhoneDB | 19895344 FAILED → **19895438** | COMPLETED | 00:02:14 |
| M20 CellChat | 19895564 CANCELLED → **19895738–41** | COMPLETED | 00:44:13 / 00:19:49 / 00:22:51 / 00:24:19 |
| NicheNet prefetch | 19895781, 19895842 → **19895947** | COMPLETED | resumable retry |
| M21 + M22 | 19896027 FAILED, ×1 more → **19896062** | COMPLETED | 00:05:42 |
| M23 NicheNet | 19896063 | COMPLETED | 00:03:13 |
| M24 LochNESS | 19896064 | COMPLETED | 00:04:14 |
| M25 + M26 | 19896065 | COMPLETED | 00:00:24 |

Peak memory anywhere: 27.43 GiB (M19) — 6% of the 450G envelope. **Every failure was diagnosed
from its log and fixed at the root cause; no resource was ever increased in response to a
failure, and no Phase 2 package was downgraded to make a tool work.**

## 26. Recommended Phase 4 spatial validation

**Not executed. Phase 4 requires separate authorization.**

Phase 3 outputs are structured so that a matched MPNST spatial dataset can test the specific
predictions:

1. **Sender/receiver proximity.** Are `MPNST-Tumor` cells physically adjacent to macrophages
   where APP→CD74 is predicted, and to pericytes/endothelium where the Notch circuit is
   predicted? The four-stream core (§17) is the priority list.
2. **Ligand localisation.** Is APP protein/transcript polarised toward CD74+ myeloid contacts?
   Is VEGFA concentrated at the vascular front?
3. **Receptor localisation.** Are KLRC1+ NK cells and LILRB1/2+ monocytes found where HLA-E and
   HLA-F are highest?
4. **Receiver-response localisation.** Do CSF1-proximal macrophages show the NicheNet-predicted
   target programme? This is the cleanest test of the CSF1 convergence.
5. **The fibroblast negative.** Spatial data could distinguish the two readings of §14 —
   non-functional LR co-expression versus a fibroblast programme driven by non-tumour inputs.
6. **The perivascular Notch circuit** is the most spatially specific prediction and the easiest
   to falsify: it requires physical tumour–pericyte–endothelial contact.

Also recommended before or alongside spatial work: **CNV inference**, which would resolve the
Phase 2 malignant-fraction uncertainty (§20.9) that propagates into every tumour-centric result
here.

**Not recommended on this dataset:** differential abundance and pseudobulk condition DE — four
patients, one sample each, no condition variable, hence no replication structure.
