# MPNST Single-Cell Tumour Ecosystem and Malignant-State Analysis

[![Phase 1 CI Workflow](https://github.com/vipinmenon1989/Pliot_Tumor_Research/actions/workflows/ci.yml/badge.svg)](https://github.com/vipinmenon1989/Pliot_Tumor_Research/actions/workflows/ci.yml)

A reproducible, configuration-driven single-cell RNA-seq study of four malignant peripheral nerve
sheath tumours (MPNST), running from independent per-sample preprocessing through integration and
annotation, tumour–microenvironment communication inference, copy-number-informed malignancy
refinement, continuous malignant transcriptional programs, and the relationship between
CNA-defined clones and transcriptional phenotype.

**Status: Phases 1–6 complete and frozen (2026-09-04). Phase 7 (spatial validation) is planned
and has not been performed.**

---

**Quick navigation** — [Overview](#1-project-overview) · [Scientific questions](#2-scientific-questions) ·
[Status](#3-project-status) · [Phase 1](#phase-1--independent-preprocessing-and-qc) ·
[Phase 2](#phase-2--integration-and-biological-annotation) ·
[Phase 3](#phase-3--tumourmicroenvironment-communication) ·
[Phase 4](#phase-4--cna-informed-malignancy-and-clone-structure) ·
[Phase 5](#phase-5--continuous-malignant-transcriptional-programs) ·
[Phase 6](#phase-6--clone--program-and-regulatory-architecture) ·
[Integrated model](#5-integrated-biological-model) ·
[Negative findings](#6-major-negative-findings) · [Methods](#7-computational-methods) ·
[Code map](#8-repository-and-code-map) · [Results map](#9-results-figures-and-tables) ·
[Reproducibility](#10-reproducibility-and-hpc-execution) · [Provenance](#11-provenance) ·
[Limitations](#12-limitations) · [Phase 7 (planned)](#13-phase-7--planned-spatial-validation-not-performed)

---

## 1. Project overview

MPNST is an aggressive soft-tissue sarcoma of Schwann-cell lineage, arising sporadically or in the
context of neurofibromatosis type 1. Its tumours are mesenchymal, ECM-rich and stroma-heavy, which
makes the boundary between malignant cells and genuine fibroblasts unusually hard to draw from
marker expression alone — and that boundary determines every downstream statement about what the
tumour compartment is and what it signals to.

**Cohort.** Four MPNST tumours from a previously processed clinical Seurat object
(`processed_mpnst.rds`), analysed as four independent datasets and carried through every phase as
**19,716 post-QC cells**:

| Sample | Cells | Note |
| --- | ---: | --- |
| MPNST_1 | 7,615 | deepest library (median 5,361 genes/cell in the malignant compartment) |
| MPNST_2 | 2,284 | |
| MPNST_3 | 2,940 | lowest depth; **SCEVAN clone inference unreliable** (Phase 4) |
| MPNST_4 | 6,877 | largest malignant compartment; carries much of the cohort-level effect |

**A structural constraint that governs the whole project:** `sample_id` is simultaneously the
dataset, the patient and the presumed technical batch, and **no clinical or condition covariate
exists**. Sample and patient are therefore the same unit throughout this README. With n = 4, a
patient-private *biological* signal cannot be fully separated from patient-level *technical*
structure, and no differential-abundance or condition-level test has any replication structure to
rest on. Both were deliberately not attempted.

**Trajectory of the study.** Phase 1 built a frozen non-integrated baseline; Phase 2 integrated the
four tumours and annotated them; Phase 3 inferred the tumour–microenvironment communication
architecture; Phase 4 brought orthogonal copy-number evidence to bear on malignant identity and
clone structure; Phase 5 asked whether malignant transcriptional *programs* recur across patients
when discrete states do not; Phase 6 asked whether those programs are constrained by CNA-defined
clones, and what regulatory and pathway architecture underlies them.

---

## 2. Scientific questions

1. Which cell populations compose the MPNST tumour ecosystem, and can they be integrated across
   patients without erasing between-tumour biology?
2. Which tumour ↔ microenvironment communication programs are reproducible across independent
   inference frameworks and across patients?
3. Does expression-derived copy-number inference refine the malignant compartment beyond what
   marker-based annotation can settle?
4. Are the fibroblast-like / ECM-rich populations partly malignant?
5. Are discrete malignant transcriptional states recurrent between patients?
6. Do *continuous* malignant programs recur where discrete states do not?
7. Are malignant transcriptional programs constrained by CNA-defined genetic clones?
8. Do the tumours converge at regulatory, pathway or microenvironmental level despite
   patient-specific transcriptional programs?

Answers, in one line each: a 17-identity ecosystem integrates well in its immune and vascular
compartments (Q1); a myeloid- and vasculature-directed communication architecture is reproducible
and survives relabelling (Q2); yes, substantially (Q3); yes — CNA-informed analysis supports a large
malignant mesenchymal/ECM-like compartment previously annotated as fibroblast-like (Q4); **no**
(Q5); **no** (Q6); **only weakly** (Q7); partially, at the regulatory, pathway and TME levels (Q8).

---

## 3. Project status

| Phase | Milestones | Scope | Status |
| --- | --- | --- | --- |
| **Phase 1** | M0–M9 | Independent preprocessing, QC, pre-integration baseline, reproducibility framework | complete / frozen |
| **Phase 2** | M10–M17 (+M15A) | Harmony integration, clustering, markers, hierarchical annotation, `annotation_ccc` | complete / frozen |
| **Phase 3** | M18–M27 | Tumour–microenvironment communication (LIANA · CellChat · CellPhoneDB · NicheNet · LochNESS) | complete / frozen |
| **Phase 4** | M28–M35 (+M35A) | SCEVAN CNA-informed malignancy refinement, clones, tumour states, CCC sensitivity | complete / frozen |
| **Phase 5** | M36–M41 | Continuous malignant programs (cNMF); malignant ECM-like vs true fibroblast | complete / frozen |
| **Phase 6** | M42–M50 | Clone ↔ program coupling, within-clone diversity, TF/pathway architecture | complete / frozen |
| **Phase 7** | — | **Spatial validation — planned, not initiated** | not started |

Frozen final objects (all 19,716 cells; each phase **adds** metadata and overwrites nothing):

| Phase | Final object | md5 |
| --- | --- | --- |
| 1 | `results/combined/pre_integration/combined_preintegration.rds` | `88a442688f912d882f6c6da01820e329` |
| 2 | `results/phase2/phase2_final_object.rds` | `153d5f6acc70f9c05aa48cabc4f4ac2d` |
| 4 | `results/phase4/phase4_final_object.rds` | `e85ba8486e456917e2483f2773bdbaf3` |
| 5 | `results/phase5/phase5_final_object.rds` | `839e5157bc7c3470ddf86746c2e719e1` |
| 6 | `results/phase6/phase6_final_object.rds` | `b8c01dd01755dda10b2be4e0f5ef7ef7` |

Phase 3 produced tables and figures rather than a new Seurat object; its CCC input object is
`results/phase3/ccc/ccc_input_object.rds` (md5 `aeb02c59df8a42a24074df509c662f7f`). sha256 digests,
sizes and per-milestone provenance live in each phase manifest and handoff.

---

## 4. Phase-by-phase summary

Each phase has a self-contained handoff document that can be read without any conversation
history. The summaries below are deliberately compressed; the handoffs carry the evidence.

### Phase 1 — Independent preprocessing and QC
*Handoff: [`reports/PHASE1_HANDOFF.md`](reports/PHASE1_HANDOFF.md) · milestones
[`reports/milestones/`](reports/milestones/)*

**Method.** Audit of the source object (M1); deterministic per-sample extraction and pre-filter QC
(M2); QC filtering with `scDblFinder` doublet removal and dataset-specific thresholds — ≥500
features, ≥1,000 counts, ≤15% mitochondrial, ≤20% ribosomal (M3); `SCTransform` v2 with 3,000
variable features (M4); independent PCA with elbow-based PC selection (M5); Louvain clustering
resolution sweep 0.1–1.0 with ARI stability assessment (M6); marker discovery at every resolution
(M7); a combined **non-integrated** baseline in a shared expression space (M8); workflow hardening,
CI, provenance and freeze (M9).

**Result.** 19,716 cells retained. Recommended PCs/resolutions: MPNST_1 1–8 / 0.6 (18 clusters),
MPNST_2 1–6 / 0.3 (9), MPNST_3 1–9 / 0.6 (13), MPNST_4 1–5 / 0.7 (14). In the shared non-integrated
PCA the mean same-dataset neighbour fraction is **96.13%**.

**Interpretation.** Strong dataset-associated structure exists before integration, so integration
had to be evaluated against this frozen baseline rather than assumed to work.

**Limitation.** Dataset identity is fully confounded with patient identity; no capture-channel or
run metadata exists. `percent.mt` is identically zero in all MPNST_1 cells (verified to produce no
non-finite values downstream).

Phase 1 is also where the reproducibility framework lives: Snakemake orchestration, synthetic-data
smoke testing, unit tests, SLURM wrappers, JSON provenance and a preflight checker. **Phases 2–6 do
not use Snakemake** — they are driven by explicit per-milestone SLURM scripts.

---

### Phase 2 — Integration and biological annotation
*Handoff: [`reports/phase2/PHASE2_HANDOFF.md`](reports/phase2/PHASE2_HANDOFF.md) · milestones
[`reports/phase2/milestones/`](reports/phase2/milestones/)*

**Method.** Harmony 1.2.4 on `sample_id` (`pca` dims 1:30 → `postint_harmony`), **every scientific
parameter at package defaults**, seed 42, bit-reproducible on re-run. Harmony was used for
**embedding only** — never as an expression basis for testing. Integration was then assessed
against the Phase 1 baseline on all 19,716 cells (M12), clustered at resolutions 0.1–1.0 with all
ten preserved (M13), characterised by markers via `PrepSCTFindMarkers()` + `FindAllMarkers` on
`SCT`/`data` (M14), and annotated hierarchically from an evidence table rather than in code —
`config/phase2/annotation_map_M15.tsv` carries every label with its positive/negative markers,
confidence and citation (M15/M16).

**Result.** Primary resolution 1.0 → **26 clusters** (alternative 0.7 → 21; every adjacent-resolution
ARI ≥ 0.89). 35,437 marker rows, 31,774 significant, markers for every cluster. Eight Level-1
compartments and **17 Level-2 cell types**. Mixing is compartment-specific: the same-sample neighbour
fraction falls from 0.933→0.527 (panleukocyte), 0.958→0.611 (myeloid), 0.953→0.689 (T/NK),
0.962→0.709 (endothelial), but only 0.990→0.900 for Schwann/neural-crest — Harmony aligned the shared
immune and vascular compartments while leaving the presumptive malignant lineage patient-private.

**The CCC-oriented layer (M15A).** A second annotation layer, `annotation_ccc`, collapses the
*evidence-supported* malignant states into a single `MPNST-Tumor` identity while keeping nine
immune, two stromal and one endothelial identity separate — 17 identities in total, of which 14 are
CCC-ready. The collapse rule is data, not code
(`config/phase2/ccc_annotation_map.tsv`).

> **`MPNST-Tumor` was NOT defined as "everything that is not immune."** It required positive
> MPNST/Schwann/neural-crest marker evidence. `MPNST-Tumor` = **3,420 cells (17.35%)**;
> `Fibroblast` = 5,064 (25.68%); `Candidate-Malignant-Unresolved` (1,231), `Uncertain` (720) and
> `Low-quality-excluded` (438) were retained but flagged **not CCC-ready**.

**Interpretation / open question.** Under this conservative annotation the malignant compartment is
23.6% (detailed) or 17.35% (`annotation_ccc`). Phase 2 recorded explicitly that if the four
fibroblast clusters are Mes-NC-like malignant, the true compartment could reach ~49%, and that
**resolving this requires CNA inference** — which became Phase 4.

**Limitation.** Batch and patient are one variable, so no between-tumour claim is supportable from
the integrated embedding; B/plasma and fibroblast cohesion each fell ~44% under Harmony (caveats C1,
C2); four of 17 cell types are >80% one patient; 16.5% of cells carry Low confidence.

---

### Phase 3 — Tumour–microenvironment communication
*Handoff: [`reports/phase3/PHASE3_HANDOFF.md`](reports/phase3/PHASE3_HANDOFF.md) · method plan
[`PHASE3_METHOD_PLAN.md`](reports/phase3/PHASE3_METHOD_PLAN.md) · milestones
[`reports/phase3/milestones/`](reports/phase3/milestones/)*

> **Interpretation rule for this entire phase.** Ligand–receptor analysis infers communication
> **potential** from expression and associated transcriptional programs. It does not establish cell
> adjacency or physical signalling. Every statement below is *predicted / inferred / supported /
> concordant with*, never *demonstrated*.

**Method — five tools with five distinct roles, not a benchmark.** 14 CCC-ready populations, 17,327
cells, computed **per sample** and only then compared; expression basis RNA / `LogNormalize`
(Harmony coordinates were never used as expression).

* **LIANA** 0.1.14 — primary multi-resource consensus ligand–receptor framework.
* **CellChat** 2.2.0.9001 — independent pathway/network framework.
* **CellPhoneDB** 5.0.1 — independent permutation-based LR framework.
* **NicheNet** 2.2.1.1 — **receiver transcriptional response**, not a fourth LR list.
* **LochNESS** (MMCA formulation) — an **orthogonal neighbourhood-phenotype analysis**, not a
  communication method.

**Raw scores were never averaged across frameworks.** Each LR framework contributes a boolean
support flag; patient recurrence and expression support enter as separate evidence streams, and
LIANA's partial non-independence from CellPhoneDB is declared rather than hidden.

**Result.** 502,846 distinct interaction keys, **36,486 supported**, of which 5,848 high-concordance
and 5,485 discordant (reported, not averaged away). Of the 1,347 interactions jointly testable by
all three LR frameworks, 763 (57%) were supported by all three. **547 tumour-centric interactions
reach ≥3 independent evidence streams**; 325 interactions reach all four.

The predicted architecture, by receiver compartment:

| Direction | Predicted axes |
| --- | --- |
| Tumour → myeloid | **APP–CD74** (macrophage, dendritic, monocyte, B, pDC — the most consistent tumour signal in the dataset) · **CD99–PILRA** · **ANXA1–FPR1** · THBS1–CD36 · **HLA-F–LILRB1/LILRB2** |
| Tumour ↔ lymphoid, deliberately **mixed** | inhibitory **HLA-E–KLRC1 (NKG2A)** and HLA-F–LILRB1/2 alongside activating **CD58–CD2** and **BAG6–NCR3 (NKp30)**; class I presentation retained (HLA-A/HLA-E–CD8A/CD8B) |
| Tumour → vascular / perivascular | **VEGFA–KDR/FLT1/NRP1** (4 streams, 4/4 patients) · **ANGPTL4–CDH5** · COL1A2–endothelial integrins · **JAG1–NOTCH3/NOTCH4** |
| Stroma/vasculature → tumour | **JAG1/JAG2/DLL4–NOTCH2** (a reciprocal perivascular Notch circuit) · COL1A1/COL1A2/COL6A2/FN1–**ITGAV_ITGB8** · CRTAM–CADM1 · TNF–TNFRSF1A · NCAM1–FGFR1 · SEMA4D–PLXNB2 |

**Receiver response (NicheNet).** Myeloid receivers converge on **CSF1**; lymphoid receivers on
**IL15**; endothelial on TGFB1/HMGB1/VEGFA/ANGPT1. **A reported divergence:** APP is *not* among the
top NicheNet ligands for macrophages, which bounds the APP–CD74 finding to *predicted engagement*
rather than *drives the macrophage state*. The **fibroblast receiver is a clear negative** (best AUPR
0.022): tumour ligands do not explain the fibroblast program.

**Negative result (LochNESS).** Sample-level context = tumour **APP** expression, with a
same-sample-excluded neighbourhood score. **No lineage shows receiver-state structure associated
with the APP context beyond chance** (best descriptive p = 0.333 against a 0.167 exact-permutation
p-floor at n = 4), and the score itself is unstable across embeddings and reference samples.
LochNESS therefore **did not independently corroborate APP–CD74 as a receiver-state driver**. This
does not refute the LR finding — LochNESS is not a communication method. The R (MMCA) and Python
implementations agree exactly (Pearson 1.0000, max abs diff 0) on equivalent inputs.

**Limitation.** Communication potential only; 53.6% of supported interactions rest on one patient;
97% of interaction keys exist in only one framework's resource, so "single-method" usually means
resource non-overlap rather than disagreement; and every tumour-centric result inherits the Phase 2
malignant-fraction uncertainty — the exposure that motivated Phase 4.

---

### Phase 4 — CNA-informed malignancy and clone structure
*Handoff: [`reports/phase4/PHASE4_HANDOFF.md`](reports/phase4/PHASE4_HANDOFF.md) · decision rules
[`MALIGNANCY_DECISION_RULES.md`](reports/phase4/MALIGNANCY_DECISION_RULES.md) · milestones
[`reports/phase4/milestones/`](reports/phase4/milestones/)*

> **SCEVAN infers copy-number architecture from gene-expression patterns. This is not DNA
> sequencing.** Inference is meaningful for broad chromosomal, arm-level and large-segment events
> only. **No single-gene CNV call is asserted anywhere in this project.** A SCEVAN non-malignant
> call is also **not proof** of non-malignancy — some malignant cells may be copy-number quiet.

**Method.** SCEVAN 1.0.3 (De Falco *et al.*, *Nat Commun* 14:1074, 2023) with yaGST 2017.8.25, run
**per patient, twice, with identical parameters** on raw integer RNA counts, seed 42:

* **Primary run — `norm_cell = NULL`**: SCEVAN's confident-normal detection ran with **no input from
  Phase 2 labels at all**. This is the run that carries inferential weight.
* **Sensitivity run** — high-confidence immune cells as the normal reference.
* Endothelium was deliberately **excluded** from the reference set so its malignancy call stayed
  independent of the Phase 3 vascular findings. `FIXED_NORMAL_CELLS = TRUE` is prohibited in this
  project because it forces every non-reference cell to malignant.
* **20 ordered decision rules with every threshold fixed a priori were committed before any SCEVAN
  result was inspected.** Two later amendments (A1, A2) are dated and labelled post-hoc.

**Result — classification and reliability.**

| Patient | cells | assessed | SCEVAN malignant | subclones | primary↔sensitivity agreement |
| --- | ---: | ---: | ---: | ---: | ---: |
| MPNST_1 | 7,615 | 7,043 | 1,952 | 7 | 0.9951 |
| MPNST_2 | 2,284 | 2,205 | 694 | 4 | 0.9663 |
| MPNST_3 | 2,940 | 2,749 | (1,414) | (3) | **0.0864 — failed** |
| MPNST_4 | 6,877 | 6,452 | 4,734 | 8 | 0.9974 |
| **total** | **19,716** | **18,449** | **8,794** | **22** | — |

Filtered cells are recorded as *not assessed*, never as non-malignant.

**MPNST_3 failed a sanity check and was handled, not hidden.** Both its runs assign large blocks of
unambiguous immune cells to the malignant class **in opposite directions**, its three "clones"
resolve to pure immune lineages, and only 25 confident normal cells were found. Malignant
*promotions* from MPNST_3 were disabled, and the cost is stated plainly: any malignant
mesenchymal population in MPNST_3 cannot be detected by this analysis. Its clone structure is
**excluded from every clone-based inference in Phases 4 and 6**.

**Result — broad CNA architecture.** 59 broad clonal segments (≥10 Mb). Recurrent among the three
reliable patients: **chr18 loss (3/3), chr2 gain (3/3), chr7 gain (3/3)**. Pairwise Jaccard of broad
event sets is only **0.20–0.23** and 1 of 59 events appears in all four samples — CNA architecture is
**predominantly patient-specific with a small shared broad-event core**. These are described as
shared or recurrent *CNA patterns*, never as "the same clone" across patients.

> **On NF1 and NF2:** MPNST_4 carries a 38.1 Mb single-copy segment on chr17 **whose interval
> contains the NF1 locus**, and MPNST_1 and MPNST_4 each carry a 34.2 Mb single-copy segment on 22q
> **containing NF2**. This analysis cannot distinguish a focal deletion from a large-segment loss
> that happens to include the locus and does not attempt to. **No NF1 or NF2 deletion is claimed.**

**Result — refined malignancy.** `malignancy_refined` ∈ {Malignant, Non-malignant, Ambiguous,
Excluded-low-quality}, each with a confidence grade and a printed rule and reason. `Ambiguous` is a
**terminal outcome and never becomes tumour**.

| Stage | Malignant | Fraction of 19,716 |
| --- | ---: | ---: |
| Phase 2 conservative | 3,420 | 17.35% |
| **Phase 4 refined** | **6,434** | **32.63%** (High 3,261 / Moderate 3,173) |

Plus 9,078 Non-malignant, **3,766 Ambiguous** and 438 Excluded. The correction is two-directional:
Phase 4 promoted 4,036 fibroblast-labelled, 836 candidate-malignant, 92 uncertain and 65 pericyte
cells, while **withdrawing confidence from 2,015 of the 3,420 cells Phase 2 called malignant** (moved
to `Ambiguous` — lineage and copy-number evidence genuinely disagree there, and the label records
that rather than picking a winner).

**Result — clones and tumour states.** 22 patient-scoped subclones; no patient is dominated by one
subclone (largest holds 24–34% of the malignant compartment in the reliable patients), so the
reliable tumours are genuinely **polyclonal**. Eight malignant transcriptional states were resolved
from a malignant-only reduction, and **0 of 8 are recurrent** across patients (≥3 patients each
contributing ≥5% and ≥10 cells); **97.2% of malignant cells sit in patient-private states**. Clone
counts per state range from 2 to 19, so state and clone are kept as separate concepts and reported
as an association, never an equivalence.

**Result — the Phase 3 architecture survives relabelling.** Phase 3's LIANA, CellChat, CellPhoneDB
and concordance scripts were reused **completely unmodified** at identical versions, resources,
thresholds and seed; only `ccc_label` differs, which is what makes this a label-sensitivity analysis
rather than a tool-drift confound. High-concordance interactions fell **−2.6%**, the median change
across 23 named axes was **−3.1%**, and **no axis was lost**. What changed is **attribution**: 3,051
interactions were sender-reassigned, concentrated in the ECM→integrin axis, so the
collagen/FN1–ITGAV_ITGB8 signalling Phase 3 read as stroma-to-tumour is substantially
**tumour-autocrine**.

#### The fibroblast-like malignant compartment — the central Phase 4 finding

Of the **5,064** cells Phase 2 labelled `Fibroblast`:

| Refined call | Cells | Share |
| --- | ---: | ---: |
| **Malignant** | **4,036** | 79.7% |
| Non-malignant | 908 | 17.9% |
| Ambiguous | 120 | 2.4% |

Four independent legs support this, and the strongest is **not** the malignant fraction itself:

1. **Clone composition — the least circular evidence.** SCEVAN builds subclones from inferred CNA
   profiles **without ever seeing a Phase 2 label**, so asking what a CNA-defined clone is made of
   is not a restatement of the malignancy call. **All 4 of MPNST_2's inferred malignant subclones
   are fibroblast-dominated** (229/234, 160/186, 138/158, 110/116). **7 of MPNST_4's 8 subclones are
   fibroblast-dominated.** In **MPNST_1, all 7 clones mix `Candidate-Malignant-Unresolved` with
   `MPNST-Tumor`, and 6 of the 7 additionally carry `Fibroblast` cells** — marker-based annotation
   split into three identities what copy-number evidence sees as one malignant compartment with
   seven subclones.
2. **Threshold independence.** Fibroblast → Malignant is **4,036 at every tested `pop_frac_low`**
   (0.15, 0.20, 0.25, 0.30, 0.40).
3. **Reference independence.** Primary/sensitivity agreement 0.966–0.997 in the reliable patients.
4. **Effect size and CNA profile shape**, reported per patient because cells within a patient are
   not independent: Cliff's delta **0.78–0.92** separating malignant-called from non-malignant
   fibroblast-labelled cells, and a genome-wide CNA profile correlating **0.934** (MPNST_1) and
   **0.972** (MPNST_4) with the malignant compartment against **0.184** with non-malignant
   fibroblasts. **This contrast is evaluable in only 2 of 4 patients** — MPNST_4 retains a single
   non-malignant fibroblast and MPNST_3 contributes no malignant fibroblasts — and the figures state
   that rather than hiding it.

**Conservative rejection was applied where the evidence pointed the wrong way.** MPNST_4's
endothelium is 52.7% SCEVAN-malignant and forms its own subclone (`MPNST_4_clone8`, 146/152
endothelial); rule R11 sent those cells to `Ambiguous/Low` on the reasoning that a
majority-malignant canonical endothelial population in **one** patient is more likely an artefact
than a discovery. MPNST_3's immune-lineage "clones" were rejected outright. 1,125 plasma cells with
an immunoglobulin-locus artefact were called `Non-malignant` with the artefact named rather than
dropped.

**Interpretation.** CNA-informed analysis supports a **substantial malignant mesenchymal/ECM-like
population that marker-based annotation had read as fibroblast-like**. Phase 2 found the
Schwann-like and cycling malignant cells and missed the ECM-programme ones, because an ECM
programme is what a fibroblast looks like. That is why copy-number evidence, rather than a better
marker panel, was the right instrument.

**Limitation, stated with the finding.** The **32.63% cohort fraction is not patient-robust** —
dropping MPNST_4 returns it to the Phase 2 value (0.214 vs 0.215). The *fibroblast* conclusion is a
separate claim and is robust in direction; its *magnitude* depends on amendment A2, which was written
after observing that the flat sanity gate excluded the sample carrying the strongest evidence. Both
gates are published side by side.

#### M35A — SCEVAN evidence consolidation (visualization only)
*Audit: [`reports/phase4/M35A_SCEVAN_FIGURE_AUDIT.md`](reports/phase4/M35A_SCEVAN_FIGURE_AUDIT.md) ·
handoff §28*

M35A is **not a new scientific analysis**. SCEVAN was not re-run; no malignancy call, decision rule,
amendment, clone assignment, tumour state, CCC result or threshold changed, and the frozen final
object was neither loaded nor modified. The CNA evidence existed in tables and native SCEVAN output
but was not visible in the figure suite.

* Native SCEVAN CNA heatmaps and consensus profiles for the three reliable patients are surfaced
  **verbatim as rasters** (figures `17`–`20`) — no CNA value was redrawn anywhere — each beside a
  companion panel built from the same frozen clone assignments.
* Eight custom figures (`21`–`28`) carry the argument: clone composition, per-patient CNA burden
  effect sizes, CNA profile correlations, the Phase 2 → Phase 4 identity transition, **MPNST_3 shown
  failing rather than omitted**, threshold robustness (robust fibroblast count vs fragile cohort
  fraction), broad CNA recurrence at segment resolution, and a one-page summary.
* Native output is indexed in `results/phase4/tables/final/SCEVAN_NATIVE_FIGURE_INDEX.tsv`.
* Cytoband onco-heatmaps were deliberately **not** surfaced: they are unreadable at page scale and
  invite gene-level over-reading. All eight `CloneTree.png` files are blank — a documented
  ggtree/ggplot2 4.x failure; ggplot2 was not downgraded and no clone data is affected.
* One discrepancy against the handoff text was found and **recorded rather than smoothed** (§28.7):
  the claim that *every* MPNST_1 clone mixes all three disputed identities holds for 6 of 7 clones.

---

### Phase 5 — Continuous malignant transcriptional programs
*Handoff: [`reports/phase5/PHASE5_HANDOFF.md`](reports/phase5/PHASE5_HANDOFF.md) · robustness
[`reports/phase5/ROBUSTNESS_REPORT.md`](reports/phase5/ROBUSTNESS_REPORT.md) · milestones
[`reports/phase5/milestones/`](reports/phase5/milestones/)*

> A cNMF program is a **continuous pattern of co-regulated expression**. It is not a cell type, a
> lineage, a discrete state or a clone, and those concepts stay separate throughout.

**The question.** Phase 4 ended on a negative result — 0 of 8 *discrete* malignant states recurrent.
Phase 5 does not re-cluster to rescue that; it asks the *continuous* version, which is a genuinely
different object: a cell belongs to one cluster but carries a loading on every program, so two
patients can share a program while their cells cluster apart.

**Method.** **cNMF 1.7.1** in the isolated `p5_cnmf_env`, on the **6,434 refined malignant cells**
only, from RNA raw integer counts (never Harmony, UMAP, PCA, SCT residuals or the CNA matrix),
19,663 genes detected in ≥0.5% of malignant cells with `^MT-` removed, K grid 4–15, 100 replicates
per K, seed 42. **Cell cycle, ECM, HLA, interferon, Schwann, neural-crest and angiogenesis genes
were all retained** — removing them would delete the signal the phase exists to find. The
**K-selection rule was declared in code before it ran** and references only measured properties of
the factorization: **K = 8**.

**The eight programs** (exact final labels; three are technical-dominated and are named for it):

| id | label | top genes | technical? |
| --- | --- | --- | :--: |
| P1 | `Translation_ribosomal` | RPL10 RPL19 RPS2 RPLP1 EEF1A1 | **yes** (54% ribosomal) |
| P2 | `Unassigned_NEURONAL_SYSTEM` | TSHZ2 AUTS2 TENM2 AFF3 CHN1 NLGN1 | no |
| P3 | `Mesenchymal_ECM_mixed` | SERPINF1 CST3 FBLN1 SFRP2 EFEMP1 MFAP5 | no |
| P4 | `Hypoxia_Angio` | GBE1 IGF1R VEGFA ADAMTS9-AS2 | no |
| P5 | `Translation_ribosomal` (pseudogene) | RPS7P1 RPL3P4 | **yes** (66% pseudogene) |
| P6 | `Schwann_like_mixed` | S100B GPM6B SERPINE2 CRYAB MIA | no |
| P7 | `Cycling` | NUSAP1 PRC1 TOP2A UBE2C MKI67 | no |
| P8 | `Myeloid_ambient_like` | C1QC C1QB C1QA CD14 TYROBP CSF1R | **yes** (18 myeloid markers) |

A first labelling pass called P1 `Mesenchymal_ECM` while its own top 15 genes were ribosomal
proteins; an explicit technical-content check was added and M38 re-run. **Only the naming changed** —
not the factorization, K, or the recurrence rule.

**Headline result — a negative one, at the pre-declared threshold.** Recurrence criteria were
declared before any program was labelled, deliberately the same shape as Phase 4's rule
(program-active at ≥0.20 relative usage; a patient *carries* a program with ≥10 active cells and ≥5%
of its malignant cells; recurrent needs ≥3 of 4 patients):

> **0 recurrent · 1 shared-limited · 7 patient-private.** The single shared-limited program is P1 —
> one of the technical ones. **Not one biologically interpretable program is carried by even two
> patients at the declared threshold.**

This survives the most aggressive available check: a full re-run with ribosomal, pseudogene/lncRNA
and ambient-myeloid genes removed across the whole K grid gives **0 recurrent programs at every K
from 5 to 15**. Reported honestly alongside it: at a relaxed **0.10** activity threshold, **P3
(`Mesenchymal_ECM_mixed`) reaches three patients** — the one biologically interpretable program that
comes close. The honest statement is exactly that: recurrent at 10%, not at the pre-declared 20%.

**An important within-patient result.** The five Phase 4 `Mesenchymal_ECM` states are **not** one
shared program — they map to three. But *within* MPNST_4, three separately-clustered ECM states
(3,383 cells) **collapse onto a single continuous program**, so Phase 4's discrete partition was
finer there than the continuous structure warrants. **Clustering can over-partition continuous
transcriptional structure**, and this is a concrete instance.

#### Malignant ECM-like cells versus true fibroblasts

The grouping variable had to be the **pre-Phase-4** annotation, so `annotation_ccc_refined` was
**barred** from that role (it already encodes the Phase 4 conclusion); every candidate column was
tested against the frozen split and `annotation_ccc_phase3` was selected **because it qualified**.

| Group | Cells | Role |
| --- | ---: | --- |
| A — historical fibroblast, refined **Malignant** | **4,036** | compared |
| B — historical fibroblast, refined **Non-malignant** | **908** | compared |
| C — historical fibroblast, **Ambiguous** | **120** | projection/description only; never reclassified |

**Evaluability was checked before any test was run** (≥30 cells on both sides): evaluable in
**MPNST_1** (512 vs 303) and **MPNST_2** (637 vs 503); **not evaluable** in MPNST_3 (0 malignant
fibroblasts, promotions disabled) or MPNST_4 (1 non-malignant fibroblast). **A four-patient paired
test was therefore not constructed.** The cross-patient Spearman correlation of the pseudobulk
log2FC is only **0.091** — most of the A-vs-B difference is patient-specific, and that number bounds
the whole comparison.

A transparent **92-gene candidate signature** was derived by declared criteria applied in *every*
evaluable patient (|pseudobulk log2FC| ≥ 1, cell-level AUC ≥ 0.65, detected in ≥10% of the higher
group) — **no classifier, no train/test split, no accuracy claim**: 61 genes up in malignant
ECM-like (CA12, SCG2, IGFBP3, TMEM176A/B, COL14A1, COL11A1, GJA1, SPOCK1 …) and 31 up in true
fibroblast (**CDH19, APOD, SCN7A, ABCA6/8/9/10**, VIT, PAMR1, SPARCL1, LAMA2 …). The fibroblast side
is coherent in a way worth naming — those genes mark **nerve-associated / endoneurial fibroblasts
and non-myelinating Schwann-lineage stroma**, exactly the resident population a peripheral-nerve
sheath tumour would retain. **The signature is descriptive and evaluable in 2 of 4 patients; it is
not a validated marker panel.**

**Robustness.** Programs are stable to rank (median cosine 0.998 at K ± 1), patient balance (median
0.940), cell-cycle removal (median 1.000) and confidence filtering; no program's usage tracks a
within-patient technical covariate above |ρ| = 0.37. Every patient-private program vanishes when its
own patient is withheld (LOO cosine 0.17–0.35), while P6 and P7 — the two whose dominant patient is
not the one defining them — survive at 0.80 and 0.88.

---

### Phase 6 — Clone ↔ program and regulatory architecture
*Handoff: [`reports/phase6/PHASE6_HANDOFF.md`](reports/phase6/PHASE6_HANDOFF.md) · robustness
[`reports/phase6/ROBUSTNESS_REPORT.md`](reports/phase6/ROBUSTNESS_REPORT.md) · milestones
[`reports/phase6/milestones/`](reports/phase6/milestones/)*

**The central question.**

> Are malignant transcriptional phenotypes constrained by CNA-defined genetic clones, is there
> substantial program diversity *within* clones, or is it a mixture?

**Method.** Patient-scoped SCEVAN clones from the frozen Phase 4 object — **19 clones across
MPNST_1/2/4, 18 evaluable** at the declared 20-cell minimum; **MPNST_3 excluded**, with the exclusion
**re-derived** from the frozen reliability flag rather than taken on trust. Clone labels are
patient-scoped: `MPNST_1_clone1` and `MPNST_4_clone1` are unrelated names and were never compared as
homologous. 178 malignant cells carry no clone label and are excluded with the number stated.

Every clone × cNMF-program association is computed **inside one patient**, against a permutation null
built by shuffling clone labels within that patient (1,000 permutations), which preserves clone sizes
and the score distribution. Effect size is **η²**. Regulatory layer: decoupleR 2.12.0 `run_ulm` over
**CollecTRI** (41,674 edges, 1,201 TFs) on the frozen log-normalised RNA layer. Pathway layer:
**PROGENy** (14 footprint pathways, `run_mlm`) and **MSigDB Hallmark** (50 sets, mean z) — reported
**side by side and never merged** into a composite score.

**Central result — transcriptional phenotype is largely decoupled from clone structure.**

| patient | clones | median η² | max η² | clone-associated programs |
| --- | ---: | ---: | ---: | ---: |
| MPNST_1 | 7 | 0.072 | 0.508 | 3 of 8 |
| MPNST_2 | 4 | 0.050 | 0.154 | 2 of 8 |
| MPNST_4 | 7 | 0.072 | 0.208 | 2 of 8 |

> **Median clone→program η² = 0.059 — roughly 94% of each program's variance sits WITHIN
> CNA-defined clones**, and between-clone divergence is only **0.7%** of within-clone dispersion. A
> patient's CNA-defined clones are transcriptionally near-interchangeable.

21 of 24 program × patient pairs reach permutation p < 0.001, which is exactly why **effect size and
not the p-value carries the conclusion**: with 651–3,623 cells almost any non-zero η² is
"significant". A patient's dominant program is not the property of one clone — MPNST_1's P2 is
active in 7 of 7 clones, MPNST_2's P3 in 4 of 4, MPNST_4's P1 in 7 of 7.

**Within-clone diversity is real but strongly patient-specific.** Effective programs per clone:
MPNST_1 **1.69–2.83**, MPNST_2 **1.00** (all four clones), MPNST_4 1.05–1.23.

> **Permitted statement:** four clones, all in MPNST_1, **contain cells spanning multiple malignant
> transcriptional programs**, which is **consistent with transcriptional plasticity**. It is **not**
> observed switching — no transition rate, direction or trajectory is claimed anywhere, because
> nothing longitudinal was measured. Program *dispersion* is additionally confounded with the
> fraction of High-confidence malignant cells in a clone (ρ = 0.72, above the declared disqualifying
> bar), so any plasticity reading of dispersion must carry that number; the two effective-number
> metrics stay below it (0.56–0.62).

**Regulatory architecture, recovered rather than imposed.** The phase brief listed AP-1, TEAD/YAP,
STAT/IRF, E2F, MYC, NF-κB and SOX-family as **hypotheses only**. E2F and MYC emerged for the cycling
program and HIF1A/ATF4/HSF1 for the translation-stress program **without being imposed**; TEAD/YAP
and SOX-family did not emerge and are not claimed.

| Program | Concordant TF associations | Independent pathway support |
| --- | --- | --- |
| P7 `Cycling` | E2F4 0.43, MYC 0.34, E2F1 0.24 | Hallmark E2F_TARGETS 0.43, MYC_TARGETS_V1 0.44 |
| P5 translation-stress | HIF1A 0.45, ATF4 0.43, HSF1 0.45, CREB1 0.43 | MTORC1 0.58, GLYCOLYSIS 0.55, PI3K-AKT 0.50 |
| P3 `Mesenchymal_ECM` | HMGA2 0.57, MYC 0.56, SP1 0.56, STAT6 0.54 | Hallmark EMT 0.33; PROGENy EGFR 0.57 |
| P4 `Hypoxia_Angio` | ID4 0.45, PGR 0.44, TBX2 0.43 | **PROGENy Hypoxia is its top pathway, 0.47** |
| P6 `Schwann_like` | ATF3 0.31, HSF2 0.32, HIVEP2 0.32 | — |
| P2 neuronal | GATA3 0.59, HOXA9 0.54; MYC −0.64 | — |

Three of those are genuine cross-checks: the labels came from the programs' **own top genes**, while
the pathway layer is an entirely separate source, and it independently supports Hypoxia for P4, E2F/MYC
targets for P7 and EMT for P3.

**Broad CNA → expression consistency.** 56 broad clonal and 308 broad subclonal segments (≥10 Mb);
within each patient, carrier clones were compared against non-carrier clones of the **same** patient.
**299 events tested, 263 evaluable, direction matches the inferred event in 180 of 263 (68%)** —
largest effects MPNST_1's chr8 gains (Cohen's d 2.65–3.23). **68% is the honest number**: well above
chance, but a third of evaluable events do not move in the inferred direction. **And this is an
internal consistency check, not validation** — SCEVAN inferred those events from expression in the
first place. Segments containing NF1 and NF2 are tabulated with the permitted and prohibited wording
printed next to each row.

**Integrated architecture — Model C, selected against criteria specified before the evidence was
assembled.**

| quantity | value | verdict |
| --- | --- | --- |
| clones concentrated in one program | 11 of 18 (61%) | clone-constrained model needs ≥70% — **fails** |
| median η² | 0.059 | needs ≥0.30 — **fails** |
| between-clone / within-clone divergence | 0.0069 | needs ≥0.50 — **fails decisively** |
| clones spanning ≥2 programs | 4 of 18 (22%) | universally-plastic model needs ≥70% — **fails** |

> **Selected: Model C — mixed architecture.** A minority of programs associate detectably with
> specific CNA-defined clones while the great majority of program variance sits within clones. How
> much program diversity a clone contains is itself patient-specific (1.00 effective programs in
> MPNST_2, up to 2.83 in MPNST_1), so this is a **mixed** architecture rather than a single rule
> holding across the cohort — neither purely clone-constrained nor universally plastic.

**Robustness.** Median η² is **0.0586 at every clone-size threshold** (20/50/100 cells); median
effective programs per clone moves only 1.222 → 1.185 across a fourfold change in the
hard-assignment margin; restricting to High-confidence cells reproduces η² at Spearman 0.901, TF
associations at 0.940 and PROGENy at 0.962. Leave-one-patient-out was deliberately **not** applied to
the clone analyses, which are inherently within-patient. 77–81% of strong pooled TF/pathway
associations are concordant in ≥3 patients — **reported with the caveat** that at n = 4,
three-of-four sign agreement occurs ≈31% of the time by chance.

**A dependency defect worked around without moving a frozen package.** `decoupleR::get_collectri()`
and `OmnipathR::collectri()` both fail in the installed stack. **OmnipathR was not upgraded** — a
frozen Phase 3/4 dependency is not moved for a convenience wrapper. The identical data was taken
from OmniPath's documented REST endpoint and cached with provenance in `external/networks/`.

---

## 5. Integrated biological model

Synthesising Phases 2–6, and stated only as far as the evidence reaches:

> MPNST does not appear to be organised around a small set of universally recurrent malignant
> transcriptional states. Malignant phenotypes are **strongly patient-specific** and only **weakly
> constrained by CNA-defined clone structure**. CNA-informed analysis supports a **substantial
> mesenchymal/ECM-rich malignant compartment that conventional expression-based annotation
> misclassifies as fibroblast-like**. Despite that transcriptional heterogeneity, the tumours show
> **partial convergence at regulatory, pathway and tumour–microenvironmental levels.**

What actually converges:

* **Myeloid-directed communication.** A dense, reproducible predicted signal set (APP–CD74 across
  five receivers, CD99–PILRA, ANXA1–FPR1, HLA-F–LILRB1/2) whose published counterparts in other
  cancers are immunosuppressive — while independently the myeloid *transcriptional* program is best
  explained by **CSF1**, an axis already in clinical development for neurofibroma. These are carried
  as **separate evidence streams and never merged**.
* **Vascular and perivascular signalling.** Canonically angiogenic (VEGFA–KDR/FLT1/NRP1,
  ANGPTL4–CDH5, collagen–endothelial integrins) plus a **reciprocal perivascular Notch circuit**
  (tumour JAG1 → pericyte NOTCH3; pericyte/endothelial JAG1/JAG2/DLL4 → tumour NOTCH2) in a tumour
  type where Notch is already implicated in Schwann-cell transformation.
* **A mixed, not uniformly suppressive, lymphoid interface.** Class I presentation is retained and
  IL15 best explains the CD8/NK/CD4 programs, yet inhibitory HLA-E–KLRC1 and HLA-F–LILRB engagement
  is among the best-supported findings in the dataset. **This dataset cannot resolve which dominates
  functionally.**
* **Selected regulatory and pathway convergence.** E2F4/E2F1/MYC for the cycling program and
  HIF1A/ATF4/HSF1 for translation-stress recur across patients even though program *identity* does
  not; PROGENy and Hallmark independently support labels derived from the programs' own genes.
* **A reinterpretation of the ECM axis after refinement.** All four collagen/FN1–ITGAV_ITGB8 axes
  changed sender: what Phase 3 read as stroma-to-tumour ECM signalling is substantially
  **tumour-autocrine**, from malignant mesenchymal/ECM-like cells to malignant cells. This is the
  firmest and most directly testable Phase 4 result.

What does **not** converge: malignant state identity, malignant program identity, the exact malignant
fraction, and CNA architecture beyond a small shared broad-event core (chr18 loss, chr2/chr7 gain in
3/3 reliable patients, pairwise Jaccard only 0.20–0.23).

---

## 6. Major negative findings

These are scientifically important and are not hidden anywhere in this repository.

| Finding | Evidence |
| --- | --- |
| **No recurrent discrete malignant transcriptional states** | 0 of 8 recurrent at ≥3/4 patients; 97.2% of malignant cells in patient-private states (Phase 4 M31/M36) |
| **No recurrent continuous malignant programs at the pre-declared threshold** | 0 recurrent / 1 shared-limited (technical) / 7 patient-private; 0 recurrent at every K from 5–15 on a technical-gene-free universe (Phase 5) |
| **LochNESS did not support APP-context-associated receiver-state structure** | best descriptive p = 0.333 against a 0.167 p-floor; score unstable across embeddings and reference samples. It therefore did **not** independently corroborate APP–CD74 as a receiver-state driver (Phase 3 M24) |
| **NicheNet did not corroborate APP for macrophages** | APP is not among the top receiver ligands; the two streams nominate different leaders, bounding APP–CD74 to predicted engagement (Phase 3 M23) |
| **The refined malignant fraction is not patient-robust** | 32.63% → 0.214 when MPNST_4 is dropped, against a Phase 2 value of 0.215 (Phase 4 §11) |
| **MPNST_3 SCEVAN clone inference failed sanity checks** | primary↔sensitivity agreement 0.0864; immune-lineage "clones"; promotions disabled and clone structure excluded from Phases 4 and 6 |
| **A tumour-state → TME communication model could not be established** | 12 of 176 evidence rows reach ≥3 patients and all 12 belong to one 290-cell state that is itself 83% one patient — an evaluability artefact, not a biological preference (Phase 4 M33) |
| **CNA-defined clones explain only a small share of program variation** | median η² 0.059; ~94% of program variance within clones; between/within divergence ratio 0.0069 (Phase 6) |
| **Tumour ligands do not explain the fibroblast program** | NicheNet fibroblast receiver AUPR 0.022 (Phase 3 M23) |
| **Newly-supported interactions after refinement are not reproducible** | 91.8% rest on a single patient; only 27 of 5,421 reach ≥3 patients (Phase 4 M34) |
| **Clustering can over-partition continuous structure** | three separately-clustered MPNST_4 ECM states collapse onto one continuous cNMF program (Phase 5 M38) |

---

## 7. Computational methods

| Analysis | Method | Version / resource |
| --- | --- | --- |
| Language / runtime | R (conda `R_env`) | R 4.4.3 |
| Single-cell framework | Seurat / SeuratObject | 5.4.0 / 5.3.0 |
| Normalisation | SCTransform v2 (`sctransform`, `glmGamPoi`) | 0.4.3 / 1.18.0 |
| Doublet detection | scDblFinder | 2.20.2 |
| Marker testing | `FindAllMarkers` (wilcox, presto) after `PrepSCTFindMarkers` | presto 1.0.0 |
| Integration **embedding** | Harmony (on `sample_id`, package defaults) | harmony 1.2.4 |
| CCC consensus (primary LR) | LIANA | 0.1.14 (Consensus resource via OmnipathR 3.14.0) |
| Independent CCC pathway/network | CellChat | 2.2.0.9001 (CellChatDB.human: 3,233 interactions, 290 pathways) |
| Independent permutation LR | CellPhoneDB | 5.0.1 (DB v5.0.0, 1,000 permutations, threshold 0.10) |
| Receiver transcriptional response | NicheNet (`nichenetr`) | 2.2.1.1 (prior model, Zenodo 7074291) |
| Neighbourhood phenotype (orthogonal) | LochNESS, official MMCA formulation | MMCA commit `af629c49…` (vendored read-only) |
| CNA inference from expression | SCEVAN (+ yaGST) | 1.0.3 / 2017.8.25 |
| Continuous programs | cNMF | 1.7.1 (isolated `p5_cnmf_env`) |
| TF activity | decoupleR `run_ulm` + CollecTRI | 2.12.0 / 41,674 edges, 1,201 TFs |
| Pathway activity | PROGENy `run_mlm` | 14 footprint pathways (top-500 model) |
| Gene-set interpretation | MSigDB Hallmark / Reactome / GO-BP | v2024.1.Hs (GMTs with recorded md5s) |
| Workflow orchestration | Snakemake (**Phase 1 only**) | see `workflow/Snakefile` |

Only versions actually recorded in the repository's environment logs are listed. Where a package was
already installed, it was **not** upgraded — no dependency was moved to make a tool work, and no
Phase 2/3 package was downgraded.

---

## 8. Repository and code map

```text
Pilot_MPNST/
├── PROJECT.md              # authoritative project specification (per-phase scope and prohibitions)
├── PROGRESS.md             # milestone-by-milestone execution record
├── CHANGELOG.md            # technical changes and scientific decisions
├── README.md               # this document
├── config/
│   ├── config.yaml         # production (real-data) Phase 1 configuration
│   ├── config.test.yaml    # synthetic test configuration
│   ├── schemas/            # config JSON/YAML schema
│   └── phase2/             # annotation_map_M15.tsv, ccc_annotation_map.tsv (annotation as DATA)
├── workflow/
│   ├── Snakefile           # Phase 1 orchestration only
│   ├── profiles/slurm/     # SLURM execution profile
│   └── envs/               # conda specs + environment snapshots (R_env_portable.yaml)
├── scripts/
│   ├── R/                  # Phase 1 core R (extraction, QC, SCT, PCA, clustering, markers, provenance)
│   │   └── phase2/         # Phase 2 R (harmony, evaluation, clustering, markers, annotation, freeze)
│   ├── python/             # preflight checker, project-state validator, report/manifest generators
│   ├── phase3/{ccc,concordance,lochness,receiver_response}/
│   ├── phase4/{scevan,malignancy,tumor_states,ccc_refinement,figures,utils}/
│   ├── phase5/{programs,malignant_ecm,validation,utils}/
│   ├── phase6/{clone_program,plasticity,regulatory,pathways,cna_expression,validation,utils}/
│   └── shell/              # SLURM wrappers: phase1 at top level, then phase2/ … phase6/
├── tests/{unit,integration}/  # Phase 1 unit tests + clean-room synthetic integration runs
├── external/               # vendored/downloaded reference material (MMCA_ref, genesets, networks)
├── reports/                # all Markdown reports, milestone records and machine-readable TSVs
└── results/                # computed objects, figures and tables (largely Git-ignored)
```

### Where each phase's code and evidence live

| Phase | Purpose | Main scripts / entry points | Primary reports | Final object / deliverable |
| --- | --- | --- | --- | --- |
| **1** (M0–M9) | Independent preprocessing, QC, pre-integration baseline | `workflow/Snakefile`; `scripts/R/{extract_dataset,filter_and_detect_doublets,normalize_and_find_features,run_pca_and_evaluation,run_clustering_sweep,discover_markers,combine_pre_integration}.R`; SLURM `scripts/shell/run_m*_workflow.sh` | `reports/PHASE1_HANDOFF.md`, `reports/milestones/M0–M9` | `results/combined/pre_integration/combined_preintegration.rds`; `results/phase1_manifest.json` |
| **2** (M10–M17, M15A) | Harmony integration, clustering, markers, annotation, `annotation_ccc` | `scripts/R/phase2/{inspect_phase1_handoff,run_harmony_integration,evaluate_harmony,cluster_sweep_harmony,discover_markers_harmony,annotate_celltypes,refine_and_compose,build_ccc_annotation,validate_and_freeze}.R`; SLURM `scripts/shell/phase2/run_m1*.sh` | `reports/phase2/PHASE2_HANDOFF.md`, `{HARMONY,CLUSTERING}_ASSESSMENT.md`, `ANNOTATION_REPORT.md`, `CCC_READINESS.md`, `milestones/M10–M17` | `results/phase2/phase2_final_object.rds`; `results/phase2/phase2_manifest.json` |
| **3** (M18–M27) | Tumour–TME communication, concordance, receiver response, LochNESS | `scripts/phase3/ccc/{prepare_ccc_inputs,run_liana,run_cellchat,prioritize_interactions}.R` + `run_cellphonedb.py`; `scripts/phase3/concordance/{build_concordance,integrate_evidence,robustness,make_ccc_figures}.R`; `scripts/phase3/receiver_response/run_nichenet.R`; `scripts/phase3/lochness/{lochness_mpnst,compare_implementations}.R`; SLURM `scripts/shell/phase3/` | `reports/phase3/PHASE3_HANDOFF.md`, `PHASE3_METHOD_PLAN.md`, `CCC_CONCORDANCE_REPORT.md`, `MPNST_INTERACTION_REPORT.md`, `RECEIVER_RESPONSE_REPORT.md`, `LOCHNESS_MPNST_{DESIGN,REPORT}.md`, `ROBUSTNESS_REPORT.md`, `milestones/M18–M27` | `results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv` (primary deliverable); `results/phase3/phase3_manifest.json` |
| **4** (M28–M35, M35A) | SCEVAN malignancy refinement, clones, tumour states, CCC sensitivity | `scripts/phase4/scevan/{m28_feasibility,m29_run_scevan,m29_aggregate}.R`; `scripts/phase4/malignancy/{m30_integrate_malignancy,m34_robustness,m30_figures}.R`; `scripts/phase4/tumor_states/{m31_tumor_states,m33_state_tme_model}.R`; `scripts/phase4/ccc_refinement/`; `scripts/phase4/figures/m35a_*.R`; `scripts/phase4/utils/m35_build_final_object.R`; SLURM `scripts/shell/phase4/` | `reports/phase4/PHASE4_HANDOFF.md`, `MALIGNANCY_DECISION_RULES.md`, `MALIGNANCY_REFINEMENT_REPORT.md`, `SCEVAN_CNV_AND_CLONE_REPORT.md`, `AMBIGUOUS_POPULATION_MALIGNANCY_REPORT.md`, `TUMOR_STATE_COMMUNICATION_REPORT.md`, `M35A_SCEVAN_FIGURE_AUDIT.md`, `ROBUSTNESS_REPORT.md`, `milestones/M28–M35` | `results/phase4/phase4_final_object.rds`; `results/phase4/phase4_manifest.json` |
| **5** (M36–M41) | Continuous malignant programs (cNMF); malignant ECM vs fibroblast | `scripts/phase5/programs/{m37_prepare_input.R,m37_build_anndata.py,m37_k_selection.py,m38_annotate_programs.R}`; `scripts/phase5/malignant_ecm/m39_malignant_ecm.R`; `scripts/phase5/validation/m40_*.R`; `scripts/phase5/utils/m41_build_final_object.R`; SLURM `scripts/shell/phase5/` | `reports/phase5/PHASE5_HANDOFF.md`, `PHASE5_FEASIBILITY.md`, `ROBUSTNESS_REPORT.md`, `milestones/M36–M41` | `results/phase5/phase5_final_object.rds`; `results/phase5/phase5_manifest.json` |
| **6** (M42–M50) | Clone ↔ program coupling, within-clone diversity, TF/pathway architecture | `scripts/phase6/clone_program/m43_clone_program.R`; `scripts/phase6/plasticity/m44_within_clone_diversity.R`; `scripts/phase6/regulatory/m45_tf_activity.R`; `scripts/phase6/pathways/m46_pathway_activity.R`; `scripts/phase6/cna_expression/m47_cna_expression.R`; `scripts/phase6/utils/{m48_integrated_architecture,m50_build_final_object}.R`; `scripts/phase6/validation/m49_robustness.R`; SLURM `scripts/shell/phase6/` | `reports/phase6/PHASE6_HANDOFF.md`, `PHASE6_FEASIBILITY.md`, `ROBUSTNESS_REPORT.md`, `milestones/M42–M50` | `results/phase6/phase6_final_object.rds`; `results/phase6/phase6_manifest.json` |

---

## 9. Results, figures and tables

Computed artefacts live under `results/`, organised by phase; **`results/` is largely Git-ignored**
(the frozen objects are ~6 GB each) and only the small per-phase `*_manifest.json` documents and the
Phase 1 manifest are tracked. The `reports/` tree — every Markdown report and every machine-readable
TSV — **is** tracked.

| Phase | Final figures | Final tables | Where to start |
| --- | --- | --- | --- |
| 2 | `results/phase2/figures/final/` (+ per-milestone `M12`…`M16`) | `results/phase2/tables/final/` | `04_harmony_broad_celltypes`, `05_harmony_detailed_celltypes`, `08_celltype_composition_by_sample` |
| 3 | `results/phase3/figures/final/` (36 files; per-milestone in `results/phase3/figures/M18..M27/`) | `results/phase3/tables/final/` and `results/phase3/ccc/prioritized/` | **`MPNST_CCC_MASTER_TABLE.tsv`** — interaction → supporting methods → supporting patients → receiver response → LochNESS context → literature → priority tier |
| 4 | `results/phase4/figures/final/` — `01`–`16` core, `17`–`20` **native SCEVAN CNA heatmaps (verbatim)**, `21`–`28` M35A evidence figures | `results/phase4/tables/final/` | **`PHASE4_MALIGNANCY_CALLS.tsv`** (19,716 × 37, per-cell rule and printed reason), `SCEVAN_CLONES.tsv`, `SCEVAN_SAMPLE_RELIABILITY.tsv`, `MALIGNANCY_THRESHOLD_SENSITIVITY.tsv`, `SCEVAN_NATIVE_FIGURE_INDEX.tsv`, `PHASE3_VS_PHASE4_CCC.tsv` |
| 5 | `results/phase5/figures/final/` (12 figures, PDF + PNG) | `results/phase5/tables/final/` (16) | `MALIGNANT_ECM_SIGNATURE.tsv`, the program recurrence and threshold-sensitivity tables |
| 6 | `results/phase6/figures/final/` (13 figures, PDF + PNG) | `results/phase6/tables/final/` (12) | the clone × program η² table, TF/PROGENy/Hallmark association tables, broad CNA→expression table |

A consolidated cross-phase figure index is tracked at **`reports/FIGURE_INDEX.tsv`** (732 rows),
recording for every figure its phase, milestone, generating script, input object and checksum,
parameters, git commit and SLURM job ID.

---

## 10. Reproducibility and HPC execution

> [!CAUTION]
> **Never load, deserialize, inspect, subset or analyse the real objects (`processed_mpnst.rds` or
> anything under `results/`) on a login node.** All production computation runs inside a SLURM
> allocation. Each frozen phase object is ~6 GB.

**There is no single command that reproduces Phases 1–6.** Phase 1 is orchestrated by Snakemake;
**Phases 2–6 are driven by explicit per-milestone SLURM wrappers** under `scripts/shell/phase{2..6}/`
and are strictly sequential, each consuming the previous milestone's object. Per-phase reproduction
command sequences are given in each handoff (Phase 2 §14, Phase 3 §23, Phase 4 §24). **All seeds are
42 throughout the project.**

### Environment

```bash
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
```

`R_env` is the **frozen** analysis environment and was treated as immutable from Phase 2 onward: it
was snapshotted before and after each phase and asserted byte-identical for Phases 5 and 6. Where a
phase needed something the frozen stack could not provide, it went into an **isolated** environment
rather than replacing anything:

| Environment | Contents | Used by |
| --- | --- | --- |
| `R_env` | R 4.4.3 + Seurat/Harmony/SCEVAN/decoupleR stack | all phases |
| `cpdb_env` | Python 3.10, cellphonedb 5.0.1, anndata, scanpy | Phase 3 (and reused unmodified in Phase 4) |
| `p4_umap_env` | Python umap-learn 0.5.12 | Phase 4 only — SCEVAN's subclone stage requires it via reticulate |
| `p5_cnmf_env` | Python 3.11, cnmf 1.7.1, numpy/scipy/scikit-learn/anndata/scanpy | Phase 5 only |

The one dependency change recorded honestly: **igraph 2.2.1 → 2.1.4**, downgraded by the conda solver
during the Phase 3 install; Phase 2 outputs are frozen artefacts and unaffected. Portable spec:
`workflow/envs/R_env_portable.yaml`. Install logs: `reports/phase3/environment/PHASE3_INSTALL_LOG.md`,
`reports/phase4/environment/PHASE4_INSTALL_LOG.md`, `reports/audits/conda_env_diff.md`.

### Phase 1 checks and synthetic testing

```bash
# environment / config / permissions / dataset visibility
python scripts/python/preflight_checker.py --mode synthetic   # or --mode real
python scripts/python/validate_project_state.py

# synthetic smoke test — no real data touched
snakemake -n --configfile config/config.test.yaml
snakemake --cores 4 --configfile config/config.test.yaml

# production dry-run
snakemake -n --configfile config/config.yaml
```

CI (`.github/workflows/ci.yml`) is restricted to the **synthetic** configuration: preflight checks,
unit tests, the project-state validator and a Snakemake DAG regression check that fails if the
synthetic dry-run references real HPC paths.

---

## 11. Provenance

Every milestone writes a machine-readable record, and the repository is designed so that any number
in a report can be traced back to the run that produced it:

* **Per-milestone provenance JSON** — input paths with md5, outputs with md5, parameters, seeds, git
  commit, git status, `sessionInfo()` and SLURM JobID (via `scripts/R/provenance_utils.R`).
* **Per-phase manifests** — `results/phase{1..6}*manifest.json`, tracked in Git as documentation.
  Phase 6's has 31 sections; Phase 4's manifest writer additionally re-reads the file after writing
  and asserts every pre-existing section is unchanged, a guard added after `jsonlite`'s default
  `digits = 4` was caught silently rounding frozen values.
* **Object checksums** — md5 and sha256 for every phase object, re-verified before each downstream
  phase loaded it and after every milestone; the frozen upstream object was opened **read-only**.
* **Preservation guards** — each freeze asserts column-by-column that prior metadata is
  byte-identical, that reductions and cell order are unchanged, and that only expected columns were
  added (Phase 4: 23/23 checks; Phase 5: 10 guards + 8 reload validations; Phase 6: 8 + 8).
* **SLURM accounting** — every job, including **every failure**, is recorded with state, elapsed
  time, requested memory and peak RSS in the phase handoffs. Peak memory anywhere in the project was
  60.87 GiB (13.5% of the 450 G envelope), and **no failure was ever addressed by increasing memory
  or walltime** — each was root-caused and fixed at source.
* **Figure index** — `reports/FIGURE_INDEX.tsv`, one row per figure with its generating script,
  input checksum, git commit and job ID.

Exact hashes, per-figure records and full SLURM tables live in the phase manifests and handoffs
rather than in this README.

---

## 12. Limitations

Applying to the project as a whole:

1. **n = 4 patients.** No population-level, epidemiological or generalising claim follows from any
   result here.
2. **`sample_id` = sample = patient = dataset = presumed batch is one variable**, and no clinical or
   condition covariate exists. Biological and technical effects cannot be fully separated, and no
   condition-stratified, differential-abundance or pseudobulk-DE analysis was possible — none was
   fabricated.
3. **Strong between-patient heterogeneity.** Malignant-cell abundance is 17.4× unequal across
   patients and confounded with sequencing depth (the largest patient is the shallowest).
4. **CNA calls are expression-derived, not DNA-derived.** SCEVAN infers broad copy number from
   expression; single-gene CNV calls are not asserted, broad segments containing NF1/NF2 are
   described as such, and confirming NF1/NF2 status requires DNA sequencing.
5. **SCEVAN's annotation covers chr1–22 only**, and it removes cell-cycle and all `HLA-*` genes
   before inference — so the CNA analysis is structurally blind to the HLA-E/HLA-F loci that Phase 3
   highlighted.
6. **MPNST_3 failed the SCEVAN reliability gate.** Its malignant promotions were disabled and its
   clone structure is excluded, so Phase 6's clone conclusions rest on **three** patients, and any
   malignant mesenchymal population in MPNST_3 is undetectable by this analysis.
7. **Copy-number-quiet malignant cells are possible**, so a SCEVAN non-malignant call is not proof
   of non-malignancy — which is why 2,015 Phase 2 `MPNST-Tumor` cells became `Ambiguous` rather than
   `Non-malignant`. `Ambiguous` (3,766 cells) means NOT EVALUABLE, not "no signalling".
8. **The exact refined malignant fraction is not patient-robust** (32.63% cohort-wide; 21.4% without
   MPNST_4). The fibroblast-compartment conclusion is a separate, more robust claim in direction —
   though its magnitude depends on the post-hoc amendment A2, which is published alongside the flat
   gate it replaced.
9. **CCC results predict communication *potential***, not physical signalling: no adjacency, no
   directionality proof, no demonstrated receptor engagement. 53.6% of supported interactions rest on
   a single patient.
10. **The CCC frameworks are not fully independent** — LIANA's method set includes a
    CellPhoneDB-style score, and 97% of interaction keys exist in only one framework's resource, so
    "single-method" usually reflects resource non-overlap rather than disagreement.
11. **Malignant states and programs show limited recurrence** (0 of 8 in both Phase 4 and Phase 5),
    and a cNMF program is a continuous expression pattern — not a cell type, lineage, discrete state
    or clone.
12. **The Phase 5 malignant-ECM vs fibroblast comparison is evaluable in only 2 of 4 patients**, with
    a cross-patient log2FC correlation of 0.091; the 92-gene signature is descriptive, with no
    classifier and no held-out accuracy.
13. **Phase 6's broad CNA→expression analysis is an internal consistency check, not validation**,
    because the CNA events were themselves inferred from expression. Within-clone program diversity
    is consistent with plasticity but does **not** demonstrate a state transition — no rate,
    direction or trajectory is claimed.
14. **Only Harmony was evaluated** for integration; CCA, RPCA, MNN and scVI were not benchmarked, and
    no second CNV method (inferCNV, CopyKAT) was run.
15. **No spatial or experimental validation exists at this stage.** Every cross-compartment and
    communication statement in this project remains an inference from dissociated single-cell
    expression.

---

## 13. Phase 7 — planned spatial validation (not performed)

**Phase 7 has not been initiated and requires separate authorization.** It is described here as
planned work, and nothing in this repository implements it.

Its purpose is to test *in situ* what dissociated scRNA-seq cannot settle, prioritising the
predictions that are both firmest and most falsifiable:

1. **Malignant mesenchymal/ECM-like cells versus true fibroblast compartments** — the central Phase
   4/5 claim. Dual staining for a malignant marker together with COL1A1/FN1, and the Phase 5
   candidate signature (malignant side CA12/IGFBP3/COL11A1…; fibroblast side
   CDH19/APOD/SCN7A/ABCA6-10), would distinguish malignant ECM-like cells from genuine
   nerve-associated fibroblasts in tissue. No further scRNA-seq analysis of these four patients can
   do this.
2. **Tumour–myeloid architecture**, including APP/CD74 where panel coverage permits, and whether APP
   is polarised toward CD74+ myeloid contacts.
3. **VEGFA / vascular niches** — whether VEGFA is concentrated at the vascular front and KDR/FLT1
   receptors localise accordingly.
4. **JAG/NOTCH perivascular architecture** — the most spatially specific prediction in the project
   and the easiest to falsify, since it requires physical tumour–pericyte–endothelial contact.
5. **Selected lymphoid interfaces** — whether KLRC1+ NK cells and LILRB1/2+ monocytes are found where
   HLA-E and HLA-F are highest.
6. **Broad malignant programs rather than exact patient-private cNMF identities.** Given 0 recurrent
   programs at the pre-declared threshold, forcing patient-private program identities onto a spatial
   panel would be unsound; broad axes (ECM/mesenchymal, cycling, hypoxia/angiogenesis, Schwann-like)
   are the appropriate targets.

Complementary evidence that would convert inference into determination, and is equally not performed:
DNA sequencing or FISH on the recurrent chr18 loss and chr2/chr7 gains, NF1/NF2 status by DNA
sequencing, and **additional independent MPNST patients** — the only thing that can address the
recurrence questions, which more analysis of these four patients cannot.

---

## 14. Reference documentation

| Document | Contents |
| --- | --- |
| [`PROJECT.md`](PROJECT.md) | Authoritative project specification: per-phase scope, method mandates and explicit prohibitions |
| [`PROGRESS.md`](PROGRESS.md) | Milestone-by-milestone execution record, M0 → M50 |
| [`CHANGELOG.md`](CHANGELOG.md) | Technical changes and scientific decisions, including every amendment and post-hoc label |
| [`reports/PHASE1_HANDOFF.md`](reports/PHASE1_HANDOFF.md) | Phase 1 handoff (M0–M9) |
| [`reports/phase2/PHASE2_HANDOFF.md`](reports/phase2/PHASE2_HANDOFF.md) | Phase 2 handoff (M10–M17, M15A) |
| [`reports/phase3/PHASE3_HANDOFF.md`](reports/phase3/PHASE3_HANDOFF.md) | Phase 3 handoff (M18–M27) |
| [`reports/phase4/PHASE4_HANDOFF.md`](reports/phase4/PHASE4_HANDOFF.md) | Phase 4 handoff (M28–M35) + §28 M35A |
| [`reports/phase5/PHASE5_HANDOFF.md`](reports/phase5/PHASE5_HANDOFF.md) | Phase 5 handoff (M36–M41) |
| [`reports/phase6/PHASE6_HANDOFF.md`](reports/phase6/PHASE6_HANDOFF.md) | Phase 6 handoff (M42–M50) |
| `reports/{milestones,phase2..6/milestones}/` | Per-milestone reports, M0 → M50 |
| `reports/FIGURE_INDEX.tsv` | Cross-phase figure provenance index (732 rows) |
| `reports/audits/` | Reconciliation audits and the conda environment diff |

### Key methods citations

SCEVAN — De Falco, Caruso, Sudmant & Ceccarelli, *A variational algorithm to detect the clonal copy
number substructure of tumors from scRNA-seq data*, **Nat Commun 14:1074 (2023)**,
[doi:10.1038/s41467-023-36790-9](https://doi.org/10.1038/s41467-023-36790-9). MPNST biological and
annotation sources, with DOIs/PMIDs, are listed in
[`reports/phase2/ANNOTATION_REPORT.md`](reports/phase2/ANNOTATION_REPORT.md); prioritised interaction
literature with DOIs is in `results/phase3/tables/CCC_LITERATURE_EVIDENCE.tsv`.
