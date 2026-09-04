# Phase 4 — Malignancy Decision Rules

**Written 2026-09-03, BEFORE any SCEVAN classification result was inspected.** Every
threshold below was fixed a priori. This document exists so that no rule producing a
`malignancy_refined` call is hidden or subjective (Phase 4 §24).

---

## 1. What this integration is for

SCEVAN is given **substantial weight but not absolute authority** (§23). A single method's
call is not the deliverable; an *integrated, auditable* call is. Every cell receives:

| Field | Values |
| --- | --- |
| `malignancy_phase2` | `malignant` · `non-malignant` · `ambiguous` · `excluded` |
| `malignancy_scevan` | `malignant` · `non-malignant` · `not-assessed` |
| `malignancy_refined` | `Malignant` · `Non-malignant` · `Ambiguous` · `Excluded-low-quality` |
| `malignancy_confidence` | `High` · `Moderate` · `Low` |
| `malignancy_reason` | the rule ID and its printed justification |

`Ambiguous` is a **legitimate terminal outcome**. Ambiguous cells never become tumour (§29).

---

## 2. Evidence streams

| # | Stream | Source | Independent of Phase 2 labels? |
| --- | --- | --- | :--: |
| E1 | Phase 2 lineage annotation | frozen `annotation_ccc` | — |
| E2 | SCEVAN classification, **primary** | `pipelineCNA(norm_cell = NULL)` | **yes, fully** |
| E3 | SCEVAN classification, **sensitivity** | `norm_cell` = high-confidence immune, `FIXED_NORMAL_CELLS = FALSE` | partly |
| E4 | Large-scale CNA burden | SCEVAN CNA matrix, per cell | yes |
| E5 | Population consistency | SCEVAN malignant fraction within (`annotation_ccc` × `sample_id`) | derived |
| E6 | Lineage-marker support | module scores, canonical panels | — |
| E7 | Patient consistency | how many of the 4 patients agree at population level | derived |
| E8 | Clone membership | SCEVAN subclone assignment | yes |

### Derived quantities, defined precisely

* **`pop_frac`** — among cells that SCEVAN actually assessed (i.e. not `filtered`) in a given
  (`annotation_ccc`, `sample_id`) stratum, the fraction called `tumor` by the **primary** run.
  Strata with fewer than **20 assessed cells** are marked `pop_frac_unreliable` and treated as
  if `pop_frac` were unavailable; such cells cannot be promoted by R13 or R16a.
* **`runs_agree`** — the primary and sensitivity runs return the same call for that cell.
* **`cnv_burden`** — the fraction of assessed genomic segments in which that cell's inferred
  CNA value exceeds **0.10** in absolute value.
* **`cnv_elevated`** — `cnv_burden` exceeds the **95th percentile** of `cnv_burden` among
  SCEVAN's own **confident-normal** cells *in the same sample*. The reference is per-sample
  because `sample_id = patient = dataset` and baseline noise is not shared across patients.

---

## 3. Population classes

| Class | Populations | Role |
| --- | --- | --- |
| `EXCLUDED` | `Low-quality-excluded` | technical artefact (C15, 96% MPNST_4, mitochondrial-dominated) |
| `PHASE2_MALIGNANT` | `MPNST-Tumor` | Phase 2 called malignant on positive lineage evidence |
| `CANONICAL_IMMUNE` | `CD4-T` `CD8-T` `NK` `T-cell-other` `B-cell` `Plasma-cell` `Macrophage` `Monocyte` `Dendritic` `Plasmacytoid-DC` | strong, unambiguous lineage identity |
| `CANONICAL_ENDO` | `Endothelial` | strong identity, but a Phase 3 partner of interest — kept out of the reference set so its call stays independent |
| `DISPUTED_STROMAL` | `Fibroblast` `Pericyte-VSMC` | **the question**: may contain Mes-NC-like malignant cells |
| `PROVISIONAL` | `Candidate-Malignant-Unresolved` `Uncertain` | Phase 2 could not resolve these; excluded from Phase 3 |

`DISPUTED_STROMAL` and `PROVISIONAL` are **never** used as fixed normal references (§16).

---

## 4. The rules — ordered, first match wins

| ID | Condition | `malignancy_refined` | Confidence | Why |
| --- | --- | --- | --- | --- |
| **R0** | class `EXCLUDED` | `Excluded-low-quality` | High | Mitochondrial-dominated technical artefact. Expression is unreliable, therefore inferred CNV is unreliable. Not promoted whatever SCEVAN says. |
| **R1a** | SCEVAN `filtered` **and** class `CANONICAL_IMMUNE` | `Non-malignant` | Moderate | No CNV evidence, but canonical immune identity is strong marker evidence on its own. |
| **R1b** | SCEVAN `filtered`, any other class | `Ambiguous` | Low | No CNV evidence and no decisive lineage evidence. Honest non-answer. |
| **R2** | `PHASE2_MALIGNANT` · SCEVAN malignant · `runs_agree` | **`Malignant`** | High | Both evidence streams agree, under both reference strategies. |
| **R3** | `PHASE2_MALIGNANT` · SCEVAN malignant · runs disagree | `Malignant` | Moderate | Agreement on the call, sensitivity to the reference choice. |
| **R4** | `PHASE2_MALIGNANT` · SCEVAN non-malignant · `pop_frac` ≥ 0.25 | `Malignant` | Moderate | The population is substantially malignant; this cell is plausibly copy-number quiet (§65C). |
| **R5** | `PHASE2_MALIGNANT` · SCEVAN non-malignant · `pop_frac` < 0.25 | `Ambiguous` | Low | Genuine disagreement, flagged for inspection rather than resolved by fiat (§25). |
| **R6** | `CANONICAL_IMMUNE` · SCEVAN non-malignant | `Non-malignant` | High | Concordant. |
| **R7** | `CANONICAL_IMMUNE` · SCEVAN malignant · `pop_frac` < 0.50 | `Non-malignant` | Moderate | A minority CNA-cluster call does not overturn canonical immune identity. Recorded as a SCEVAN disagreement. |
| **R8** | `CANONICAL_IMMUNE` · SCEVAN malignant · `pop_frac` ≥ 0.50 | `Ambiguous` | Low | A majority-malignant immune population is more likely a SCEVAN artefact than a real finding, but it is **flagged, not silently overridden**. |
| **R9** | `CANONICAL_ENDO` · SCEVAN non-malignant | `Non-malignant` | High | Concordant. |
| **R10** | `CANONICAL_ENDO` · SCEVAN malignant · `pop_frac` < 0.50 | `Non-malignant` | Moderate | As R7. |
| **R11** | `CANONICAL_ENDO` · SCEVAN malignant · `pop_frac` ≥ 0.50 | `Ambiguous` | Low | As R8. |
| **R12** | `DISPUTED_STROMAL` · SCEVAN non-malignant · `pop_frac` < 0.25 | `Non-malignant` | High | **Supports a true fibroblast / pericyte** — the population is predominantly CNV-normal. |
| **R13** | `DISPUTED_STROMAL` · SCEVAN malignant · `runs_agree` · `pop_frac` ≥ 0.50 · `cnv_elevated` | **`Malignant`** | High | **The key promotion rule.** A majority-malignant stromal-like stratum, agreeing under both reference strategies, with CNA burden above the confident-normal 95th percentile → candidate Mes-NC-like malignant population. |
| **R14** | `DISPUTED_STROMAL` · SCEVAN malignant · (`pop_frac` ≥ 0.25) | `Malignant` | Moderate | Promoted, but one of `runs_agree` / `pop_frac ≥ 0.50` / `cnv_elevated` is missing. |
| **R15** | `DISPUTED_STROMAL` · SCEVAN malignant · `pop_frac` < 0.25 | `Ambiguous` | Low | Isolated malignant calls inside a predominantly normal stromal population. Not promoted. |
| **R16** | `DISPUTED_STROMAL` · SCEVAN non-malignant · `pop_frac` ≥ 0.25 | `Non-malignant` | Moderate | Called normal, but inside a stratum with meaningful malignant signal — confidence reduced. |
| **R17** | `PROVISIONAL` · SCEVAN malignant · `runs_agree` · `pop_frac` ≥ 0.50 | **`Malignant`** | High | Resolves what Phase 2 explicitly left open, on evidence Phase 2 did not have. |
| **R18** | `PROVISIONAL` · SCEVAN malignant | `Malignant` | Moderate | Promoted on weaker support. |
| **R19** | `PROVISIONAL` · SCEVAN non-malignant · `pop_frac` < 0.25 | `Non-malignant` | Moderate | Never `High`: Phase 2 could not establish these cells' identity, so CNV evidence alone earns `Moderate`. |
| **R20** | `PROVISIONAL` · SCEVAN non-malignant · `pop_frac` ≥ 0.25 | `Ambiguous` | Low | Mixed stratum, unresolved population. Stays unresolved. |
| **R99** | anything unmatched | `Ambiguous` | Low | Fail-safe. A cell that reaches R99 is a rule-coverage bug and is reported as such. |

### Marker evidence as a confidence modifier (E6)

Marker evidence does not create a call; it can only **reduce confidence by one level**. After
the table above, for any cell called `Malignant` at `High`, if its population's malignant-panel
module score is *not* above the object-wide median **and** its canonical non-malignant panel
score *is*, confidence drops to `Moderate` and the reason string records
`marker_evidence_discordant`. The reverse downgrade applies to `Non-malignant` / `High`.

This direction is deliberate: markers are what Phase 2 already used, so letting them *raise*
confidence in a Phase 4 call would smuggle circularity back in (§16).

---

## 5. What these rules deliberately do NOT do

* They never execute `if (!immune) malignant` in any form. The promotion of a stromal-like cell
  requires positive CNA evidence — rules R13/R14 — and cannot be reached by exclusion.
* They never force Phase 2 and SCEVAN to agree. R5, R8, R11, R15 and R20 exist precisely to
  keep disagreement visible (§25).
* They never treat a SCEVAN `non-malignant` call as proof of non-malignancy: R4 exists because
  some MPNST malignant cells may be copy-number quiet (§65C), and no cell is ever assigned
  `Non-malignant` at `High` confidence on a SCEVAN call alone — `High` always requires
  concordant lineage evidence as well (R6, R9, R12).
* They never promote `Ambiguous` cells to tumour, in this milestone or in the refined CCC
  annotation (§29).

---

## 6. Refined CCC annotation (`annotation_ccc_refined`)

Built from `malignancy_refined`, preserving `annotation_ccc_phase3` verbatim alongside it.

```text
malignancy_refined == "Malignant"                  -> "MPNST-Tumor"
malignancy_refined == "Non-malignant"              -> the Phase 3 annotation_ccc label, unchanged
malignancy_refined == "Ambiguous"                  -> "Ambiguous-unresolved"      (NOT tumour)
malignancy_refined == "Excluded-low-quality"       -> "Low-quality-excluded"      (unchanged)
```

Tumour-state labels are **not** erased by this collapse: `tumor_state_phase4` and
`tumor_clone_phase4` carry the heterogeneity within `MPNST-Tumor` (§29, and Phase 2 §P2.4's
"do not delete tumour heterogeneity" constraint, which remains in force).

---

# Amendment A1 — per-sample SCEVAN immune sanity gate

**Added 2026-09-03, AFTER M29 completed. This was NOT fixed a priori, and saying otherwise
would be dishonest.** It was added in response to an observed, specific failure.

## What was observed

In **MPNST_3**, the two SCEVAN runs disagreed almost completely — primary/sensitivity call
agreement **0.0864** over 2,940 cells — and *both* runs assigned large blocks of unambiguous
immune cells to the malignant class:

| Population | primary → tumour | sensitivity → tumour | of assessed |
| --- | ---: | ---: | ---: |
| B-cell | **89** | 0 | 89 |
| CD4-T | **327** | 4 | 328 |
| CD8-T | **135** | 0 | 136 |
| NK | **32** | 0 | 32 |
| Plasma-cell | **415** | 2 | 418 |
| Plasmacytoid-DC | **165** | 0 | 165 |
| Macrophage | 4 | **678** | 680 |
| Dendritic | 8 | **142** | 144 |
| Monocyte | 1 | **72** | 72 |

The primary run calls the lymphoid compartment malignant and the myeloid compartment normal;
the sensitivity run does the exact opposite. Meanwhile only **91 of 212** Phase 2 `MPNST-Tumor`
cells were called malignant by the primary run.

## Why it happens

SCEVAN found only **25 confident normal cells** in MPNST_3, and that sample has the lowest
depth of the four (median 4,265 counts and **1,594 genes** per cell, versus 3,139 in MPNST_1).
With the baseline estimated from 25 cells, `classifyTumorCells()` → `classifyCluster()`
partitions the top-level CNA hierarchical clustering along whichever axis dominates the
residual expression structure — here **lymphoid versus myeloid**, not malignant versus normal.
The two runs then label those two clusters oppositely depending on which reference set they
started from.

This is a known limitation of expression-derived CNV inference, not a coding error: SCEVAN
assumes a detectable confident-normal population and a substantial malignant population, and
MPNST_3 satisfies neither well (Phase 2 called 231 of 2,940 cells malignant, 7.9%).

## The gate

For each sample and each run, compute

```text
immune_malignant_rate = (canonical immune cells called malignant) / (canonical immune cells assessed)
```

A run **FAILS the immune sanity check** when `immune_malignant_rate > 0.25`. A run that calls
more than a quarter of unambiguous T, NK, B, plasma, macrophage, monocyte and dendritic cells
malignant is not producing a usable malignant/normal partition, whatever its internal
statistics say.

`scevan_sample_reliable` is TRUE only when the **primary** run passes. The primary run is used
because it is the non-circular one; **the passing run is never selected per sample to suit the
answer**, which would be cherry-picking.

## Consequences for the rules

| Rule | Change when `scevan_sample_reliable` is FALSE |
| --- | --- |
| **R13, R14** (`DISPUTED_STROMAL` → Malignant) | **disabled.** New rule **R13x** → `Ambiguous` / `Low` |
| **R17, R18** (`PROVISIONAL` → Malignant) | **disabled.** New rule **R17x** → `Ambiguous` / `Low` |
| **R2** (`PHASE2_MALIGNANT` + SCEVAN agree) | new rule **R2u** caps confidence at `Moderate` |
| R6–R12, R16, R19, R20 (`Non-malignant` outcomes) | unchanged |

The asymmetry is deliberate and conservative. An unreliable CNV run must not be allowed to
**promote** a cell into the malignant compartment — that is the inference Phase 4 exists to
make carefully. It may still contribute to a `Non-malignant` call, because those calls always
require concordant lineage evidence as well (§4, rules R6/R9/R12), and because Phase 2's
marker-based annotation remains valid evidence in its own right.

Cells affected by R13x/R17x carry `sample_scevan_reliable=FALSE` in `malignancy_reason`, so
every downgrade is traceable to this amendment rather than buried.

## What this costs, stated plainly

Any Mes-NC-like malignant fibroblast population that exists **in MPNST_3** cannot be detected
by this analysis. That is a real loss of sensitivity in one of four patients, and the
fibroblast conclusion in §69 must therefore be read as resting on the samples that passed.
`results/phase4/tables/SCEVAN_SAMPLE_RELIABILITY.tsv` records the rate for every sample and
both runs so a reader can see exactly which patients contribute.

---

# Amendment A2 — refining the sanity gate, and the plasma-cell artefact

**Added 2026-09-03, after M29 and after Amendment A1 had already been applied once.
The ordering matters and is stated plainly, because the revision made my strongest
result usable and that is exactly the circumstance in which a reader should be
suspicious.**

## What happened, in order

1. Amendment A1 introduced a flat gate: a sample fails if more than 25% of its canonical
   immune cells are called malignant.
2. Applied to the completed M29 runs, it failed **two** samples:

   | Sample | rate over all canonical immune | verdict under A1 |
   | --- | ---: | --- |
   | MPNST_1 | 0.036 | PASS |
   | MPNST_2 | 0.055 | PASS |
   | MPNST_3 | **0.590** | FAIL |
   | MPNST_4 | **0.293** | FAIL |

3. Failing MPNST_4 disabled promotions there, sending **2,967 cells to `Ambiguous` via R13x**
   — including the 2,887 MPNST_4 fibroblasts that SCEVAN called malignant, the single largest
   piece of evidence in the analysis.
4. That prompted me to ask **why** MPNST_4 failed. The two failures turn out not to be the same
   phenomenon at all.

## The two failures are different

**MPNST_3 is a global inversion.** The primary run calls the lymphoid compartment malignant
and the myeloid compartment normal; the sensitivity run does the reverse; agreement is 0.0864.
**Six of nine** testable immune populations individually exceed 25%: B-cell 1.00, CD4-T 1.00,
CD8-T 0.99, NK 1.00, Plasmacytoid-DC 1.00, T-cell-other 1.00. The malignant/normal partition
is structurally wrong.

**MPNST_4 is one artefact-prone population.** Its failure is driven almost entirely by
**Plasma-cell: 608 of 619 assessed (98.2%)**. Every other immune population is clean —
B-cell 0.036, CD4-T 0.007, CD8-T 0.000, Dendritic 0.000, Macrophage 0.000, Monocyte 0.000,
NK 0.024, Plasmacytoid-DC 0.000, T-cell-other 0.031. **Zero of nine** populations fail. And
primary/sensitivity agreement is **0.9974** — the highest of the four samples.

## Why plasma cells, mechanistically

Plasma cells express immunoglobulin loci — **IGH at 14q32, IGK at 2p11, IGL at 22q11** — at
extraordinary levels. Expression-derived CNV inference reads that as large segmental gain at
those loci. This is a **documented artefact class** for inferCNV-, CopyKAT- and SCEVAN-style
methods, it follows from how the methods work rather than from anything about this dataset, and
plasma cells are not part of the malignancy question Phase 4 exists to answer. The same
signature appears at lower magnitude in the passing samples too (MPNST_1 plasma 0.197,
MPNST_2 plasma 0.635), which is itself evidence that this is a population property rather than
a per-sample failure.

## The revised gate

A sample's SCEVAN run **fails** if either:

* **rate criterion** — more than **25%** of its canonical immune cells *excluding Plasma-cell*
  are called malignant; **or**
* **breadth criterion** — **3 or more** canonical immune populations (excluding Plasma-cell,
  and requiring ≥20 assessed cells to be judged at all) *individually* exceed 25%.

The breadth criterion exists so that a global inversion cannot pass by diluting itself, and it
is what independently condemns MPNST_3.

| Sample | rate excl. plasma | populations failing | A1 flat gate | **A2 verdict** |
| --- | ---: | ---: | --- | --- |
| MPNST_1 | 0.000 | 0 / 8 | PASS | **PASS** |
| MPNST_2 | 0.000 | 0 / 5 | PASS | **PASS** |
| MPNST_3 | **0.494** | **6 / 9** | FAIL | **FAIL (rate and breadth)** |
| MPNST_4 | 0.004 | 0 / 9 | FAIL | **PASS** |

`SCEVAN_SAMPLE_RELIABILITY.tsv` records **both** rates, the per-population failures, and a
`flat_gate_would_have_failed` column, so a reader who prefers the stricter A1 gate can apply it
and see exactly what changes.

## New rule R7p

| ID | Condition | `malignancy_refined` | Confidence | Why |
| --- | --- | --- | --- | --- |
| **R7p** | `annotation_ccc == "Plasma-cell"` and SCEVAN malignant | `Non-malignant` | Moderate | Attributed to the immunoglobulin-locus artefact above. `Moderate`, not `High`, because the CNV evidence genuinely disagrees. |

R7p is placed before R6/R7/R8. Without it those 1,709 plasma cells would land in `Ambiguous`
and be dropped from the refined CCC analysis — losing a Phase 3 population for a reason we can
actually name.

## What a sceptical reader should check

The revision helps the fibroblast conclusion, so it must not be the only thing holding it up.
Two independent checks are reported:

1. **The finding replicates in MPNST_1 and MPNST_2**, which passed under **both** gates:
   fibroblasts are 62.8% and 55.9% malignant there. MPNST_4 strengthens an already-replicated
   result rather than creating it.
2. **M34 reports leave-one-patient-out explicitly** (`MALIGNANCY_LEAVE_ONE_PATIENT_OUT.tsv`,
   `MALIGNANCY_PATIENT_CONSISTENCY.tsv`), so the effect of dropping MPNST_4 is quantified
   rather than argued.

If the fibroblast conclusion held only in MPNST_4, it would rest entirely on this amendment,
and the report says so.
