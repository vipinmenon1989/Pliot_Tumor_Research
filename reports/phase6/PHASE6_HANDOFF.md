# MPNST Phase 6 Handoff — Clonal, Regulatory and Transcriptional-Plasticity Architecture

*Milestones M42–M50 · executed and frozen 2026-09-04 · readable without conversation history*

> **MPNST_3's SCEVAN clone structure failed the Phase 4 immune sanity gate and is EXCLUDED from
> every clone-based inference.** Clone conclusions rest on three patients.
> **Clone labels are patient-scoped.** `MPNST_1_clone1` and `MPNST_4_clone1` are unrelated names and
> were never compared as homologous.
> **SCEVAN infers copy number from expression.** The CNA→expression analysis is an internal
> consistency check, **not** independent validation.
> **Within-clone program diversity is not observed state switching.** No transition rate, direction
> or trajectory is claimed anywhere.
> **Phases 1–5 are unmodified**; both upstream md5s were re-verified after every milestone.

---

## 1. The questions

> Is malignant transcriptional heterogeneity primarily associated with distinct CNA-defined clones,
> with substantial program diversity *within* clones, or a mixture?

> Which transcription-factor and pathway activities distinguish the Phase 5 malignant programs?

## 2. Input

`results/phase5/phase5_final_object.rds`, verified against its own manifest by md5
`839e5157bc7c3470ddf86746c2e719e1` and sha256 `fe99ecf5…`. The frozen Phase 4 refined counts were
re-asserted inside it and hold. The MPNST_3 exclusion was **re-derived** from the frozen
`scevan_sample_reliable` flag rather than taken on trust.

19 clones across MPNST_1/2/4; **18 evaluable** at the declared 20-cell minimum
(`MPNST_4_clone8` = 5 cells, reported NOT EVALUABLE). 178 malignant cells carry no clone label and
are excluded from clone analyses with the number stated.

## 3. Clone ↔ program coupling — the central result

Every association is computed **inside one patient**, with a permutation null built by shuffling
clone labels within that patient (1,000 permutations), which preserves clone sizes and the score
distribution.

| patient | clones | median η² | max η² | clone-associated programs |
| --- | ---: | ---: | ---: | ---: |
| MPNST_1 | 7 | 0.072 | 0.508 | 3 of 8 |
| MPNST_2 | 4 | 0.050 | 0.154 | 2 of 8 |
| MPNST_4 | 7 | 0.072 | 0.208 | 2 of 8 |

**Median η² ≈ 0.06 — roughly 94% of each program's variance sits WITHIN clones.** 7 of 24
program × patient pairs are clone-associated. 21 of 24 reach permutation p < 0.001, which is exactly
why effect size and not the p-value carries the conclusion: with 651–3,623 cells almost any non-zero
η² is "significant".

**A patient's dominant program is not a property of one clone.** MPNST_1's P2 is active in 7 of 7
clones (median usage 0.26–0.81), MPNST_2's P3 in 4 of 4 (0.84–0.93), MPNST_4's P1 in 7 of 7
(0.85–0.89).

## 4. Within-clone diversity — real, and strongly patient-specific

Measured on both continuous scores and a hard assignment whose threshold is printed (a cell keeps
its top program only when it leads the runner-up by ≥ 0.10 relative usage).

| patient | effective programs per clone | dominant-program share | Mixed cells |
| --- | --- | --- | --- |
| MPNST_1 | **1.69 – 2.83** | 0.43 – 0.78 | 3 – 25% |
| MPNST_2 | **1.00** (all four clones) | **1.000** | 0% |
| MPNST_4 | 1.05 – 1.23 | 0.95 – 0.99 | 0 – 1% |

**Between-clone divergence against within-clone spread** — the quantity that decides the
architecture:

| patient | between-clone JSD | within-clone dispersion | ratio |
| --- | ---: | ---: | ---: |
| MPNST_1 | 0.0130 | 0.384 | 0.034 |
| MPNST_2 | 0.00018 | 0.129 | 0.0014 |
| MPNST_4 | 0.00085 | 0.123 | 0.0069 |

**A patient's CNA-defined clones are transcriptionally near-interchangeable.**

**Confounders, assessed before any plasticity wording.** Program dispersion tracks the fraction of
High-confidence malignant cells in a clone at **ρ = 0.72**, above the declared |ρ| ≥ 0.7 disqualifying
bar. The two effective-number metrics stay below it (0.56–0.62), so the multi-program finding does
not rest on the confounded metric — but any plasticity reading of *dispersion* must carry this
number. Clone size (ρ ≤ 0.17) and sequencing depth (ρ ≤ 0.23 for nCount) are not drivers;
nFeature_RNA reaches 0.53 and the cell-cycle program 0.51.

**Permitted statement**: four clones, all in MPNST_1, **contain cells spanning multiple malignant
transcriptional programs**. That is consistent with transcriptional plasticity. It is **not**
observed switching.

## 5. Regulatory architecture

decoupleR `run_ulm` over **CollecTRI** (41,674 edges, 1,201 TFs) on the frozen RNA log-normalised
layer. `*` = direction concordant in ≥ 3 evaluable patients **and** |pooled ρ| ≥ 0.20.

* **P7 Cycling → E2F4 0.43\*, MYC 0.34\*, E2F1 0.24\***
* **P5 Translation → HIF1A 0.45\*, ATF4 0.43\*, HSF1 0.45\*, CREB1 0.43\***
* P3 Mesenchymal_ECM → HMGA2 0.57\*, MYC 0.56\*, SP1 0.56\*, STAT6 0.54\*
* P6 Schwann_like → ATF3 0.31\*, HSF2 0.32\*, HIVEP2 0.32\*
* P2 Neuronal → GATA3 0.59\*, HOXA9 0.54\*; MYC −0.64\*
* P4 Hypoxia_Angio → ID4 0.45\*, PGR 0.44\*, TBX2 0.43\*

The phase brief listed AP-1, TEAD/YAP, STAT/IRF, E2F, MYC, NF-kB and SOX-family as **hypotheses
only**. E2F and MYC were recovered for the cycling program and HIF1A/ATF4/HSF1 for the
translation-stress program **without being imposed**; TEAD/YAP and SOX-family did not emerge and are
not claimed.

**A dependency defect worked around without moving a frozen package.**
`decoupleR::get_collectri()` and `OmnipathR::collectri()` both fail in the installed stack
(OmnipathR 3.14.0's `unnest_evidences()` errors on the CollecTRI static table). **OmnipathR was not
upgraded** — a frozen Phase 3/4 dependency is not moved for a convenience wrapper. The identical
data was taken from OmniPath's documented REST endpoint and cached with its provenance in
`external/networks/`.

## 6. Pathway architecture — independent support for the program labels

PROGENy (14 footprint pathways, `run_mlm`) and MSigDB Hallmark (50 sets, mean z), **reported side by
side and never merged into a composite score or a single numeric rank**.

Three cross-checks matter, because the labels came from the programs' own top genes and the pathway
layer is an entirely separate source:

* **P4 was labelled `Hypoxia_Angio` from its genes; the PROGENy Hypoxia footprint is its top pathway
  (0.47\*).**
* **P7 was labelled `Cycling`; Hallmark E2F_TARGETS (0.43\*) and MYC_TARGETS_V1 (0.44\*) are its top
  sets, matching M45's E2F4/E2F1/MYC.**
* **P3 was labelled `Mesenchymal_ECM`; Hallmark EMT (0.33\*) is among its top sets, with EGFR
  (0.57\*) leading PROGENy.**

P5's MTORC1 (0.58\*) / GLYCOLYSIS (0.55\*) / PI3K-AKT (0.50\*) signature shows that the
ribosomal-pseudogene factor is a coherent translation-and-growth axis — which does not make it a
distinct malignant identity, and it remains flagged technical-dominated.

## 7. Broad CNA → transcriptional consequences

56 broad clonal and 308 broad subclonal segments (≥ 10 Mb) in the reliable patients; gene
coordinates from the stored SCEVAN annotation. Within each patient, carrier clones were compared
with non-carrier clones of the **same** patient on the mean expression of all genes in the interval.

**299 events tested, 263 evaluable, direction matches the inferred event in 180 of 263 (68%).**
Largest effects are MPNST_1's chr8 gains (Cohen's d 3.23, 3.09, 2.97, 2.65 across four adjacent
segments), with chr19 and chr11 losses at d ≈ −1.6 and −1.3.

**68% is the honest number**: well above chance, confirming that the clone assignments carry real
chromosome-scale expression structure, but a third of evaluable events do not move in the inferred
direction. **And this is not validation** — SCEVAN inferred those events from expression in the
first place.

Segments containing NF1 and NF2 are tabulated with the permitted and prohibited wording printed next
to each row. Permitted: *"a broad inferred loss segment whose interval contains the NF1 locus"*.
Prohibited and not used anywhere: *"NF1 deletion"*, *"NF1-deleted cells"*.

## 8. Integrated architecture — Model C

Criteria were specified before the evidence was assembled.

| quantity | value | verdict |
| --- | --- | --- |
| clones concentrated in one program | 11 of 18 (61%) | A needs ≥ 70% — **fails** |
| median η² | 0.059 | A needs ≥ 0.30 — **fails** |
| between-clone / within-clone | 0.0069 | A needs ≥ 0.50 — **fails decisively** |
| clones spanning ≥ 2 programs | 4 of 18 (22%) | B needs ≥ 70% — **fails** |

> **SELECTED — Model C, mixed architecture.** A minority of transcriptional programs associate
> detectably with specific CNA-defined clones while the great majority of program variance sits
> within clones. Median η² is 0.059, so roughly 94% of each program's variance is within-clone, and
> between-clone divergence is only 0.7% of within-clone dispersion — a patient's CNA-defined clones
> are transcriptionally near-interchangeable. How much program diversity a clone contains is itself
> patient-specific, so this is a mixed architecture rather than a single rule holding across the
> cohort.

## 9. Per-patient tumour portraits

| | MPNST_1 | MPNST_2 | MPNST_4 | MPNST_3 |
| --- | --- | --- | --- | --- |
| clone structure usable | yes | yes | yes | **NO — excluded** |
| malignant cells | 1,886 | 651 | 3,685 | 212 |
| clones / evaluable | 7 / 7 | 4 / 4 | 8 / 7 | 3 / — |
| dominant program | P2 Neuronal | P3 **Mesenchymal_ECM** | P1 Translation_ribosomal | P6 Schwann_like |
| effective programs per clone | **2.15** | 1.00 | 1.19 | — |
| median η² | 0.072 | 0.050 | 0.072 | — |
| clone-associated programs | 3 | 2 | 2 | — |

**MPNST_3**: Phase 5 programs described, **no clone-based interpretation offered or invented**.

## 10. Robustness

Full detail in `ROBUSTNESS_REPORT.md`. Median η² is **0.0586 at every clone-size threshold**
(20 / 50 / 100 cells); median effective programs per clone moves only 1.222 → 1.185 across a
fourfold change in the hard-assignment margin; and restricting to High-confidence cells reproduces
η² at Spearman **0.901**, TF associations at **0.940** and PROGENy at **0.962**. Leave-one-patient-out
was deliberately not applied to the clone analyses, which are inherently within-patient.

77–81% of strong pooled TF/pathway associations are concordant in ≥ 3 patients — **a number reported
with its caveat**, since with n = 4 three-of-four sign agreement occurs ≈ 31% of the time by chance.

## 11. Final object

```text
path      results/phase6/phase6_final_object.rds
size      6,060,757,871 bytes
md5       b8c01dd01755dda10b2be4e0f5ef7ef7
sha256    bde592d461d22253b646da1426db8545404b9d31b7f1b948e23deed03aa643f6
cells     19,716        metadata  196 Phase 1-5 columns + 46 Phase 6 columns
```

Added: 14 `progeny_*` pathway activities · 25 `tfact_*` TF activities (the strongest
cross-patient-concordant regulons) · `clone_hard_effective_n`, `clone_continuous_effective_n`,
`clone_dominant_program_share_hard`, `clone_program_dispersion`, `clone_dominant_program_label`,
`clone_evaluable_for_diversity`, `clone_structure_reliable_phase6`. The **full** 690-regulon TF
matrix and 50-set Hallmark matrix stay as separate artefacts referenced by the manifest rather than
bloating the object.

**8 preservation guards and 8 reload validations passed**, asserted column by column, including that
every Phase 5 `program_P*_score` is byte-identical and the frozen Phase 4 counts still hold.

## 12. Limitations

**A.** n = 4 patients; `sample_id` = patient = dataset. **B.** MPNST_3's clone structure is unreliable
and excluded, so clone conclusions rest on **3** patients. **C.** SCEVAN infers broad copy number
from expression — not DNA sequencing. **D.** Broad CNA–expression association is an **internal
consistency** analysis, not orthogonal validation, because the CNA call was itself derived from
expression. **E.** cNMF programs are continuous transcriptional patterns, not discrete states, cell
types, lineages or clones. **F.** Within-clone program diversity is consistent with plasticity but
does not demonstrate a transition; no rate, direction or trajectory is claimed. **G.** Cell-level
p-values do not represent biological replication. **H.** Clone labels are patient-scoped and were
never compared across patients. **I.** No single-gene CNV claim is made; loci are described only as
contained within a broad inferred segment. **J.** Program dispersion is confounded with the
fraction of High-confidence cells in a clone (ρ = 0.72). **K.** True generalization requires
independent MPNST patients.

## 13. SLURM accounting

| JobID | stage | State | Elapsed | CPUs | ReqMem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| 19899746 | Phase 6 M42–M50, first attempt | **FAILED** (M43) | 00:02:14 | 8 | 128 G | 22.99 GiB |
| 19899748 | **M42 · M43 · M44** | **FAILED** (M45) | 00:02:42 | 8 | 128 G | 21.92 GiB |
| 19899749 | **M45 – M50** | COMPLETED | 00:12:11 | 8 | 128 G | **8.21 GiB** |
| 19899752 | figure re-render + re-finalize | COMPLETED | 00:01:00 | 2 | 32 G | 0.19 GiB |

**Peak memory 22.99 GiB against the 450 G envelope — 5.1%.** No failure was addressed by increasing
memory or walltime.

**Two failures, both root-caused and fixed at source, both preserved above.**
*19899746* — `m43_clone_program.R` grouped by a column `clone` that only exists after renaming
`tumor_clone_phase4`. Fixed by renaming explicitly.
*19899748* — `decoupleR::get_collectri()` errored inside OmnipathR 3.14.0's `unnest_evidences()`.
**OmnipathR was not upgraded**; CollecTRI was fetched from OmniPath's documented REST endpoint and
cached with provenance instead. M42/M43/M44 had completed in that job and were not repeated.

## 14. Environment

`R_env` **unchanged**: `R_env_PRE_PHASE6.yml` and `R_env_POST_PHASE6.yml` byte-identical (asserted in
the job script). decoupleR 2.12.0 and OmnipathR 3.14.0 were **already installed**, so nothing was
added for Phase 6. Networks cached in `external/networks/` (CollecTRI raw + processed, PROGENy
top-500).

## 15. Deliverables

**13 final figures** in `results/phase6/figures/final/` (PDF + PNG) · **12 final tables** in
`results/phase6/tables/final/` · manifest `results/phase6/phase6_manifest.json` (31 sections) ·
`reports/phase6/{PHASE6_FEASIBILITY.md, ROBUSTNESS_REPORT.md, milestones/M42–M50}` ·
`reports/FIGURE_INDEX.tsv` 706 → 732 rows.

## 16. What Phase 6 answers

1. **Are programs strongly associated with CNA clones?** **No** — median η² 0.059.
2. **Do individual clones occupy multiple programs?** **Yes in MPNST_1** (4 of 7 clones span ≥ 2),
   **no in MPNST_2 and MPNST_4**.
3. **How much within-clone diversity?** ~94% of program variance is within-clone; effective programs
   per clone 1.00–2.83 depending on patient.
4. **Is it robust to clone size, depth and cell cycle?** Yes to clone size (η² identical at 20/50/100
   cells), depth (ρ ≤ 0.23) and clone size (ρ ≤ 0.17); **dispersion is confounded with
   High-confidence fraction at ρ = 0.72** and that is stated.
5. **Clone-constrained, plastic, or mixed?** **Mixed (Model C)**, with clones transcriptionally
   near-interchangeable within a patient.
6. **Which TFs distinguish the programs?** E2F4/E2F1/MYC (Cycling); HIF1A/ATF4/HSF1 (translation);
   HMGA2/MYC/SP1/STAT6 (ECM); ATF3/HSF2 (Schwann-like); GATA3/HOXA9 (neuronal).
7. **Which pathways?** PROGENy Hypoxia for P4, EGFR for P3; Hallmark E2F/MYC targets for P7, EMT for
   P3, MTORC1/glycolysis for P5.
8. **Do broad CNA events show coordinated transcriptomic consequences?** **Yes, in 68% of 263
   evaluable events** — as internal consistency, not validation.
9. **What recurs across patients?** The regulatory associations of Cycling and translation-stress;
   the clone-decoupling of transcriptional phenotype. **Program identity itself does not recur**
   (Phase 5).
10. **What is limited by MPNST_3?** Every clone-based conclusion, which rests on 3 patients.
