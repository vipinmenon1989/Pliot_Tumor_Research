# Phase 4 — Ambiguous Population Malignancy Report

**The dedicated analysis required by §22, answering the question Phase 4 exists to settle.**

> Inferred copy number, not DNA sequencing. Reliable for broad chromosomal, arm-level and
> large-segment events only. A SCEVAN normal call is not proof of non-malignancy.

---

## 1. The five populations in question

`Fibroblast`, `Candidate-Malignant-Unresolved`, `Pericyte-VSMC`, `MPNST-Tumor` and `Uncertain`.
None of them was ever used as a fixed normal reference — that is the whole point (§16).

| Population | n | SCEVAN malignant | SCEVAN normal | not assessed | **malignant fraction of assessed** | patients with majority malignant | median CNA burden | fraction CNV-elevated |
| --- | ---: | ---: | ---: | ---: | ---: | :--: | ---: | ---: |
| **Fibroblast** | 5,064 | **4,045** | 912 | 107 | **0.816** | **3 / 4** | 0.269 | 0.278 |
| **Candidate-Malignant-Unresolved** | 1,231 | 851 | 316 | 64 | **0.729** | 2 / 4 | **0.328** | **0.658** |
| **MPNST-Tumor** | 3,420 | 1,227 | 2,111 | 82 | 0.367 | 1 / 4 | 0.200 | 0.299 |
| `Pericyte-VSMC` | 435 | 77 | 348 | 10 | 0.181 | 1 / 4 | 0.168 | 0.103 |
| `Uncertain` | 720 | 99 | 487 | 134 | 0.169 | 1 / 4 | 0.125 | 0.088 |

Per-patient malignant fraction (SCEVAN primary):

| Population | MPNST_1 | MPNST_2 | MPNST_3 † | MPNST_4 |
| --- | ---: | ---: | ---: | ---: |
| **Fibroblast** | **0.621** | **0.547** | 0.078 | **0.976** |
| Candidate-Malignant-Unresolved | **0.696** | 0.167 | 0.714 | 0.375 |
| MPNST-Tumor | 0.222 | 0.072 | 0.394 | **0.889** |
| Pericyte-VSMC | 0.016 | 0.069 | 0.077 | **0.890** |
| Uncertain | 0.006 | 0.056 | 0.333 | **0.933** |

† MPNST_3's SCEVAN partition failed the immune sanity gate (Amendments A1/A2). Its column is
shown for completeness but carries no weight, and malignant promotions from it were disabled.

---

## 2. Fibroblast — the central question (§22, §43)

### Were the Phase 2 fibroblast-like populations true fibroblasts, malignant Mes-NC-like cells, or mixed?

**Mixed, and predominantly malignant.** Of 5,064 `Fibroblast` cells:

```text
4,036  refined Malignant       (79.7%)   -> 1,217 High confidence, 2,819 Moderate
  908  refined Non-malignant   (17.9%)   -> 101 High confidence (rule R12), 807 Moderate
  120  refined Ambiguous        (2.4%)
```

The evidence is not a single method's say-so:

1. **Replication across patients.** Majority-malignant in **three of four** patients —
   MPNST_1 0.621, MPNST_2 0.547, MPNST_4 0.976. The one exception is the sample whose SCEVAN
   run failed independently. **It is not a single-patient finding.**
2. **Threshold independence.** 4,036 malignant fibroblasts at *every* `pop_frac_low` from 0.15
   to 0.40. The conclusion does not depend on a cut-off choice.
3. **Reference-strategy independence.** Primary↔sensitivity agreement is 0.995, 0.966 and 0.997
   in the three contributing patients. It does not depend on how the CNV baseline was set.
4. **Clone composition, which is the most direct evidence.** SCEVAN's subclones are *built
   from* these cells. **All four MPNST_2 clones** are fibroblast-dominated (229/234, 160/186,
   138/158, 110/116), as are **seven of eight MPNST_4 clones** (873/1158, 585/907, 534/772,
   304/582, 205/549, 200/326, 184/288). A CNV-defined subclone composed of cells that
   marker-based annotation read as fibroblasts is precisely the signature of a malignant
   Mes-NC-like population.
5. **CNA burden.** Median 0.269 with 27.8% of cells above their own sample's confident-normal
   95th percentile.

**What survives as genuinely non-malignant fibroblast**: 908 cells, concentrated in MPNST_1
(303) and MPNST_2 (503), of which 101 reach High confidence via R12 (predominantly-normal
stratum plus concordant lineage evidence). A real fibroblast compartment exists; it is simply
much smaller than Phase 2's 5,064.

**The honest caveat.** MPNST_4 contributes 2,887 of the 4,036 promoted cells, and MPNST_4 is
the sample rescued by Amendment A2. If that amendment is rejected and the stricter flat gate
applied, the promoted count falls to roughly 1,149 from MPNST_1 and MPNST_2 — **still a
majority-malignant fibroblast compartment in two independent patients, but a much smaller
absolute number.** The direction of the conclusion is robust to the amendment; its magnitude
is not.

---

## 3. Candidate-Malignant-Unresolved — resolved

Phase 2 created this label deliberately provisional: "the only evidence is neural-programme
expression in clusters that are 93–99% one patient", and it was **excluded from Phase 3**
entirely. Phase 4 was written to settle it.

```text
836  refined Malignant      (67.9%)  -> 825 High confidence (R17), 11 Moderate
395  refined Ambiguous      (32.1%)
  0  refined Non-malignant
```

**Not one cell in this population is called non-malignant.** It carries the **highest CNA
burden of any population** (median 0.328) and the **highest fraction of CNV-elevated cells
(65.8%)** — higher than `MPNST-Tumor` itself. Majority-malignant in MPNST_1 (0.696) and
MPNST_3 (0.714).

Phase 2's caution was appropriate on the evidence it had, and Phase 4 vindicates the decision
to keep these cells visible and promotable rather than either absorbing or discarding them.
1,196 of the 1,231 cells are MPNST_1, so this remains **substantially a one-patient
population** — the malignancy call is resolved, its generality is not.

---

## 4. MPNST-Tumor — the uncomfortable result

SCEVAN corroborates only **36.7%** of the cells Phase 2 called malignant. 1,405 become refined
Malignant; **2,015 become Ambiguous** through rule R5.

This is the one place where Phase 4 *reduces* confidence in a Phase 2 conclusion, and it is
reported rather than buried. Three readings are available and this dataset cannot separate
them:

1. **Copy-number-quiet malignant cells.** §65C anticipates exactly this, which is why a SCEVAN
   normal call is never treated as proof of non-malignancy and why these cells became
   `Ambiguous` rather than `Non-malignant`.
2. **SCEVAN under-calling where the malignant compartment is small.** MPNST_2 has 97
   Phase 2 `MPNST-Tumor` cells and the lowest corroboration (0.072).
3. **Genuine Phase 2 over-assignment** in some clusters — the M16 review had already downgraded
   the "Schwann-lineage tumour-like" label to Low confidence.

Supporting the first reading: the population's **marker evidence is `malignant-consistent`** —
its Schwann/neural-crest panel score is above the object-wide median while its non-malignant
panel score is not. Lineage evidence and copy-number evidence genuinely disagree here, and the
`Ambiguous` call records that rather than picking a winner.

MPNST_4 is the exception at 0.889 corroboration, and it is also where clone membership shows
`MPNST-Tumor` cells distributed across seven clones alongside the fibroblast-derived ones — a
single coherent malignant compartment.

---

## 5. Pericyte-VSMC — largely genuine

`348 / 435` refined Non-malignant (344 at High confidence via R12); 65 Malignant; 22 Ambiguous.
Malignant fraction is 0.016–0.077 in three patients and 0.890 in MPNST_4 alone — **one-patient,
so not a recurrent finding**. Phase 2's decision to keep pericytes out of the malignant
compartment was right, and the Phase 3 perivascular Notch circuit rests on a compartment that
copy-number evidence supports as non-malignant in MPNST_1, MPNST_2 and MPNST_3.

## 6. Uncertain — mostly non-malignant

`476 / 720` refined Non-malignant (R19, all Moderate — never High, because Phase 2 could not
establish these cells' identity); 92 Malignant; 152 Ambiguous. Malignant fraction 0.006–0.056
in two patients, 0.933 in MPNST_4. Phase 2's two `Uncertain` clusters (C12
hypoxic-versus-perineurial, C23 low-complexity) are largely **not** malignant, so excluding
them from Phase 3 cost little.

## 7. Endothelial — flagged, not promoted

Not one of the five audited populations, but it must be recorded: MPNST_4's endothelium is
52.7% SCEVAN-malignant and forms its **own subclone** (`MPNST_4_clone8`, 146 of 152 cells
endothelial). Rule R11 sent those 168 cells to `Ambiguous / Low` rather than promoting them,
because a majority-malignant canonical endothelial population in one patient is more likely a
CNV-inference artefact than a real finding — but it is **flagged, not silently overridden**.
Endothelium was deliberately excluded from the normal-reference set precisely so this call
would be independent (§16). In MPNST_1, MPNST_2 and MPNST_3 endothelium is 0.000–0.046
malignant.

---

## 8. Summary table

| Population | Phase 2 status entering Phase 4 | Phase 4 verdict | Confidence |
| --- | --- | --- | --- |
| `Fibroblast` | non-malignant | **predominantly malignant Mes-NC-like (4,036 / 5,064)**, with a real 908-cell non-malignant remainder | High — replicated in 3 patients, threshold- and reference-independent |
| `Candidate-Malignant-Unresolved` | provisional, excluded from Phase 3 | **malignant (836 / 1,231, none non-malignant)**, highest CNA burden of any population | High for the call; limited by being 97% one patient |
| `MPNST-Tumor` | malignant on marker evidence | **only 1,405 / 3,420 corroborated; 2,015 Ambiguous** | Low — the most fragile number in Phase 4, threshold-dependent |
| `Pericyte-VSMC` | non-malignant | **confirmed non-malignant (348 / 435)** | High in 3 patients |
| `Uncertain` | unresolved, excluded from Phase 3 | **mostly non-malignant (476 / 720)** | Moderate |
| `Endothelial` | non-malignant | non-malignant in 3 patients; **flagged Ambiguous in MPNST_4** | Moderate |
