# Phase 4 — Inferred CNA Architecture and Clone Structure

> **SCEVAN infers copy number from gene-expression patterns. This is not DNA sequencing.**
> Inference is reliable for **broad chromosomal, arm-level and large-segment** events only.
> **Single-gene CNV calls are not asserted anywhere in this report** (§51, §65B). Where a
> segment happens to contain a gene of interest, the statement is that the *segment* is altered
> and the locus falls inside it — not that the gene was found deleted.

---

## 1. Broad copy-number events

59 broad clonal segments (≥10 Mb, non-neutral) across the four patients:

| Patient | broad gains | broad losses |
| --- | ---: | ---: |
| MPNST_1 | 15 | 15 |
| MPNST_2 | 9 | 5 |
| MPNST_3 † | 1 | 2 |
| MPNST_4 | 5 | 7 |

† MPNST_3's SCEVAN partition failed the immune sanity check, so its "tumour" set is largely
immune cells and its CNA profile carries no weight. Recurrence is therefore reported **twice**.

## 2. Recurrent broad events

| Chromosome | direction | patients (reliable 3) | patients (all 4) |
| --- | --- | :--: | :--: |
| **18** | **loss** | **3 / 3** | 4 / 4 |
| **2** | **gain** | **3 / 3** | 3 / 4 |
| **7** | **gain** | **3 / 3** | 3 / 4 |
| 3 | gain | 2 / 3 | 2 / 4 |
| 3 | loss | 2 / 3 | 2 / 4 |
| 5 | gain | 2 / 3 | 2 / 4 |
| 6 | loss | 2 / 3 | 2 / 4 |
| 11 | loss | 2 / 3 | 2 / 4 |
| 15 | gain | 2 / 3 | 2 / 4 |
| 19 | loss | 2 / 3 | 2 / 4 |
| 22 | loss | 2 / 3 | 2 / 4 |

**Three events are shared by every reliable patient**: whole-chromosome-scale **loss on
chromosome 18** and **gains on chromosomes 2 and 7**. With n = 3 contributing patients this is
a *consistent pattern*, not an established recurrence rate, and no frequency estimate for MPNST
as a disease follows from it.

Note the wording rule (§50): these are **shared or recurrent CNA patterns across patients**.
They are **not** evidence of a shared clone — clones are patient-scoped by construction and
similar architecture between two different tumours is not the same clone.

## 3. Segments containing NF1 and NF2

Reported with deliberate caution.

| Patient | Chr | segment | CN | segm.mean | width | contains |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| MPNST_4 | 17 | 0.44–38.51 Mb | **1** | −0.078 | **38.1 Mb** | **NF1 (17q11.2, ≈29.4 Mb)** |
| MPNST_4 | 17 | 38.67–50.87 Mb | 4 | +0.163 | 12.2 Mb | — |
| MPNST_4 | 17 | 50.96–83.10 Mb | 3 | +0.065 | 32.1 Mb | — |
| MPNST_1 | 22 | 16.60–50.80 Mb | **1** | −0.130 | **34.2 Mb** | **NF2 (22q12.2, ≈29.6 Mb)** |
| MPNST_4 | 22 | 16.60–50.78 Mb | **1** | −0.056 | 34.2 Mb | **NF2 (22q12.2)** |

MPNST_4 shows a **38.1 Mb single-copy segment spanning 17p and proximal 17q**, whose interval
contains the NF1 locus, alongside gain across the remainder of chromosome 17. Two patients show
a **34.2 Mb single-copy segment on 22q** containing the NF2 locus.

These are consistent with what is expected in a nerve-sheath malignancy, and that consistency
is worth stating — but the analysis **cannot** distinguish a focal NF1 deletion from
large-segment loss that happens to include NF1, and it does not attempt to. Confirming NF1 or
NF2 status requires DNA sequencing or a targeted assay, neither of which is in scope.

The recurrent **chromosome 7 gain** likewise contains EGFR (7p11.2), and MPNST is reported to
overexpress EGFR — but the same caution applies, and Phase 3 did not find EGFR among its
supported interactions.

## 4. Clone and subclone structure

22 subclones. Labels are **patient-scoped** (`<PATIENT>_cloneN`).

| Patient | subclones | largest clone (share of that patient's malignant cells) |
| --- | ---: | --- |
| MPNST_1 | 7 | clone4, 516 cells (26.4%) |
| MPNST_2 | 4 | clone3, 234 cells (33.7%) |
| MPNST_3 † | 3 | clone3, 674 cells (47.7%) — **immune, not tumour** |
| MPNST_4 | 8 | clone4, 1,158 cells (24.5%) |

No patient is dominated by a single subclone: the largest clone holds 24–34% of the malignant
compartment in the three reliable patients. Subclonal structure is genuinely polyclonal in each.

## 5. What the clones are made of — the most direct evidence in Phase 4

Clone composition in Phase 2 terms is the least interpretive way to answer the fibroblast
question, because SCEVAN's clones are constructed from CNA profiles without reference to any
Phase 2 label.

**MPNST_2 — all four clones are fibroblast-dominated:**

```text
clone3  234 cells   Fibroblast=229  Plasma-cell=3  MPNST-Tumor=2
clone2  186 cells   Fibroblast=160  Plasma-cell=21 Uncertain=3  ...
clone1  158 cells   Fibroblast=138  Plasma-cell=13 MPNST-Tumor=2 ...
clone4  116 cells   Fibroblast=110  Plasma-cell=3  MPNST-Tumor=2 ...
```

**MPNST_4 — seven of eight clones are fibroblast-dominated:**

```text
clone4 1158  Fibroblast=873  Plasma-cell=141 MPNST-Tumor=116 ...
clone2  907  Fibroblast=585  Plasma-cell=166 MPNST-Tumor=63  ...
clone6  772  Fibroblast=534  Plasma-cell=105 MPNST-Tumor=67  ...
clone3  582  Fibroblast=304  MPNST-Tumor=155 Low-quality=90  ...
clone7  549  Fibroblast=205  Low-quality=120 MPNST-Tumor=98  ...
clone5  326  Fibroblast=200  Plasma-cell=56  MPNST-Tumor=37  ...
clone1  288  Fibroblast=184  MPNST-Tumor=51  Low-quality=25  ...
clone8  152  Endothelial=146 MPNST-Tumor=3   Fibroblast=2       <- flagged, see §7
```

**MPNST_1 — clones span all three disputed labels together:**

```text
clone4  516  Candidate-Malignant=235  Fibroblast=219  MPNST-Tumor=62
clone1  442  Fibroblast=235           Candidate-Malignant=170 MPNST-Tumor=36 ...
clone5  370  Candidate-Malignant=137  MPNST-Tumor=134 Plasma-cell=61 Fibroblast=33 ...
clone2  340  MPNST-Tumor=209          Candidate-Malignant=123 Fibroblast=7 ...
clone6  176  Candidate-Malignant=108  MPNST-Tumor=68
clone3   54  Candidate-Malignant=26   Fibroblast=17   MPNST-Tumor=10 ...
clone7   54  Candidate-Malignant=33   MPNST-Tumor=20  Fibroblast=1
```

Every MPNST_1 clone mixes `Candidate-Malignant-Unresolved`, `Fibroblast` and `MPNST-Tumor`.
Marker-based annotation split into three identities what copy-number evidence sees as **one
malignant compartment with seven subclones**. That is the clearest statement of what Phase 4
found.

## 6. MPNST_3's clones are an artefact, and are reported as one

```text
clone3  674  CD4-T=327 CD8-T=135 T-cell-other=102 NK=29 ... (14 annotations)
clone2  402  Plasma-cell=170 Plasmacytoid-DC=153 B-cell=76 ...
clone1  338  Plasma-cell=230 MPNST-Tumor=72 B-cell=13 ...
```

These are **immune lineages**, not tumour subclones: a T/NK cluster, a plasma/pDC/B cluster,
and a plasma-dominated cluster. This independently corroborates the immune sanity failure —
with only 25 confident normal cells, SCEVAN's top-level CNA clustering partitioned MPNST_3
along immune lineage rather than malignant/normal. MPNST_3 contributes **no clone structure**
to Phase 4's conclusions.

## 7. MPNST_4 clone8 — endothelium, flagged not promoted

`MPNST_4_clone8` is 146/152 `Endothelial`. MPNST_4's endothelium is 52.7% SCEVAN-malignant
while endothelium in the other three patients is 0.0–4.6%. Rule R11 sent those 168 cells to
`Ambiguous / Low` rather than promoting them: a majority-malignant canonical endothelial
population in **one** patient is more likely a CNV-inference artefact than a discovery. It is
recorded here rather than dropped, and endothelium was deliberately kept out of the normal
reference set so that this call would be independent of the Phase 3 findings that concern it.

## 8. Cytoband-level events

`SCEVAN_ARM_LEVEL_EVENTS.tsv` holds SCEVAN's own cytoband × subclone onco-heatmap codes:
MPNST_1 103 gains / 86 losses over 530 neutral calls, MPNST_2 39 / 6 over 105, MPNST_4 49 / 48
over 453, MPNST_3 3 / 6 over 38. Used for figure `04_scevan_cna_heatmap` and for describing
which subclones differ from one another; **not** used to make single-locus claims.

## 9. Native SCEVAN figures

SCEVAN's own outputs are indexed rather than redrawn, in
`results/phase4/tables/SCEVAN_NATIVE_FIGURE_INDEX.tsv`: CNA heatmaps with a malignant/normal
annotation track, tumour-cell-only heatmaps, subclone-annotated heatmaps, clonal CN profiles,
consensus profiles, cytoband onco-heatmaps and SCEVAN's own CNA-space and expression-space
UMAPs, for all eight runs.

The **clone phylogeny plot is absent** for every run: SCEVAN's `plotCloneTree` calls ggtree,
which calls `ggplot2:::is.waive()`, removed in ggplot2 4.x. SCEVAN catches the error internally
so no data was lost. ggplot2 was not downgraded — the clone/state relationship is presented as
an explicit cross-tabulation instead (limitation K).
