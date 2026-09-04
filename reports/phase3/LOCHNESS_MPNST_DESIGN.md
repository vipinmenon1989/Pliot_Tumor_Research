# LochNESS MPNST Design

**Phase 3 · Milestone M24 · MPNST project**
*Written 2026-09-03, BEFORE any LochNESS run, as required by Phase 3 §31.*

Prerequisite reading: `reports/phase3/LOCHNESS_IMPLEMENTATION_AUDIT.md` (publication vs
official MMCA R vs our perturb-seq Python implementation).

---

## 1. Scientific question

> Within a given receiver lineage, is there local transcriptional-neighbourhood structure
> associated with the tumour-derived signalling context of the sample a cell comes from?

Concretely, for the macrophage compartment: do macrophages from tumours with a **high**
tumour-derived signalling context for a prioritised pathway occupy different neighbourhoods
of macrophage state space than macrophages from **low**-context tumours?

**What this is not.** LochNESS is not a ligand–receptor method and this design does not
treat it as one. It contributes an orthogonal statement about *receiver-state structure*,
downstream of CCC:

```text
CCC (M20-M22)  ->  prioritised tumour-derived ligand / pathway
                        |
                        v
             sample-level signalling CONTEXT label (from TUMOUR cells)
                        |
                        v
       LochNESS in RECEIVER transcriptional space (this design)
                        |
                        v
        does receiver-state structure track the signalling context?
```

---

## 2. Mapping the original design onto MPNST

The original contrast is *mutant embryo vs wildtype embryo*, with neighbourhoods computed
per developmental trajectory and same-embryo neighbours excluded. MPNST has **no genotype
contrast and no biological-condition variable at all** (Phase 2 M10 audited all 93 metadata
columns). The mapping is therefore:

| MMCA (Huang et al.) | MPNST adaptation | Justification |
| --- | --- | --- |
| Embryo (`RT_group`) — the unit of exclusion | **`sample_id`** — the unit of exclusion | Both are the independent biological/technical replicate unit. In MPNST, sample = patient = batch. |
| Genotype (`Mutant` vs `WT`) — the target/reference contrast | **Sample-level tumour-derived signalling context** (`high` vs `low`) | Provides a non-trivial group label that is *not* identical to the exclusion unit, exactly as genotype is not identical to embryo. |
| Developmental trajectory — the stratification | **Receiver lineage** (Macrophage, CD8-T, Fibroblast, Endothelial …) | Both restrict the manifold to a biologically coherent compartment so that "local neighbourhood" is meaningful. Trajectory inference is a Phase 3 prohibition; lineage is the available analogue. |
| Aligned trajectory PCA | **L2-normalised `postint_harmony` embedding**, restricted to the lineage | See §4. |

**The group label is defined at the sample level from TUMOUR cells; the score is computed in
RECEIVER cell space.** That separation is what keeps the design non-circular (§7).

---

## 3. Target and reference populations

For a receiver lineage *L* and a context label `context ∈ {high, low}` assigned per sample:

- **Analysis subset**: all cells of lineage *L*, across all four samples.
- **Query cells**: for each sample *j* in turn, the cells of *L* belonging to sample *j*.
- **Reference cells**: the cells of *L* belonging to **every other sample**.
- **Target group within the reference**: reference cells whose sample is labelled `high`.

For each query cell *i*:

```
lochNESS_i = ( n_high_context_neighbours_i / k ) / global_fraction_high_in_reference − 1
```

- `> 0` → the cell sits in a neighbourhood over-represented for high-context samples
- `≈ 0` → neighbourhood composition matches the reference composition
- `< 0` → under-represented

The global fraction is computed **over the reference set** (the cells that were eligible to
be neighbours), following MMCA exactly.

---

## 4. Neighbourhood representation

**`postint_harmony` (30 dims), L2-normalised with `Seurat::L2Dim`, subset to lineage *L*.**

- L2 normalisation and exact `FNN::get.knnx(..., algorithm="kd_tree")` follow MMCA.
- The Harmony-corrected embedding is used rather than a fresh within-lineage PCA because a
  within-lineage PCA would carry the full uncorrected batch effect, and with sample = batch
  the score would then be a pure batch detector. Using the corrected space means LochNESS
  measures **residual, post-correction** sample-associated structure.

**Trade-off, stated plainly.** Phase 2 M12 measured that Harmony left 78% of nearest
neighbours same-sample overall. Harmony was *fitted on `sample_id`*, which is also the
exclusion unit here. Any positive LochNESS signal is therefore structure that survived an
explicit attempt to remove it — conservative in that direction — but the corrected space is
not sample-neutral, and residual signal cannot be cleanly attributed to biology rather than
incompletely-corrected batch. This limitation is repeated in the results report.

A within-lineage PCA sensitivity analysis is run as a secondary representation so the
dependence on this choice is visible rather than assumed away.

---

## 5. Same-sample exclusion — adopted

**Yes. Neighbours from the query cell's own sample are excluded**, implemented exactly as
MMCA does it: the query and reference matrices are disjoint by construction, so no
same-sample neighbour can be selected.

This is not optional here. `sample_id` is simultaneously the dataset, the patient and the
batch. Without exclusion, each cell's neighbourhood would be dominated by its own sample and
the score would largely re-detect batch structure. Our perturb-seq Python implementation has
**no** such exclusion, which is the single most important reason the R implementation is
authoritative for Phase 3.

---

## 6. Value of *k*

**`k = round(0.5 · sqrt(N_L))`**, where `N_L` is the number of cells in the lineage analysis
subset — the published MMCA rule, adopted unchanged.

Worked values for the CCC-ready lineages:

| Lineage | Cells | k |
| --- | ---: | ---: |
| Fibroblast | 5,064 | 36 |
| MPNST-Tumor | 3,420 | 29 |
| Macrophage | 3,065 | 28 |
| Plasma-cell | 1,709 | 21 |
| Endothelial | 960 | 15 |
| Dendritic | 612 | 12 |
| CD4-T | 606 | 12 |
| Pericyte-VSMC | 435 | 10 |
| CD8-T | 382 | 10 |
| Plasmacytoid-DC | 287 | 8 |
| T-cell-other | 258 | 8 |
| B-cell | 233 | 8 |
| Monocyte | 190 | 7 |
| NK | 106 | 5 |

`k` must also not exceed the size of the smallest reference set. Lineages where
`k > min_j(N_L − N_{L,j})` or where any sample contributes fewer than 10 lineage cells are
declared **not analysable** and reported as such. On these numbers only lineages with
≥ ~200 cells will give stable estimates; **NK (k = 5) and Monocyte (k = 7) are pre-flagged
as unstable** and will be reported with that caveat rather than silently included.

Primary analyses are restricted to the four largest CCC-ready receiver lineages —
**Macrophage, Fibroblast, Endothelial, CD8-T** — chosen for cell count, not for outcome.

---

## 7. Circularity — the risk and how it is handled

Phase 3 §33 forbids defining a cell as "high pathway X" using the very genes later claimed
as independent validation of pathway X.

**Design safeguards:**

1. The context label is computed from **tumour cells only** (mean expression of the
   prioritised tumour-derived *ligand* in `MPNST-Tumor` cells of that sample). The LochNESS
   score is computed in **receiver** cell space. Sender and receiver cell sets are disjoint.
2. The label is assigned at the **sample** level — one value per tumour — never per receiver
   cell, so no receiver cell's own expression enters its own label.
3. The receiver-state characterisation that follows compares **receiver** genes (receptor and
   downstream targets) across LochNESS strata. Those genes are not used to build the label.

**Residual circularity that cannot be removed, and is declared:** the prioritised ligand was
itself nominated by CCC analyses that used receiver receptor expression. So the *choice of
which pathway to test* is not independent of the receiver data, even though the label and the
score are. The correct reading is therefore "receiver-state structure is consistent with the
nominated context", not "LochNESS independently validates the interaction".

---

## 8. Permutation / null strategy

Two nulls, because neither alone is adequate at n = 4 samples.

**Null A — sample-level exact permutation (primary; biologically correct).**
Enumerate **all** assignments of the `high`/`low` context label to the four samples that
preserve the observed group sizes, recompute LochNESS for each, and compare the observed
statistic to that exact distribution. This preserves the sample hierarchy, as Phase 3 §35
requires.

> **Hard limitation, stated up front.** With 4 samples split 2/2 there are only
> `C(4,2) = 6` labelings, and 3 after collapsing complements. **The smallest attainable
> one-sided exact p-value is ≈ 1/6 ≈ 0.17.** No sample-level permutation test can reach
> conventional significance with four samples. Null A is therefore reported as a
> **descriptive rank of the observed statistic within the exact null**, never as a
> significance test.

**Null B — cell-label shuffle preserving per-sample sizes (secondary; finer but weaker).**
Shuffle the sample→context mapping at the cell level while preserving the number of cells per
sample, 100 iterations, seeded per iteration via a derived seed. Neighbourhood geometry and
the same-sample exclusion structure are held fixed, exactly as MMCA's permutation script does.
This has more resolution but does not respect the biological hierarchy, so it can be
anti-conservative.

MMCA's own null shuffles the group label globally across cells, which is closest to Null B;
Null A is an addition made necessary by having four samples rather than dozens of embryos.
**This deviation from MMCA is deliberate and is recorded here rather than applied silently.**

---

## 9. Additional robustness output

Following MMCA's *alternative implementation* (`demo_lochness.R`), a **per-reference-sample
variant** `lochNESS_i` is computed: hold in one reference sample at a time and correlate the
resulting score vectors. High correlation means the score does not hinge on which tumour
served as reference. This doubles as the leave-one-sample-out check M26 requires.

---

## 10. Expected interpretation

| Observation | Reading |
| --- | --- |
| LochNESS ≈ 0 throughout a lineage | Receiver states are not organised by signalling context; the CCC signal, if real, is not accompanied by detectable state restructuring. **A legitimate and reportable negative.** |
| Strong positive/negative regions co-localising with a receiver programme | Receiver-state structure is *consistent with* the nominated context. Suggestive, not confirmatory. |
| Strong structure but scores that flip depending on which sample is the reference | Driven by one patient. Reported as such, not as a general finding. |
| Observed statistic unremarkable within the exact null | No support beyond chance at the resolution four samples allow. |

---

## 11. Limitations

1. **Four samples.** The exact sample-level null cannot yield p < 0.17. Everything is
   descriptive.
2. **Sample = patient = batch.** Residual structure cannot be cleanly separated into biology
   versus incompletely-corrected batch.
3. **The embedding was corrected on the exclusion variable.** Conservative in direction, but
   not sample-neutral (§4).
4. **The context label is coarse** — a 2/2 split of four tumours on a continuous quantity.
5. **Small lineages are unstable.** NK (k = 5) and Monocyte (k = 7) are pre-flagged.
6. **Partial circularity in pathway selection** is unavoidable and is declared (§7).
7. **LochNESS is not evidence of communication.** It describes neighbourhood composition.

---

## 12. Implementation

`scripts/phase3/lochness/lochness_mpnss.R` — native R, following MMCA:
`k = round(0.5·√N)`, `Seurat::L2Dim`, exact `FNN::get.knnx(algorithm="kd_tree")`, disjoint
query/reference sets, global fraction over the reference, both nulls above, and the
per-reference-sample variant.

The perturb-seq Python implementation is run on a deterministic subset as an independent
cross-check; agreement is reported in
`results/phase3/lochness/LOCHNESS_IMPLEMENTATION_COMPARISON.tsv` and
`reports/phase3/LOCHNESS_IMPLEMENTATION_COMPARISON.md`.
