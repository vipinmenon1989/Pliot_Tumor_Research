# LochNESS MPNST Report

**Phase 3 · Milestone M24 · MPNST** · *Generated 2026-09-03 · SLURM JobID 19896064*

> **LochNESS is not a ligand–receptor method.** Nothing below is evidence of cell–cell
> communication. LochNESS describes local composition in receiver-state space.

Design written before execution: `reports/phase3/LOCHNESS_MPNST_DESIGN.md`.
Implementation audit: `reports/phase3/LOCHNESS_IMPLEMENTATION_AUDIT.md`.

---

## 1. What was run

Official MMCA formulation, implemented natively in R:
`lochNESS_i = (n_target_neighbours / k) / global_fraction_in_reference − 1`, with
`k = round(0.5·√N)`, L2-normalised `postint_harmony` restricted to each lineage, exact
`FNN::get.knnx`, **disjoint query/reference sets so same-sample neighbours are excluded**, and
the global fraction taken over the reference set.

**Context label** (from M22, computed from **tumour cells only**): the prioritised
tumour-derived ligand was **APP**, and samples were split on mean APP expression in
`MPNST-Tumor` cells:

| Sample | Mean APP in tumour | Detected in | Context |
| --- | ---: | ---: | --- |
| MPNST_1 | 1.081 | 76.6% of 2,428 tumour cells | **high** |
| MPNST_4 | 0.985 | 62.8% of 664 | **high** |
| MPNST_2 | 0.682 | 43.3% of 97 | low |
| MPNST_3 | 0.666 | 51.9% of 231 | low |

Receiver cells never contribute to this label — that separation is what keeps the design
non-circular (design §7).

---

## 2. Results

| Lineage | N | k | mean lochNESS | high-context | low-context | frac > 0 | Null A p | Null B p | Harmony-vs-PCA ρ | LOSO ρ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Macrophage | 3,065 | 28 | −0.163 | −0.195 | −0.104 | 0.391 | **0.667** | 0.730 | 0.419 | 0.497 |
| Fibroblast | 5,064 | 36 | +0.104 | +0.130 | +0.028 | 0.685 | **0.833** | 0.830 | 0.023 | 0.347 |
| Endothelial | 960 | 15 | −0.176 | −0.158 | −0.222 | 0.352 | **0.333** | 0.360 | 0.219 | 0.448 |
| CD8-T | 382 | 10 | −0.311 | −0.574 | +0.020 | 0.440 | **0.333** | 0.420 | −0.265 | 0.310 |

Null A = exact sample-level permutation (6 labelings, **p-floor 0.167**).
Null B = 100 cell-level context shuffles preserving per-sample sizes.
LOSO ρ = mean pairwise Spearman across leave-one-reference-sample-out variants.

---

## 3. Interpretation — a negative result

**No receiver lineage shows neighbourhood structure associated with the tumour-derived APP
signalling context beyond what the null allows.** The best descriptive p-value is 0.333
(Endothelial and CD8-T), twice the 0.167 floor, and the two largest lineages are further from
significance than that.

This is a **legitimate and reportable negative** (Phase 3 §68). Its meaning, stated
precisely: *the CCC analysis nominated APP as the leading tumour-derived signal, but receiver
cells are not organised in transcriptional-neighbourhood space according to how much APP their
tumour of origin expresses.* Candidate reasons, none testable here:

1. Four samples split 2/2 give almost no resolution — the design cannot detect a modest effect.
2. APP expression varies only ~1.6-fold across tumours; the contrast is weak.
3. Harmony was fitted on `sample_id`, so any sample-associated receiver structure was
   explicitly targeted for removal before this analysis (design §4).
4. There may genuinely be no such structure.

## 4. Two additional findings that undercut confidence in the score itself

**Representation sensitivity is poor.** Spearman correlation between the Harmony-space score
and the within-lineage-PCA score is 0.42 (Macrophage), 0.22 (Endothelial), 0.02 (Fibroblast)
and **−0.27 (CD8-T)**. The score depends heavily on the embedding, and for CD8-T the two
representations disagree in direction. Any positive result would have needed to survive this;
none did.

**Leave-one-sample-out stability is moderate at best** (ρ 0.31–0.50). Which tumour serves as
reference materially changes the score — expected with four samples, but it means per-cell
LochNESS values here should not be treated as stable quantities.

## 5. The CD8-T signal, and why it is not claimed

CD8-T shows the largest separation by context (high −0.574 vs low +0.020). It is nonetheless
**not** reported as a finding: Null A p = 0.333, k = 10 (small), n = 382 with only 26 cells in
MPNST_2, representation correlation is negative (−0.265), and LOSO ρ is the lowest of the four
(0.310). Every stability check argues against it.

## 6. Limitations

1. Four samples: the exact null cannot go below p = 0.167. Everything is descriptive.
2. Sample = patient = batch; residual structure is not attributable to biology.
3. The embedding was corrected on the exclusion variable (conservative, but not neutral).
4. The context label is a 2/2 split of a continuous quantity.
5. Partial circularity in *pathway selection* remains — APP was nominated by CCC analyses
   that used receiver data — although the label and the score are independent (design §7).
6. **LochNESS is not evidence of communication**, positive or negative. This negative result
   does **not** refute the APP–CD74 LR finding; it says receiver-state *architecture* does not
   track the context.

## 7. Files

`results/phase3/lochness/lochness_scores.tsv` · `lochness_summary.tsv` ·
`lochness_permutation_summary.tsv` · `lochness_leave_one_sample_out.tsv` ·
`lochness_context_labels.tsv` · figures `results/phase3/figures/M24/M24_01..05`.
