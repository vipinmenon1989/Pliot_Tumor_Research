# CCC Method Concordance Report

**Phase 3 · Milestone M21 · MPNST** · *Generated 2026-09-03 · SLURM JobID 19896062*

> scRNA-seq ligand–receptor analysis infers communication **potential** from expression and
> associated transcriptional programmes. It does not establish physical cell adjacency or
> direct signalling.

---

## 1. What was combined, and how

Three independent LR frameworks were run **per sample** and pooled:

| Framework | Version | Resource | Support rule |
| --- | --- | --- | --- |
| LIANA | 0.1.14 | Consensus (OmnipathR 3.14.0); methods natmi, connectome, logfc, sca, cellphonedb | `aggregate_rank ≤ 0.05` in ≥1 sample |
| CellChat | 2.2.0.9001 | CellChatDB.human (3,233 interactions, 290 pathways) | present in `subsetCommunication()`, which returns only significant links |
| CellPhoneDB | 5.0.1 | CellPhoneDB v5.0.0, 1,000 permutations | permutation p < 0.05 **and** mean expression > 0 in ≥1 sample |

**Raw scores were never averaged.** A CellChat probability, a CellPhoneDB p-value and a LIANA
rank are different quantities on different scales; each contributes only a boolean support
flag. Canonical key: `sender | receiver | ligand | receptor`, with gene symbols upper-cased
and multi-subunit complexes normalised to a sorted, underscore-joined set.

**NicheNet is not counted here.** It supports a ligand→receiver-*programme* claim, not the
identical receptor pair, so it enters as orthogonal receiver-response evidence in M23/M25.

### Declared partial non-independence

LIANA's method set **includes a CellPhoneDB-style score**. LIANA and the standalone
CellPhoneDB run are therefore not fully independent, and a LIANA+CellPhoneDB agreement is
weaker evidence than a LIANA+CellChat agreement. This is stated rather than glossed over.

---

## 2. The resource-overlap problem, and why every row carries `testable_*`

The three frameworks use different LR resources, so "supported by 1 of 3" is meaningless
unless one knows whether the other two ever *contained* the interaction. Every row therefore
carries `testable_LIANA`, `testable_CellChat`, `testable_CellPhoneDB` alongside the support
flags.

| Testable by | Interaction keys |
| --- | ---: |
| 3 frameworks | **1,347** |
| 2 frameworks | 13,064 |
| 1 framework | **488,435** |
| **Total distinct keys** | **502,846** |

**97% of keys exist in only one framework's resource.** This single fact governs the
interpretation of everything below: most "single-method" results reflect resource
non-overlap, not scientific disagreement.

---

## 3. Concordance classes

| Class | Definition | Interactions |
| --- | --- | ---: |
| **High concordance** | supported by all frameworks that could test it (≥2) | **5,848** |
| **Moderate concordance** | supported by 2 of 3 testable frameworks | 499 |
| **Single-method** | supported by the **only** framework whose resource contained it | 24,654 |
| **Discordant/ambiguous** | testable by ≥2 frameworks but supported by only 1 — **genuine disagreement** | **5,485** |
| **Total supported** | | **36,486** |

## 4. Method overlap

| Support combination | Interactions |
| --- | ---: |
| LIANA only | 15,407 |
| CellPhoneDB only | 10,899 |
| CellChat only | 3,833 |
| LIANA + CellChat | 3,450 |
| CellChat + CellPhoneDB | 1,171 |
| LIANA + CellPhoneDB | 963 |
| **All three** | **763** |

The 763 three-framework agreements are the highest-confidence core of this analysis. Note
they are drawn from a pool of only 1,347 interactions that all three could test — a **57%
agreement rate among jointly testable interactions**, which is high for CCC inference.

## 5. Directional breakdown

| Direction | Supported interactions |
| --- | ---: |
| MPNST-Tumor → TME | 3,673 |
| TME → MPNST-Tumor | 3,407 |
| Tumour ↔ stromal/endothelial | 3,899 |
| TME ↔ TME | remainder |

## 6. Sample recurrence

| Supported in | Interactions |
| --- | ---: |
| 4/4 patients | 4,083 |
| ≥3/4 | 8,966 |
| ≥2/4 | 16,928 |
| **1 patient only** | **19,558 (53.6%)** |

More than half of all supported interactions rest on a single patient. Those are retained
and flagged, never silently promoted (see `ROBUSTNESS_REPORT.md`).

## 7. Negative and disagreement results, kept

- **5,485 discordant interactions** where frameworks that could all test the pair disagreed.
  Full list: `results/phase3/tables/CCC_METHOD_DISAGREEMENT.tsv`.
- **TGFB1 → TGFBR1/TGFBR2** — a pathway one would expect *a priori* in a mesenchymal
  tumour — reached only **1 framework and P5_low**. It is not well supported here. The
  better-supported route to the same biology in this dataset is stromal collagen/FN1 →
  tumour `ITGAV_ITGB8`, an integrin known to activate latent TGF-β.
- **SPP1 → CD44**, prominent in the published NF1/MPNST literature, reached only 2
  frameworks and 3/4 samples here — present but weaker than APP → CD74.

## 8. Limitations

1. Resource non-overlap dominates (§2); the concordance denominator is small.
2. LIANA and CellPhoneDB are partially non-independent (§1).
3. CellChat returns only significant links, so "testable but not supported" cannot be
   distinguished from "not in the resource" for CellChat as cleanly as for the other two.
4. All support flags are per-sample maxima; an interaction supported in one patient counts
   as supported for that framework, with recurrence tracked separately.

## 9. Files

`results/phase3/ccc/concordance/CCC_CONCORDANCE.tsv` ·
`method_overlap_counts.tsv` · `m21_concordance_record.json` ·
figures `results/phase3/figures/M20/M20_01..05`, `results/phase3/figures/M21/M21_06..10`.
