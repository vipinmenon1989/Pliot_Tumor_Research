# MPNST Tumour–Microenvironment Interaction Report

**Phase 3 · Milestone M22 · MPNST** · *Generated 2026-09-03 · SLURM JobID 19896062*

> scRNA-seq ligand–receptor analysis infers communication **potential**. It does not
> establish physical adjacency or direct signalling. All statements below are *predicted*,
> *inferred* or *concordant*, never demonstrated.

---

## 1. Prioritisation rule — printed, not hidden

Lexicographic tiers. **No weighted composite score was invented** (Phase 3 §24, §69):

| Tier | Rule | n |
| --- | --- | ---: |
| **P1_highest** | all testable frameworks agree **and** 4/4 patients **and** expression support | 1,587 |
| **P2_high** | all testable frameworks agree **and** ≥3/4 patients **and** expression support | 1,435 |
| **P3_moderate** | ≥2 frameworks **and** ≥3/4 patients **and** expression support | 163 |
| P4_moderate_low | ≥2 frameworks **and** ≥2/4 patients **and** expression support | 1,513 |
| P5_low | ≥2/4 patients **and** expression support | 8,697 |
| P6_exploratory | retained but flagged | 23,091 |

*Expression support* = ligand **and** receptor each detected in ≥10% of the relevant
population, with complexes scored as the **minimum** across subunits.

Of 7,080 tumour-centric supported interactions, **297 tumour→TME** and **224 TME→tumour**
reach P1–P3.

---

## 2. MPNST-Tumor → immune (Question 1)

The dominant picture is a **myeloid-directed, immunosuppressive, antigen-presenting tumour**.

### 2.1 Three-framework, 4/4-patient tumour→myeloid signals

| Ligand → Receptor | Receivers | Reading |
| --- | --- | --- |
| **APP → CD74** | Macrophage, Dendritic, Monocyte, B-cell, pDC | The single most consistent tumour-derived signal in this dataset — five receivers, three frameworks, all four tumours |
| **CD99 → PILRA** | Monocyte, Dendritic, Macrophage | PILRA is a myeloid **inhibitory** receptor |
| **ANXA1 → FPR1** | Monocyte, Dendritic, Macrophage | Myeloid chemotaxis and polarisation |
| **THBS1 → CD36** | Monocyte | Scavenger-receptor engagement |
| **HLA-F → LILRB1 / LILRB2** | Monocyte | Inhibitory LILR engagement; 4 evidence streams |

### 2.2 Tumour → lymphoid

| Ligand → Receptor | Receiver | Reading |
| --- | --- | --- |
| **HLA-A / HLA-E → CD8A / CD8B** | CD8-T | Class I presentation retained |
| **HLA-E → KLRC1 (NKG2A)** | NK | **Inhibitory checkpoint**; 4 evidence streams, 4/4 patients |
| **BAG6 → NCR3 (NKp30)** | NK | NK-activating-receptor ligand |
| **CD58 → CD2** | NK | Adhesion/co-stimulation |

**The lymphoid picture is mixed and worth stating plainly:** the tumour presents class I and
supplies an NKp30 ligand (potentially activating) while simultaneously engaging NKG2A and
LILRB1/2 (inhibitory). This is a balance, not a one-way suppression, and this dataset cannot
resolve which dominates functionally.

### 2.3 Tumour → stromal / endothelial (Question requiring separate treatment)

| Ligand → Receptor | Receiver | Reading |
| --- | --- | --- |
| **VEGFA → KDR / FLT1 / NRP1** | Endothelial | Canonical angiogenesis; 4 evidence streams, 4/4 patients |
| **JAG1 → NOTCH4** | Endothelial | Notch-mediated vascular signalling |
| **JAG1 → NOTCH3** | Pericyte-VSMC | Canonical mural recruitment |
| **COL1A2 → ITGA1/2/3/9_ITGB1** | Endothelial | Tumour-derived ECM engaging endothelial integrins |
| **FGF2 → FGFR1** | Fibroblast | 4 evidence streams — the only tumour→fibroblast signal with full support |
| **SLIT2 → ROBO1** | Fibroblast | Axon-guidance family repurposed in stroma |

---

## 3. Immune / TME → MPNST-Tumor (Question 2)

Not omitted, and biologically informative:

| Ligand → Receptor | Sender | Reading |
| --- | --- | --- |
| **JAG1 / JAG2 / DLL4 → NOTCH2** | Pericyte-VSMC, Endothelial | A **perivascular Notch niche** feeding tumour NOTCH2 |
| **CRTAM → CADM1** | CD8-T | T-cell adhesion to tumour |
| **TNF → TNFRSF1A** | CD4-T | Direct inflammatory input to tumour |
| **NCAM1 → FGFR1** | NK | Adhesion-molecule-to-RTK crosstalk |
| **SEMA4D → PLXNB2** | NK | Semaphorin signalling into tumour |
| **COL1A1/COL1A2/COL6A2/FN1 → ITGAV_ITGB8** | Fibroblast, Pericyte | Stromal ECM engaging a latent-TGF-β-activating integrin |

**The Notch axis runs in both directions** — tumour JAG1 → pericyte NOTCH3 and pericyte/
endothelial JAG1/JAG2/DLL4 → tumour NOTCH2 — suggesting a reciprocal perivascular circuit
rather than one-way signalling.

---

## 4. Literature interpretation (Question 4 context)

Full table with DOIs: `results/phase3/tables/CCC_LITERATURE_EVIDENCE.tsv`.

### Established in MPNST/NF1, and recovered here

- **NOTCH (JAG1/NOTCH2/NOTCH3).** NOTCH3, NOTCH4, DLL1, DLL3, JAG1 are among genes
  upregulated on PRC2 loss in NF1-deficient Schwann-lineage cells and MPNST; cleaved NOTCH
  with higher HES1 is found in PRC2-deficient MPNST, and NICD transforms primary Schwann
  cells ([Neuro-Oncol Adv 2024, doi:10.1093/noajnl/vdae188](https://academic.oup.com/noa/article/6/1/vdae188/7888919);
  [Oncogene 2004](https://www.nature.com/articles/1207068)). Our data add a **cellular
  source**: the perivascular compartment.
- **CSF1 → CSF1R.** CSF1R inhibition (pexidartinib) is in clinical development for
  neurofibroma. Recovered at P1 with 2 frameworks and 4/4 patients — but with **pericyte and
  fibroblast, not tumour, as the strongest senders**, a specific and testable prediction.
- **SPP1 → CD44.** Published NF1 work reports that SPP1-CD44 is tumour-autocrine in
  plexiform neurofibroma while in MPNST a macrophage subset becomes the dominant SPP1 source
  ([npj Precis Oncol 2025, doi:10.1038/s41698-025-01078-2](https://www.nature.com/articles/s41698-025-01078-2)).
  We recover SPP1-CD44 with macrophage and endothelial senders, consistent with that reported
  shift, though at 2 frameworks / 3-of-4 patients it is weaker here than APP-CD74.

### Established in other cancers, apparently not yet characterised in MPNST

- **APP → CD74.** Well documented as a TAM-immunosuppression axis in gastric cancer,
  testicular tumours and glioblastoma, where APP on tumour cells engages macrophage CD74,
  suppresses phagocytosis and drives M2-like polarisation
  ([J Pathol 2024, doi:10.1002/path.6343](https://pathsocjournals.onlinelibrary.wiley.com/doi/full/10.1002/path.6343);
  [npj Precis Oncol 2025](https://www.nature.com/articles/s41698-025-01268-y)).
  We found no MPNST-specific report. Given that it is the strongest and most reproducible
  tumour-derived signal here, this is the most interesting candidate the analysis produced.
- **ANXA1 → FPR1**, **HLA-E → KLRC1**, **HLA-F → LILRB1/2**, **BAG6 → NCR3**: established
  immunology, not characterised in MPNST.

### Apparently less characterised anywhere in this context

- **CD99 → PILRA.** Supported by all three frameworks in all four tumours; no MPNST-specific
  literature found.

> **We do not claim novelty on the basis of a search that failed to find a paper.** Each row
> in the literature table records *characterisation status*, not a novelty claim.

---

## 5. Limitations

1. Communication **potential** only — no adjacency, no directionality proof.
2. `sample_id` = dataset = patient = batch. No between-tumour claim is supportable.
3. The Phase 2 malignant fraction is itself uncertain (17.35% conservatively, up to ~49% if
   the fibroblast clusters are Mes-NC-like malignant). **Every tumour-centric result inherits
   that uncertainty**, and the tumour↔fibroblast axis is the most exposed.
4. 53.6% of supported interactions rest on one patient (see `ROBUSTNESS_REPORT.md`).
5. Expression support uses a 10% detection floor; complexes scored by their weakest subunit.
6. No condition variable exists, so no condition-stratified interaction analysis was possible.

## 6. Files

`results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv` ·
`results/phase3/tables/{tumor_to_immune,immune_to_tumor,tumor_to_stromal}_interactions.tsv` ·
`CCC_LITERATURE_EVIDENCE.tsv` · figures `results/phase3/figures/M22/`.
