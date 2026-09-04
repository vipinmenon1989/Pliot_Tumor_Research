# Cell-Type Annotation Report

**Phase 2 · Milestone M15 · MPNST single-cell integration**
*Generated: 2026-09-02 · SLURM JobID 19886713*

Machine-readable evidence: **`results/phase2/annotation/ANNOTATION_EVIDENCE.tsv`**
Auditable annotation map: **`config/phase2/annotation_map_M15.tsv`**

---

## 1. Method

Annotation is applied **per cluster** to the M13 primary solution (Harmony resolution 1.0,
26 clusters, 19,716 cells) and is hierarchical:

| Level | Content |
| --- | --- |
| **Level 1** | Broad compartment (Malignant/tumour, T/NK, Myeloid, B/Plasma, Endothelial, Fibroblast/Stromal, Other, Uncertain) |
| **Level 2** | Canonical cell type, only where supported |
| **Level 3** | State/subtype, only where marker evidence is strong |
| **Confidence** | High / Moderate / Low / Uncertain |

Every label, its positive markers, its negative markers, its conflicting evidence, its
confidence and its citation live in one reviewable TSV (`config/phase2/annotation_map_M15.tsv`).
The annotation script applies that map; it does not infer labels. This keeps the
scientific judgement visible and auditable rather than buried in code.

**No CNV inference was performed** — it is outside the authorised Phase 2 scope. Every
malignant call is therefore expression- and provenance-based, and is labelled
conservatively.

---

## 2. Result

### Level 1 — broad compartments

| Compartment | Cells | % |
| --- | ---: | ---: |
| Fibroblast/Stromal | 5,499 | 27.9% |
| **Malignant / tumour** | **4,651** | **23.6%** |
| Myeloid | 4,154 | 21.1% |
| B/Plasma | 1,942 | 9.8% |
| T/NK | 1,352 | 6.9% |
| Endothelial | 960 | 4.9% |
| Uncertain | 720 | 3.7% |
| Other (technical) | 438 | 2.2% |
| **Total** | **19,716** | 100% |

### Level 2 — cell types

| Cell type | Cells | Clusters |
| --- | ---: | --- |
| Fibroblast | 5,064 | C0, C3, C5, C25 |
| Macrophage | 3,065 | C1, C6 |
| Plasma cell | 1,709 | C2 |
| MPNST-like malignant (SCP-like) | 1,680 | C8, C9 |
| T cell | 1,352 | C4 |
| Candidate malignant | 1,231 | C13, C17, C18 |
| Schwann-lineage tumour-like | 1,011 | C7 |
| Endothelial cell | 960 | C11, C19 |
| Uncertain | 720 | C12, C23 |
| Dendritic cell (cDC2) | 612 | C10 |
| MPNST-like malignant (NC-like) | 440 | C14 |
| Low-quality / mitochondrial-high | 438 | C15 |
| Mural (pericyte/VSMC) | 435 | C16 |
| Cycling tumour-like | 289 | C20 |
| Plasmacytoid dendritic cell | 287 | C21 |
| B cell / plasmablast | 233 | C22 |
| Monocyte | 190 | C24 |

Confidence: **High 10,115 (51.3%)**, Moderate 7,650 (38.8%), Low 1,951 (9.9%).

---

## 3. Malignant / MPNST annotation — the conservative call

Phase 2 forbids the inference "non-immune ⇒ tumour". The malignant assignment here
integrates five independent lines of evidence, and the resulting labels are graded by how
much of that evidence is present.

### 3.1 Evidence used

1. **Schwann/neural-crest lineage markers** in the M14 marker tables.
2. **MPNST-specific literature.** Spatial transcriptomics of peripheral nerve sheath
   tumours resolves malignant cells into **SCP-like** (high **L1CAM**, S100),
   **NC-like** (high **APOD**, S100, **p75/NGFR**) and **Mes-NC-like** (vimentin, MYH9)
   states, and identifies an **MPNST-G1** subgroup with **SHH** pathway activation and
   NC-like character (worse prognosis) versus an **MPNST-G2** subgroup with WNT/β-catenin
   and SCP-like character ([Neuro-Oncology, doi:10.1093/neuonc/noaf016](https://academic.oup.com/neuro-oncology/advance-article-abstract/doi/10.1093/neuonc/noaf016/7976860)).
3. **Schwann-cell dedifferentiation.** MPNST down-regulates SOX10, CNP, PMP22 and NGFR
   relative to neurofibroma while up-regulating neural-crest stem programmes, so the
   *absence* of classical Schwann markers does not exclude a nerve-sheath lineage
   ([review, PMC11763529](https://pmc.ncbi.nlm.nih.gov/articles/PMC11763529/);
   [Science Advances, doi:10.1126/sciadv.abo5442](https://www.science.org/doi/10.1126/sciadv.abo5442)).
4. **Absence of convincing immune, endothelial or mural identity** (PTPRC−, PECAM1−,
   COL1A1−, LYZ−, CD3D−).
5. **Sample provenance** — used as supporting context only, never on its own.

### 3.2 Graded malignant labels

| Label | Cells | Clusters | Evidence |
| --- | ---: | --- | --- |
| **MPNST-like malignant (SCP-like)** | 1,680 | C8, C9 | **C9 retains MPZ**, the canonical Schwann myelin protein, with GFRA3, CRYAB, NOV(CCN3), RSPO3 and the neural-crest stem marker ABCB5 — the strongest Schwann-lineage evidence in the dataset (canonical-programme fraction 0.986). **C8 expresses L1CAM**, the reported SCP-like marker, with TGFA, PTHLH, MAL and OLFM1. |
| **MPNST-like malignant (NC-like)** | 440 | C14 | **SHH** expression with the Schwann-associated GAL3ST1 (sulfatide synthesis) and KLK6 (myelin-associated), matching the reported SHH-activated NC-like MPNST-G1 subgroup. |
| **Schwann-lineage tumour-like** | 1,011 | C7 | A coherent neural adhesion / axon-guidance programme (NLGN1, TENM2, NDST4, EPHB1, COL25A1, ROBO2, ZNF423) with no immune, stromal or endothelial identity — but no classical Schwann marker. Deliberately weaker wording. |
| **Candidate malignant** | 1,231 | C13, C17, C18 | Neural programmes (NRXN1, NCAM2, KIRREL3, SORCS3; KCNJ6, GPC5, UNC5D, CHRM3) in clusters that are 93–99% one patient, with no immune/stromal identity. **Confidence Low** — patient-private expression alone is weak evidence for malignancy. |
| **Cycling tumour-like** | 289 | C20 | A dominant G2/M programme (TOP2A, UBE2C, NEK2, KIF20A, DLGAP5, HJURP, ASPM) that masks lineage identity. Confidence Moderate; cycling cells of any lineage can collapse into a shared proliferation cluster. |

**Nothing is asserted as definitively malignant.** The strongest available statement is
"MPNST-like malignant", and three clusters carry only "candidate malignant" at Low
confidence.

### 3.3 The largest open question — fibroblast versus MES-like malignant

The four fibroblast clusters (C0, C3, C5, C25; 5,064 cells, 25.7%) carry canonical
fibroblast markers (SFRP2, CXCL14, C7, FBLN1, DCN-family collagens, PI16, COMP, CILP).
But MPNST also contains **Mes-NC-like** malignant cells that express mesenchymal and ECM
programmes and would look fibroblast-like by marker expression alone.

- **For fibroblast:** C3 carries the PI16+ epineurial fibroblast signature, reported as
  distinct from GLUT1+ perineurial cells ([PNAS, doi:10.1073/pnas.1913444117](https://www.pnas.org/doi/10.1073/pnas.1913444117)).
  C25 expresses **CDKN2A**, which is characteristically *deleted* in MPNST — evidence
  against a malignant identity for that cluster.
- **Against:** C5 is 78% MPNST_4 and C0 is 56% MPNST_4; single-patient concentration is
  compatible with patient-private malignant biology.

**Resolving this requires CNV inference, which is explicitly outside Phase 2 scope.**
The consequence for the malignant fraction is stated plainly:

> **The malignant fraction is 23.6% under the conservative annotation. If some or all of
> the fibroblast-programme clusters are in fact Mes-NC-like malignant cells, it could be as
> high as ~49%.** This is the single largest quantitative uncertainty in the Phase 2
> annotation and the highest-priority target for Phase 3 CNV analysis.

For context only, and explicitly not used to make any call: the legacy `orig.anno` labels
inherited from the original analysis assigned ~49% of cells to malignant categories, which
is consistent with the upper bound above.

---

## 4. Immune, stromal and endothelial annotation

| Compartment | Clusters | Decisive markers | Confidence |
| --- | --- | --- | --- |
| **Macrophage, FOLR2+ resident-like** | C1 | RNASE1, CD163, FOLR2, F13A1, C1QA, SIGLEC1, STAB1 | High |
| **Macrophage, inflammatory** | C6 | IL1A, IL1B, CCL3, CCL4, CXCL2, CXCL3, BCL2A1 | High |
| **Monocyte, classical** | C24 | S100A8, S100A9, S100A12, FCN1, CD300E, AQP9 | High |
| **cDC2** | C10 | CD1C, CD1E, FCER1A, CLEC10A, FLT3 | High |
| **Plasmacytoid dendritic cell** | C21 | LILRA4, CLEC4C(BDCA2), SPIB, GZMB | High |
| **T cell, cytotoxic/effector-memory** | C4 | CD3D/E/G, CD2, LCK, GZMA, GZMK, CCL5, KLRB1 | High |
| **Plasma cell** | C2 | JCHAIN(IGJ), MZB1, XBP1, FKBP11, SSR4, SEC11C, ERN1 | High |
| **B cell / plasmablast** | C22 | CD79A, POU2AF1, MZB1, FCRL5 | Moderate |
| **Endothelial, tip/angiogenic** | C11 | ESM1, ANGPT2, APLN, DLL4, FLT1 | High |
| **Endothelial, venous/lymphatic** | C19 | ACKR1(DARC), SELP, CCL14, MMRN1, SLCO2A1 | Moderate |
| **Mural (pericyte/VSMC)** | C16 | MYH11, ACTG2, MYOCD, PLN, RGS5, NOTCH3, HIGD1B | High |

CD163+ tumour-associated macrophages are the dominant myeloid population reported in
NF1-associated peripheral nerve sheath tumours, and they are the dominant myeloid
population here too — 3,065 macrophages, 15.5% of all cells
([immunotyping study, PMID 39321200](https://pubmed.ncbi.nlm.nih.gov/39321200/);
[single-cell TME atlas, PMC12204358](https://pmc.ncbi.nlm.nih.gov/articles/PMC12204358/)).

### A case where marker evidence overrode a prior

The coarse M12 canonical-programme score labelled **C21** "B/plasma", because pDCs share
IGJ and MZB1 with plasma cells. The M14 markers (LILRA4, CLEC4C, SPIB, GZMB) are decisive
for plasmacytoid dendritic cells, and the annotation follows the markers. This is recorded
in the evidence table's `conflicting_evidence` column.

---

## 5. Deliberate non-assignments

Two clusters (720 cells, 3.7%) are labelled **Uncertain** rather than forced into a
compartment:

- **C12 (500 cells).** Carries *both* a hypoxia programme (SLC2A1, VEGFA, NDRG1, ADM,
  PHLDA3, ENO2) and the reported **perineurial** signature — perineurial cells are
  ITGB4+/SLC2A1(GLUT1)+, and GLUT1 is the standard perineurial IHC marker
  ([PMID 24719203](https://pubmed.ncbi.nlm.nih.gov/24719203/);
  [Modern Pathology](https://www.nature.com/articles/3880761)). It also shares PTHLH,
  TGFA and AP000462.2 with the malignant cluster C8. Malignant-hypoxic and perineurial are
  both defensible; distinguishing them needs CNV or GLUT1/EMA IHC. Left Uncertain.
- **C23 (220 cells).** Top markers are ribosomal pseudogenes and lncRNAs (RPL13AP5,
  RPS7P10, SNHG16) — a low-complexity/ambient-RNA signature, not a lineage.

**C15 (438 cells)** is labelled `Other → Low-quality / mitochondrial-high`: its top markers
are almost entirely mitochondrial transcripts, and it is 96% MPNST_4, matching the Phase 1
flag on that sample's `percent.mt`-correlated cluster (R² = 0.47). It is a technical
artefact and should be excluded from biological interpretation in Phase 3.

---

## 6. M12 caveats carried into the annotation

- **C1 (B/plasma).** M12 found B/plasma cohesion fell 44% under Harmony. C22 is explicitly
  caveated. C2 (plasma cells) was also flagged but is distributed across samples (43%
  MPNST_4), so it is less exposed.
- **C2 (fibroblast).** M12 found fibroblast cohesion fell 44%. C0, C3, C5 and C25 all carry
  the flag, and C5 and C25 are additionally single-sample dominated. Combined with §3.3,
  the fibroblast compartment is the least secure part of this annotation.

---

## 7. Limitations

1. **No CNV evidence.** Malignant calls rest on lineage markers, absence of alternative
   identity, and provenance. CNV inference would resolve §3.3 and is the top Phase 3
   priority.
2. **Dataset = patient = batch.** Eight of 26 clusters are >90% one patient. Patient-private
   expression is expected for malignant cells but cannot be distinguished from residual
   uncorrected batch effect.
3. **Annotation is per cluster**, so it inherits the resolution-1.0 boundaries. C19
   (venous + lymphatic endothelium) and C4 (T and NK together) are visibly under-resolved.
4. **Level 3 states are coarse** and were assigned only where markers clearly supported
   them; no exhaustion, IFN-response or polarisation sub-state was attempted.
5. **9.9% of cells carry Low confidence** and 3.7% are Uncertain. These were retained, not
   forced.
6. Literature was triangulated across MPNST spatial transcriptomics, MPNST single-cell
   multiomics, NF1 nerve sheath tumour immunotyping, peripheral-nerve fibroblast biology
   and canonical immunology — no single paper was copied.

---

## Sources

- [Spatially resolved transcriptomics of benign and malignant peripheral nerve sheath tumors — Neuro-Oncology, doi:10.1093/neuonc/noaf016](https://academic.oup.com/neuro-oncology/advance-article-abstract/doi/10.1093/neuonc/noaf016/7976860)
- [Single-cell multiomics identifies clinically relevant mesenchymal stem-like cells and key regulators for MPNST malignancy — Science Advances, doi:10.1126/sciadv.abo5442](https://www.science.org/doi/10.1126/sciadv.abo5442)
- [A Sequencing Overview of Malignant Peripheral Nerve Sheath Tumors — PMC11763529](https://pmc.ncbi.nlm.nih.gov/articles/PMC11763529/)
- [Single-cell tumor microenvironment profiling in NF1 nerve sheath tumors — PMC12204358](https://pmc.ncbi.nlm.nih.gov/articles/PMC12204358/)
- [Multidimensional Immunotyping of Human NF1-Associated Peripheral Nerve Sheath Tumors — PMID 39321200](https://pubmed.ncbi.nlm.nih.gov/39321200/)
- [Single-cell transcriptomic profiling of malignant peripheral nerve sheath tumors — Neuro-Oncology Advances](https://academic.oup.com/noa/article/8/Supplement_1/i43/8499782)
- [The fibroblast-derived protein PI16 controls neuropathic pain — PNAS, doi:10.1073/pnas.1913444117](https://www.pnas.org/doi/10.1073/pnas.1913444117)
- [Glut-1, best immunohistochemical marker for perineurial cells — PMID 24719203](https://pubmed.ncbi.nlm.nih.gov/24719203/)
- [Immunohistochemical Demonstration of EMA/Glut1-Positive Perineurial Cells — Modern Pathology](https://www.nature.com/articles/3880761)
- [Glial-to-mesenchymal transition of tumor Schwann cells in MPNST — Science Advances, doi:10.1126/sciadv.adt9210](https://www.science.org/doi/10.1126/sciadv.adt9210)
