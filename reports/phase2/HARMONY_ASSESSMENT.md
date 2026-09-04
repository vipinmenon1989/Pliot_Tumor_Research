# Harmony Integration Assessment

**Phase 2 · Milestone M12 · MPNST single-cell integration**

*Generated: 2026-09-02 · SLURM JobID 19886628*

**Recommendation: `ACCEPT DEFAULT HARMONY WITH CAVEATS`** — see §J.

---

## A. Integration setup

| Item | Value |
| --- | --- |
| Evaluated object | `results/phase2/harmony/phase2_harmony_integrated.rds` |
| md5 / sha256 (re-verified at runtime) | `cf63e84313a91de25fe2e41660f78f06` / `6085976f78…e69d` ✅ both match M11 |
| Cells evaluated | **19,716 — all of them. No subsampling anywhere in M12.** |
| PRE space | reduction `pca`, dims 1:30 (frozen Phase 1 M8 baseline, assay `SCT`) |
| POST space | reduction `postint_harmony`, dims 1:30 (M11 output) |
| Harmony grouping variable | `sample_id` (4 levels) |
| Harmony version | 1.2.4; all scientific parameters at package defaults (`theta` 2, `sigma` 0.1, `lambda` 1, `nclust` 100, `max_iter` 10, `early_stop` TRUE) |
| Harmony re-run or re-tuned in M12? | **No.** The M11 default result is evaluated exactly as produced. |
| Random seed | 42 throughout |

### Fairness of the comparison

Both spaces use the **same 19,716 cells, the same 30 dimensions, the same metadata and
the same seed**. The Phase 1 UMAP `umap_preintegration` was reused unmodified; the
Harmony UMAP `umap_harmony_m12` was built with an identical `RunUMAP` call in which only
the input reduction changed (dims 1:30, `seed.use = 42`, `n.neighbors = 30`,
`min.dist = 0.3`, `metric = "cosine"` — all Seurat defaults). **No UMAP or neighbour
parameter was optimised to make either side look better.** Digests of the Phase 1 PCA,
the Phase 1 UMAP and the M11 Harmony embedding were re-checked after all computation and
were unchanged.

### Metric implementations

`lisi`, `kBET` and `clustree` are not installed and were **not** installed. Substitutes,
all computed natively on exact k-nearest neighbours (`RANN::nn2`) and full pairwise
distances (`cluster::silhouette`):

| Unavailable | Substitute used |
| --- | --- |
| iLISI | Inverse Simpson index `1 / Σp²` over `sample_id` proportions among the exact k nearest neighbours. Same estimator as `lisi`, computed on exact kNN rather than a Gaussian-kernel neighbourhood. |
| kBET | Per-cell **dominance ratio** = observed same-sample neighbour fraction ÷ that sample's global abundance, reported descriptively against the value 1.0. |
| cLISI | Two non-legacy biological references (§E). |

---

## B. Metadata structure

| Concept | Field | Levels | Note |
| --- | --- | ---: | --- |
| Dataset | `sample_id` | 4 | |
| Sample | `sample_id` | 4 | identical to dataset |
| Patient | — | — | 1 sample = 1 patient; no separate field exists |
| Technical batch | — | — | never recorded |
| Biological condition | — | — | **does not exist** — no treatment, stage, site or NF1 status |
| Harmony grouping variable | `sample_id` | 4 | identical to `orig.ident` cell-for-cell |

**Dataset = sample = patient = presumed technical batch is a single 4-level variable.**
Figures A, B, C and E of the M12 specification therefore collapse into one comparison;
this is documented on the figure itself rather than padded out with four identical
copies. Figure D (biological condition) cannot be produced because no such field exists,
and none was fabricated.

Global composition: MPNST_1 38.62% · MPNST_2 11.58% · MPNST_3 14.91% · MPNST_4 34.88%.

### Reference values for a "fully mixed" embedding

Because the four samples are unequal in size, perfect mixing does **not** mean maximum
entropy. If every neighbourhood matched the global composition exactly:

| Metric | Fully-mixed reference | Interpretation |
| --- | ---: | --- |
| Same-sample neighbour fraction | **0.3065** (= Σp²) | not 0 |
| Neighbourhood entropy | **1.2683 nats** | not log K = 1.3863 |
| Inverse Simpson | **3.2627** | not 4 |
| Dominance ratio | **1.00** | |

Every statement below is made against these composition-adjusted references, not against
the naive maxima.

### Reconciliation with the frozen Phase 1 baseline

The M12 PRE-Harmony metrics reproduce Phase 1 M8 **exactly**, which validates the
pipeline:

| Sample | Phase 1 M8 same-dataset fraction | M12 PRE same-sample fraction |
| --- | ---: | ---: |
| MPNST_1 | 0.988321295688334 | 0.988321295688334 |
| MPNST_2 | 0.910245183887916 | 0.910245183887916 |
| MPNST_3 | 0.969024943310658 | 0.969024943310658 |
| MPNST_4 | 0.977907033105521 | 0.977907033105521 |

Identical to 15 significant figures. The two headline numbers differ only in
presentation, not in substance:

- Phase 1 quoted **0.9613** = the *unweighted* mean of the four per-sample means. The
  *cell-weighted* mean of the same four numbers is **0.9728**, which is what M12 reports.
- Phase 1 quoted entropy **0.0966** in **log2 (bits)**. M12 reports natural log. Phase 1's
  cell-weighted entropy is 0.0707 bits = 0.0707 × ln2 = **0.0490 nats**, exactly the M12
  PRE value.

---

## C. Visual evidence

Figures in `reports/phase2/figures/m12/` (PDF + PNG). Principal observations:

- **`01_pre_post_umap_by_sample_id`** — pre-Harmony, the four samples occupy essentially
  disjoint territories. Post-Harmony a large shared central region appears in which all
  four samples interleave, while several substantial lobes remain visibly dominated by a
  single sample.
- **`02_pre_post_umap_faceted_by_sample_id`** — each sample's occupancy of the shared
  embedding. Post-Harmony, all four samples contribute to the same central territory;
  MPNST_1 additionally retains a large private lobe, and MPNST_3 a smaller one.
- **`04_pre_post_umap_by_phase1_independent_clusters`** — the 54 independently derived
  Phase 1 per-sample clusters remain visible as compact, contiguous islands after
  Harmony. They are not smeared out.
- **`11_pre_post_umap_by_canonical_program`** — post-Harmony, cells scoring highest for
  the T/NK, myeloid, endothelial and mural programs converge into shared, cross-sample
  territories. Cells scoring highest for the Schwann/neural-crest program remain largely
  within one region dominated by MPNST_1.
- **`03_pre_post_umap_by_qc_covariates`** — no QC covariate (`nCount_RNA`,
  `nFeature_RNA`, `percent.mt`, `percent.ribo`) tracks the major axes of either
  embedding, so the residual structure is not a library-depth artefact.

**Visual impressions are recorded here as observations only.** No conclusion in §J rests
on them; the quantitative results in §D and §E carry the argument.

---

## D. Quantitative technical mixing

All values are cell-weighted means over all 19,716 cells. Machine-readable:
`results/phase2/harmony/evaluation/pre_post_mixing_summary.tsv`.

### D.1 Global

| Metric (k = 15) | PRE | POST | Fully-mixed reference | Gap closed | Desired direction |
| --- | ---: | ---: | ---: | ---: | --- |
| Same-sample neighbour fraction | 0.9728 | **0.7784** | 0.3065 | **29.2%** | lower |
| Neighbourhood entropy (nats) | 0.0490 | **0.3565** | 1.2683 | **25.2%** | higher |
| Inverse Simpson (iLISI-equivalent) | 1.0549 | **1.4408** | 3.2627 | **17.5%** | higher |
| Dominance ratio | 3.30 | **2.53** | 1.00 | — | toward 1.0 |

| Metric (k = 50, robustness) | PRE | POST | Gap closed |
| --- | ---: | ---: | ---: |
| Same-sample neighbour fraction | 0.9443 | 0.7053 | 37.5% |
| Neighbourhood entropy (nats) | 0.1098 | 0.5207 | 35.5% |

The direction is unambiguous and consistent at both neighbourhood sizes: local
neighbourhoods became substantially more sample-diverse. In absolute terms the embedding
remains far from composition-matched — 78% of a cell's 15 nearest neighbours still come
from its own sample. §D.3 shows that this residual is not uniform, and §F/§G explain why
the uniform reading would be the wrong one.

### D.2 Per sample

| Sample | n | Global share | Same-sample fraction PRE → POST | Entropy PRE → POST | Dominance ratio PRE → POST |
| --- | ---: | ---: | --- | --- | --- |
| MPNST_1 | 7,615 | 38.6% | 0.9883 → **0.8490** | 0.0227 → 0.2647 | 2.56 → 2.20 |
| MPNST_2 | 2,284 | 11.6% | 0.9102 → **0.5037** | 0.1496 → 0.6458 | 7.86 → 4.35 |
| MPNST_3 | 2,940 | 14.9% | 0.9690 → **0.7383** | 0.0521 → 0.4264 | 6.50 → 4.95 |
| MPNST_4 | 6,877 | 34.9% | 0.9779 → **0.8086** | 0.0434 → 0.3324 | 2.80 → 2.32 |

Every sample moved in the desired direction. MPNST_2 mixed the most (same-sample
fraction more than halved); MPNST_1 the least. The two small samples retain the highest
dominance ratios, consistent with the 3.3× size imbalance and the uniform default
`theta = 2` noted in M11.

### D.3 Technical separation — silhouette width, and why it is uninformative here

| Label | PRE | POST | Change |
| --- | ---: | ---: | ---: |
| `sample_id` (all cells) | **−0.0104** | **−0.0076** | +0.0028 |

**This metric failed to detect the batch structure and must not be read as evidence.**
The pre-Harmony space has 97% same-sample nearest neighbours yet a `sample_id`
silhouette of essentially zero. The reason is structural: silhouette compares each cell's
mean distance to its own group against its mean distance to other groups, and each
sample here spans the full range of malignant, immune, stromal and vascular states. The
mean within-sample distance is therefore about as large as the mean between-sample
distance, and the global statistic cancels out. Batch structure in this dataset is
**local**, and only the local kNN metrics in §D.1–D.2 capture it.

Reported for completeness and as a methodological caveat, not as support for any
conclusion. The same silhouette machinery *is* informative for the compact biological
program groups in §E.

### D.4 Compartment-resolved mixing — the decisive result

Same-sample neighbour fraction (k = 15) broken down by canonical program (§E.2).
Source: `sample_restricted_population_check.tsv`.

| Program | n | Dominant sample (share) | Same-sample fraction PRE → POST | Δ |
| --- | ---: | --- | --- | ---: |
| Panleukocyte (PTPRC) | 1,396 | MPNST_1 (45.4%) | 0.933 → **0.527** | **−0.406** |
| Myeloid | 2,441 | MPNST_4 (34.9%) | 0.958 → **0.611** | **−0.347** |
| T/NK | 1,093 | MPNST_3 (53.2%) | 0.953 → **0.689** | **−0.264** |
| Endothelial | 966 | MPNST_1 (39.3%) | 0.962 → **0.709** | **−0.253** |
| Mural | 932 | MPNST_1 (42.7%) | 0.951 → **0.724** | **−0.227** |
| Mast | 40 | MPNST_3 (75.0%) | 0.885 → 0.678 | −0.207 |
| Proliferation | 668 | MPNST_4 (64.4%) | 0.944 → 0.786 | −0.158 |
| Fibroblast | 4,075 | MPNST_4 (70.6%) | 0.990 → 0.839 | −0.151 |
| **Schwann/neural-crest** | **2,087** | **MPNST_1 (86.5%)** | **0.990 → 0.900** | **−0.090** |
| B/plasma | 836 | MPNST_3 (93.5%) | 0.986 → 0.918 | −0.067 |

The residual sample structure is **not uniform**. Harmony mixed the shared
immune and vascular compartments hard (−0.23 to −0.41) and left the
Schwann/neural-crest compartment — the presumptive malignant lineage, 86.5% confined to
MPNST_1 — almost as sample-private as it began (−0.09). The B/plasma population, 93.5%
confined to a single sample and therefore with almost no counterpart to mix into, also
stayed private (−0.07).

This is the behaviour one would want. It also explains the modest global figure in §D.1:
the cell-weighted average is dominated by the large Schwann-lineage, fibroblast and
unassigned populations for which patient-private structure is expected, not by a failure
to align the shared compartments.

---

## E. Biological preservation

Machine-readable: `biological_preservation_summary.tsv`, `technical_silhouette_summary.tsv`.

### E.0 What was deliberately NOT used

Legacy `orig.anno` was **excluded**. It derives from the original integrated analysis;
M12 forbids validating Harmony with legacy integration-derived annotation, and the M10
gate left its use as a reference label unanswered with a default of "no". Legacy
clusters and legacy Harmony/CCA/RPCA/MNN embeddings were likewise not used. Two
independent, non-legacy references were used instead.

### E.1 Reference 1 — independent Phase 1 per-sample structure

Phase 1 clustered each dataset **separately, before any integration existed**. Those 54
labels (`preint_recommended_cluster`) are an integration-free description of real
within-sample transcriptional structure.

| Measure | PRE | POST | Relative change |
| --- | ---: | ---: | ---: |
| Within-sample kNN retention (k = 15) | — | **0.8206** | — |
| Phase 1 cluster coherence | 0.8644 | **0.8427** | **−2.5%** |
| Global kNN retention | — | 0.6598 | — |

Per sample:

| Sample | Within-sample kNN retention | Phase 1 cluster coherence PRE → POST |
| --- | ---: | --- |
| MPNST_1 | 0.8186 | 0.9202 → 0.8983 (−2.4%) |
| MPNST_2 | 0.8755 | 0.9465 → 0.9381 (−0.9%) |
| MPNST_3 | 0.8759 | 0.9150 → 0.8959 (−2.1%) |
| MPNST_4 | 0.7810 | 0.7538 → 0.7266 (−3.6%) |

**Interpretation.** Global kNN retention fell to 0.66 — that fall *is* the correction and
is not damage. The discriminating number is the within-sample one: **82% of each cell's
own-sample nearest neighbours survived Harmony**, and agreement with independently
derived pre-integration clusters fell by only 2.5% relative. Within-sample biology was
largely left alone. MPNST_4 shows the largest loss on both measures, consistent with the
Phase 1 finding that its cluster C09 correlates with `percent.mt` (R² = 0.47).

### E.2 Reference 2 — canonical broad-lineage program scores

Ten canonical broad-lineage gene sets were scored per cell with `AddModuleScore`
(`nbin = 24`, `ctrl = 100`, both package defaults; seed 42) on the default `SCT` assay.
Genes found are recorded in `biological_program_gene_sets.tsv`:

| Program | Genes used |
| --- | --- |
| Panleukocyte | PTPRC |
| T_NK | CD3D CD3E CD2 IL7R TRAC NKG7 GNLY KLRD1 |
| Myeloid | LYZ CD68 AIF1 CSF1R ITGAM C1QA C1QB FCGR3A |
| B_Plasma | 4/5 (CD79B absent from the assay) |
| Mast | 3/4 |
| Endothelial | 3/4 |
| Fibroblast | COL1A1 COL1A2 COL3A1 DCN LUM PDGFRA |
| Mural | ACTA2 RGS5 PDGFRB MYH11 NOTCH3 |
| SchwannNC | SOX10 S100B PLP1 MPZ NGFR PMP22 ERBB3 |
| Proliferation | MKI67 TOP2A CCNB1 CDK1 PCNA |

*Rationale and provenance:* these are long-established, textbook-level lineage markers
used here purely to ask whether cells sharing a broad transcriptional program stay
together. **They are not a cell-type annotation** — annotation is M15, requires
literature triangulation, and is not performed here. Each cell was given the program with
the highest z-scored value, or `Unassigned` when no program reached z ≥ 1 (5,182 cells,
26.3%). The labels are identical for the PRE and POST comparison, so any normalisation
quirk affects both sides equally and cannot bias the contrast. They were never written
into any saved object.

Silhouette width of the program label:

| Program | n | PRE | POST | Change |
| --- | ---: | ---: | ---: | ---: |
| **T/NK** | 1,093 | 0.252 | **0.461** | **+0.209** |
| **Endothelial** | 966 | 0.070 | **0.288** | **+0.218** |
| **Schwann/neural-crest** | 2,087 | 0.201 | **0.372** | **+0.170** |
| **Myeloid** | 2,441 | 0.045 | **0.149** | **+0.104** |
| Mural | 932 | −0.123 | −0.083 | +0.040 |
| Panleukocyte | 1,396 | −0.143 | −0.152 | −0.008 |
| Proliferation | 668 | −0.109 | −0.151 | −0.042 |
| **Fibroblast** | 4,075 | 0.174 | **0.097** | **−0.077** |
| **B/plasma** | 836 | 0.437 | **0.243** | **−0.194** |
| Mast | 40 | −0.160 | −0.296 | −0.136 |
| **Overall** | 19,716 | **0.0627** | **0.0729** | **+0.0103** |

**Interpretation.** Overall program cohesion did not fall — it rose slightly. Six of ten
programs became *more* cohesive, and the four largest gains are in exactly the
compartments that exist in several tumours (T/NK, endothelial, myeloid) plus the
Schwann-lineage compartment. Cells of the same broad lineage from different patients were
brought together, which is what integration is for. Two programs lost cohesion
materially (§F).

---

## F. Evidence for potential over-correction

Real, localised, and limited to two populations:

1. **B/plasma cohesion fell 44%** (silhouette 0.437 → 0.243, n = 836). This was the most
   cohesive program in the pre-Harmony space and is 93.5% confined to MPNST_3. Its
   same-sample neighbour fraction barely moved (0.986 → 0.918), so these cells were not
   dispersed *across* samples; they became less separable from adjacent lymphoid and
   myeloid territory *within* the embedding. Because there is essentially no B-cell
   population in the other three tumours, there was nothing legitimate for Harmony to
   align them with, and the diversity penalty acted on them anyway. This is the textbook
   over-correction failure mode for a batch-private population, and it is present here at
   moderate severity.
2. **Fibroblast cohesion fell 44%** (silhouette 0.174 → 0.097, n = 4,075; 70.6% MPNST_4).
   Same pattern, weaker starting cohesion. Its within-sample kNN retention is also the
   lowest of any program (0.795).
3. **Mast cells** (n = 40) show the largest relative silhouette drop, but at n = 40 the
   estimate is unstable and no weight is placed on it.
4. **MPNST_4** shows the largest within-sample damage of the four samples (retention
   0.781, Phase 1 cluster coherence −3.6%).

Evidence **against** widespread over-correction:
- Overall program silhouette *rose* (+0.0103); six of ten programs improved.
- Within-sample kNN retention is 0.82 and Phase 1 cluster coherence fell only 2.5%.
- The Schwann/neural-crest compartment — the population whose erasure would be most
  damaging — both stayed sample-private (0.990 → 0.900) and became *more* internally
  cohesive (0.201 → 0.372). It was not collapsed into the other tumours.
- The 54 independent Phase 1 clusters remain visible as discrete islands
  (`04_pre_post_umap_by_phase1_independent_clusters`).

**Conclusion: over-correction is present but bounded**, affecting the B/plasma and
fibroblast programs rather than the embedding as a whole.

---

## G. Evidence for potential under-correction

1. **In absolute terms the embedding is far from composition-matched.** At k = 15 the
   same-sample neighbour fraction is 0.7784 against a fully-mixed reference of 0.3065 —
   only 29.2% of the achievable gap was closed. Entropy closed 25.2%, inverse Simpson
   17.5%.
2. **Dominance ratios remain well above 1.0 for every sample** (2.20 / 4.35 / 4.95 /
   2.32), and highest for the two small samples, so the default uniform `theta = 2`
   corrected the minority samples proportionally less than the majority ones.
3. **Sample-dominated regions persist**, chiefly the large MPNST_1 lobe visible in
   `02_pre_post_umap_faceted_by_sample_id`.

**Why this is largely appropriate rather than a failure.** §D.4 shows the residual is
concentrated in the Schwann/neural-crest (0.900), B/plasma (0.918) and fibroblast (0.839)
compartments, while the shared immune and vascular compartments reached 0.53–0.72. In a
dataset where dataset = patient and roughly half the cells are tumour, patient-private
malignant structure is *expected to survive integration* — forcing it to mix would be the
error, not the fix. The persistent MPNST_1 lobe coincides with the compartment holding
86.5% of the Schwann/neural-crest signal.

**Conclusion: the shared compartments are adequately corrected. The apparent global
under-correction is dominated by structure that should not be corrected.** The one
genuine concern is the weaker relative correction of the small samples (MPNST_2,
MPNST_3), which follows from uniform `theta` under a 3.3× size imbalance.

---

## H. Confounding

`sample_id` is simultaneously the dataset, the patient and the only available proxy for
technical batch. No clinical or technical covariate of any kind exists in the object.
Consequences that cannot be argued away:

1. **Integration quality cannot be cleanly separated from biological preservation.** Every
   unit of residual `sample_id` structure is some unknown mixture of uncorrected batch
   effect and genuine between-patient tumour biology. There is no design, covariate or
   statistic in this dataset that can apportion it.
2. **No amount of metric improvement proves the correction was "right".** M12 can show
   that neighbourhoods became more sample-diverse while independent within-sample
   structure survived. It cannot show that the specific variance Harmony removed was
   technical.
3. **The compartment-resolved result in §D.4 is the strongest available evidence**
   precisely because it is *differential*: Harmony behaved differently in compartments
   expected to be shared versus compartments expected to be patient-private. That pattern
   is hard to produce by accident, but it is a consistency argument, not proof.
4. **Downstream constraint, to be restated in every later Phase 2 document:** the
   integrated embedding cannot support any between-tumour, between-patient,
   between-condition or differential-abundance claim. The between-tumour axis is the axis
   that was deliberately removed.

---

## I. Limitations

1. `lisi` and `kBET` are not installed; native inverse-Simpson and dominance-ratio
   substitutes were used. They are close analogues but not the published implementations.
2. Global silhouette of `sample_id` proved insensitive to this dataset's local batch
   structure (§D.3) and contributes nothing to the conclusion.
3. The canonical-program labels are a coarse sanity check: 26.3% of cells are
   `Unassigned`, `Panleukocyte` rests on a single gene, `Mast` has n = 40, and the
   `SCT` `data` layer spans four SCT models (M10 Finding 1). They are used only for the
   internally fair PRE/POST contrast.
4. No statistical test is reported. With 19,716 cells any difference would be
   "significant"; cells are not independent biological replicates, and effect sizes and
   distributions are reported instead. Per-sample summaries are given wherever the
   quantity supports them.
5. Only one integration method and one parameter set were evaluated. CCA, RPCA and MNN
   (proposed in Phase 1 `INTEGRATION_PREPARATION.md`) were not benchmarked; that was not
   in scope.
6. Some canonical genes were absent from the assay (CD79B, and one gene each from the
   Mast and Endothelial sets), slightly weakening those scores.
7. The UMAPs are used for description only. No conclusion depends on them.

---

## J. Recommendation

# `ACCEPT DEFAULT HARMONY WITH CAVEATS`

**Both required criteria are met.**

*Technical mixing improved substantially and in the right places.* Same-sample nearest
neighbours fell 0.9728 → 0.7784 (k = 15) and 0.9443 → 0.7053 (k = 50); entropy rose
7.3-fold; inverse Simpson rose 1.05 → 1.44. Critically, the improvement is
compartment-specific: the shared immune and vascular compartments reached 0.53–0.72
same-sample neighbours (from 0.93–0.96), which is the alignment needed for joint
cell-type annotation.

*Biological structure was preserved.* Within-sample kNN retention is 0.82; agreement with
the 54 independently derived, integration-free Phase 1 clusters fell only 2.5% relative;
overall canonical-program cohesion rose slightly (+0.0103), with the four largest gains in
T/NK, endothelial, Schwann/neural-crest and myeloid. The Schwann/neural-crest compartment
was neither dispersed across patients nor collapsed.

**Why not `ACCEPT` outright** — two caveats must be carried forward:

- **C1. B/plasma (n = 836, 93.5% MPNST_3) lost 44% of its cohesion.** Treat any B/plasma
  cluster arising in M13–M15 as potentially distorted, and cross-check it against the
  preserved non-integrated baseline before annotating it.
- **C2. Fibroblast (n = 4,075, 70.6% MPNST_4) lost 44% of its cohesion**, and MPNST_4 has
  the weakest within-sample preservation overall (retention 0.781). Apply the same
  cross-check to fibroblast/stromal clusters.

**Why not `REQUIRES SENSITIVITY ANALYSIS`.** The obvious knob is `theta`, and the
evidence argues against turning it up. Raising `theta` would push harder on exactly the
Schwann/neural-crest and B/plasma populations that are already the most fragile, in a
dataset where patient-private malignant structure is expected to survive. Lowering it
would undo the immune/vascular alignment that is the whole purpose of integrating. The
default result already sits in the desired regime.

**If the researcher nonetheless wants a sensitivity analysis**, the two experiments worth
running — neither executed, both requiring explicit approval — are:

1. **`theta` sweep {1, 2, 4}** with everything else fixed, re-running this exact M12
   metric suite on each, to map the trade-off curve between the immune/vascular mixing
   gain and the B/plasma and fibroblast cohesion loss.
2. **Per-variable `lambda` or a `theta` weighting that accounts for the 3.3× sample-size
   imbalance**, targeting the under-correction of MPNST_2 and MPNST_3 identified in §G.2
   without increasing pressure on the malignant compartment.

**Conditional on C1 and C2 being carried into M13–M15, the M11 default Harmony embedding
is fit for post-integration clustering, marker discovery and joint annotation.**

---

## Files

Tables — `results/phase2/harmony/evaluation/`:
`pre_post_mixing_summary.tsv` · `neighborhood_mixing_metrics_by_sample.tsv` ·
`technical_silhouette_summary.tsv` · `biological_preservation_summary.tsv` ·
`sample_restricted_population_check.tsv` · `biological_program_gene_sets.tsv` ·
`program_by_sample_composition.tsv` · `harmony_evaluation_parameters.tsv` ·
`m12_headline_metrics.json` · `prov_m12_evaluation.json` ·
`umap_harmony_m12_embedding.{tsv,rds}` · `figure_index_m12.tsv`

Figures — `reports/phase2/figures/m12/` (12 figures, PDF + PNG), all registered in
`reports/FIGURE_INDEX.tsv`.
