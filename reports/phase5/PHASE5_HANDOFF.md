# MPNST Phase 5 Handoff — Continuous Malignant Transcriptional Programs

*Milestones M36–M41 · executed and frozen 2026-09-03/04 · readable without conversation history*

> **A cNMF program is a continuous pattern of co-regulated expression.** It is **not** a cell type,
> a lineage, a discrete state or a clone, and those four remain separate concepts throughout.
> **n = 4 patients** and `sample_id` = patient = dataset, so cells are **not** biological replicates
> and no cell-level p-value carries inferential weight here.
> **Phases 1–4 are unmodified.** The Phase 4 object was opened read-only and its md5
> `e85ba8486e456917e2483f2773bdbaf3` was re-verified after every milestone.

---

## 1. The question, and why it is not the Phase 4 question again

Phase 4 ended on a negative result: **0 of 8 discrete malignant transcriptional states were
recurrent** across ≥ 3 patients, and **97.2% of malignant cells sat in patient-private states**.
Phase 5 does not re-run clustering to rescue that. It asks the *continuous* version, which is a
different object — a cell belongs to exactly one cluster, but carries a loading on every program,
and two patients can share a program while their cells cluster apart.

> **Do recurrent continuous malignant transcriptional programs exist across patients even when
> discrete malignant clusters do not?**

**Answer: no, not in this cohort.** The result is reported as the negative it is.

## 2. Input and the field that was verified rather than assumed

`results/phase4/phase4_final_object.rds`, verified by md5 **and** sha256 before loading. 19,716
cells; 13 frozen Phase 4 counts re-asserted and all passed.

M39 compares historical fibroblasts Phase 4 called malignant against historical fibroblasts that
stayed non-malignant, so the grouping variable had to be the *pre-Phase-4* annotation.
`annotation_ccc_refined` is **barred** from that role — it holds 908 fibroblasts, all
Non-malignant, because it already encodes the Phase 4 conclusion. Every `annotation*` column was
tested against the frozen split 5,064 → 4,036 / 908 / 120; `annotation_ccc_phase3` and
`annotation_ccc` both reproduce it, and **`annotation_ccc_phase3` was selected because it
qualified**, not because it was expected.

## 3. Patient imbalance, which shaped the whole design

| patient | malignant cells | share | median nCount | median nFeature | High / Moderate | clone-reliable |
| --- | ---: | ---: | ---: | ---: | ---: | :--: |
| MPNST_4 | 3,685 | **57.3%** | 5,250 | 2,249 | 936 / 2,749 | yes |
| MPNST_1 | 1,886 | 29.3% | 16,034 | 5,361 | 1,861 / 25 | yes |
| MPNST_2 | 651 | 10.1% | 18,901 | 3,978 | 464 / 187 | yes |
| MPNST_3 | 212 | 3.3% | 9,327 | 3,299 | **0** / 212 | no |

**17.4× imbalance, confounded with depth** — the largest patient is also the shallowest. This is
why the patient-balanced analysis is part of the design rather than an afterthought.

## 4. Method

cNMF 1.7.1 in the **isolated `p5_cnmf_env`**. `R_env` was not modified: `conda env export` before
and after Phase 5 are byte-identical. MSigDB v2024.1.Hs gene sets were downloaded as versioned GMTs
with recorded md5s rather than installing `msigdbr` into the frozen stack.

Declared **before any factor was inspected**: malignant cells only (6,434); RNA raw integer counts
(never Harmony, UMAP, PCA scores, SCT residuals or the SCEVAN CNA matrix); genes detected in ≥ 0.5%
of malignant cells with `^MT-` removed (19,663); K grid 4–15; 100 replicates per K; seed 42;
`numgenes` 2000; consensus at local-density-threshold 0.10.

**Cell cycle, ECM, HLA, interferon, Schwann, neural-crest and angiogenesis genes were all retained**
in the primary run — removing them would delete the signal the phase exists to find.

**K-selection rule, declared in code before it ran**: the largest K with (i) max pairwise cosine
between consensus spectra ≤ 0.75, (ii) every program dominant in ≥ 1% of cells, (iii) stability ≥
the median of the K values satisfying both. It references only measured properties of the
factorization — no gene, no pathway, no label. Eligible K = {4,5,6,7,8}; **K = 8**.

## 5. The eight programs

| id | label | top genes | top gene set | technical? |
| --- | --- | --- | --- | :--: |
| P1 | Translation_ribosomal | RPL10 RPL19 RPS2 RPLP1 EEF1A1 | REACTOME_EUKARYOTIC_TRANSLATION_INITIATION | **yes** (54% ribosomal) |
| P2 | Unassigned_NEURONAL_SYSTEM | TSHZ2 AUTS2 TENM2 AFF3 CHN1 NLGN1 | REACTOME_NEURONAL_SYSTEM | no |
| P3 | **Mesenchymal_ECM_mixed** | SERPINF1 CST3 FBLN1 SFRP2 EFEMP1 MFAP5 | HALLMARK_TNFA_SIGNALING_VIA_NFKB | no |
| P4 | Hypoxia_Angio | GBE1 IGF1R VEGFA ADAMTS9-AS2 | HALLMARK_UV_RESPONSE_DN | no |
| P5 | Translation_ribosomal | RPS7P1 RPL3P4 (pseudogenes) | REACTOME_TRANSLATION | **yes** (66% pseudogene) |
| P6 | Schwann_like_mixed | S100B GPM6B SERPINE2 CRYAB MIA | HALLMARK_INTERFERON_GAMMA_RESPONSE | no |
| P7 | **Cycling** | NUSAP1 PRC1 TOP2A UBE2C MKI67 | REACTOME_CELL_CYCLE_MITOTIC | no |
| P8 | Myeloid_ambient_like | C1QC C1QB C1QA CD14 TYROBP CSF1R | REACTOME_NEUTROPHIL_DEGRANULATION | **yes** (18 myeloid markers) |

**Three of eight are technical-dominated and are named for it.** The first labelling pass called P1
`Mesenchymal_ECM` on a 0.20 marker-family overlap while its own top 15 genes were ribosomal
proteins; that label would have misled every downstream reader, so an explicit technical-content
check was added and M38 re-run. Only the *naming* changed — not the factorization, K, or the
recurrence rule.

## 6. Recurrence — the headline

Criteria declared before any program was labelled, deliberately the same shape as Phase 4's rule:
program-active at ≥ 0.20 relative usage; a patient *carries* a program with ≥ 10 active cells **and**
≥ 5% of its malignant cells; recurrent needs ≥ 3 of 4 patients.

**0 recurrent · 1 shared-limited · 7 patient-private.** The single shared-limited program is P1
(MPNST_2 + MPNST_4) — one of the technical ones. **Not one biologically interpretable program is
carried by even two patients at the declared threshold.**

Dominant-program assignment is equally stark: P1 → 3,490 MPNST_4; P2 → 1,002 MPNST_1; P3 → 651
MPNST_2 (every malignant cell in that patient); P4 → 581 MPNST_1; P5 → 301 MPNST_1; P6 → 122
MPNST_3 + 55 MPNST_4; P7 → 139 MPNST_4 + 3 MPNST_3 + 2 MPNST_1; P8 → 83 MPNST_3.

**Threshold sensitivity, across a grid rather than at one cut.** At a relaxed 0.10 activity
threshold, two programs reach three patients: P1 (technical) and **P3, the ECM/mesenchymal
program** — per-patient active fractions MPNST_2 1.000, MPNST_4 0.108, MPNST_3 0.085, MPNST_1 0.007.
**P3 is the one biologically interpretable program that comes close to recurrence**, and the honest
statement is exactly that: recurrent at 10%, not at the pre-declared 20%.

## 7. The five Mesenchymal_ECM states are not one shared program

| Phase 4 state | patient | dominant Phase 5 program | median relative usage |
| --- | --- | --- | ---: |
| Mesenchymal_ECM-1 | MPNST_1 (1,851) | P2 Neuronal | 0.460 |
| Mesenchymal_ECM-2 | MPNST_4 (1,587) | **P1** | 0.867 |
| Mesenchymal_ECM-3 | MPNST_4 (1,291) | **P1** | 0.894 |
| Mesenchymal_ECM-4 | MPNST_2 (645) | P3 Mesenchymal_ECM | 0.887 |
| Mesenchymal_ECM-5 | MPNST_4 (505) | **P1** | 0.877 |

The structure is informative in **both** directions. *Within* MPNST_4, three separately-clustered ECM
states (3,383 cells) collapse onto **one** continuous program — Phase 4's discrete partition was
finer there than the continuous structure warrants. *Across* patients, the five states map to three
different programs. **This is possible result C of the three the phase brief anticipated: ECM
programs remain substantially patient-private.**

## 8. Malignant ECM-like cells versus true fibroblasts

Groups from the verified historical field: A = 4,036 malignant, B = 908 non-malignant, C = 120
Ambiguous (projection only). **Evaluability was checked before any test was run**, with a declared
minimum of 30 cells on both sides:

| patient | A | B | evaluable | reason |
| --- | ---: | ---: | :--: | --- |
| MPNST_1 | 512 | 303 | **yes** | — |
| MPNST_2 | 637 | 503 | **yes** | — |
| MPNST_3 | 0 | 101 | no | Phase 4 disabled its malignant promotions |
| MPNST_4 | 2,887 | **1** | no | too few non-malignant fibroblasts to compare against |

**A four-patient paired test was not constructed** — the groups do not coexist in 2 of 4 patients.

**The cross-patient Spearman correlation of the pseudobulk log2FC is only 0.091.** Most of the
A-vs-B difference is patient-specific. That number bounds the whole comparison and is stated first
rather than buried.

**Signature** (declared criteria applied in *every* evaluable patient: |pseudobulk log2FC| ≥ 1,
cell-level AUC ≥ 0.65, detected in ≥ 10% of the higher group; **no classifier, no train/test
split**):

* **61 up in malignant ECM-like** — CA12, SCG2, IGFBP3, TMEM176A/B, TUBB2B, COL14A1, COL11A1, GJA1,
  IGFBP2, SPOCK1, PTGIS, BASP1, ACKR3, ANGPT1, ENPP2 …
* **31 up in true fibroblast** — CDH19, APOD, SCN7A, ABCA6, ABCA8, ABCA9, ABCA10, VIT, PAMR1,
  SPARCL1, GPC3, A2M, LAMA2, CFH, NID1 …

The fibroblast side is coherent in a way worth naming: **CDH19, SCN7A, APOD and the ABCA6/8/9/10
cluster mark nerve-associated / endoneurial fibroblasts and non-myelinating Schwann-lineage
stroma** — exactly the resident population a peripheral-nerve-sheath tumour would retain.

The 120 Ambiguous fibroblasts were projected onto the fixed spectra for description and **were not
reclassified**; an assertion confirms every one still carries `malignancy_refined == "Ambiguous"`.

## 9. Robustness

Full detail in `ROBUSTNESS_REPORT.md`. The two results that matter most:

* **Removing the technical gene classes changes nothing.** A fourth full cNMF run with ribosomal,
  pseudogene/lncRNA and canonical myeloid genes removed, across the entire K grid, gives **0
  recurrent programs at every K from 5 to 15**. The patient-private result is not an artefact of
  those genes.
* **Every patient-private program vanishes when its own patient is withheld** (LOO cosine 0.17–0.35),
  while P6 and P7 — the two whose dominant patient is not the one that defines them — survive at
  0.80 and 0.88.

Programs are stable to rank (median cosine 0.998 at K ± 1), patient balance (median 0.940, minimum
0.504, all supported), cell-cycle removal (median 1.000) and confidence filtering, and no program's
usage tracks a within-patient technical covariate above |ρ| = 0.37.

## 10. Final object

```text
path      results/phase5/phase5_final_object.rds
size      6,058,471,006 bytes
md5       839e5157bc7c3470ddf86746c2e719e1
sha256    fe99ecf51154046145cf21f9c6960664d10205b90736bb13c61cccf409c40f00
cells     19,716        metadata  183 Phase 1-4 columns + 13 Phase 5 columns
```

Phase 5 fields added, nothing overwritten: `program_P1_score` … `program_P8_score` ·
`dominant_malignant_program` · `dominant_malignant_program_label` · `dominant_program_score` ·
`dominant_program_margin` · `malignant_program_confidence` · `program_score_source`.

Non-malignant, Ambiguous and Excluded cells carry **projected** scores (NNLS with the spectra held
fixed) and are flagged `program_score_source == "projected"`; they never influenced a program.

**10 preservation guards and 8 reload validations passed**, asserted column by column against a
copy taken before anything was added: cells and order unchanged, all Phase 1–4 columns identical,
only the expected columns added, reductions and assays unchanged, `malignancy_refined`,
`scevan_clone` and `tumor_state_phase4` unchanged, frozen counts hold.

## 11. Limitations

**A.** n = 4 patients; `sample_id` = patient = dataset, so biological and technical effects cannot be
fully separated. **B.** Cells are not biological replicates. **C.** Malignant-cell abundance is
17.4× unequal and confounded with depth. **D.** A cNMF program is a continuous expression pattern,
not a discrete state, cell type, lineage or clone. **E.** The malignant-ECM comparison is evaluable
in **2 of 4 patients**, and the cross-patient log2FC correlation is 0.091. **F.** The signature is
descriptive; no classifier was trained and no pseudo-independent accuracy is reported.
**G.** Ambiguous fibroblasts were projected for description only and never reclassified.
**H.** Phase 4 CNA metrics reported alongside the groups are context, **not** independent validation
of the split that defined them. **I.** cNMF ran on unintegrated counts, so a patient-private
*biological* program cannot be distinguished from patient-level *technical* structure — no analysis
of these four patients can remove that. **J.** True generalization requires independent MPNST
patients.

## 12. SLURM accounting

| JobID | stage | State | Elapsed | CPUs | ReqMem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| 19899333 | build `p5_cnmf_env` | COMPLETED | 00:07:39 | 4 | 16 G | 0.52 GiB |
| 19899335 | **M36** feasibility | COMPLETED | 00:02:05 | 4 | 96 G | 29.28 GiB |
| 19899339 | **M37** prepare + factorize (3 runs) | COMPLETED | 00:31:01 | 8 | 48 G | 1.92 GiB |
| 19899353 | M37 consensus (12 K × 2 dt × 3 runs) | COMPLETED | 00:21:29 | 4 | 32 G | — |
| 19899354 | M37d K selection + **M38** + **M39** | COMPLETED | 00:02:37 | 4 | 64 G | 2.19 GiB |
| 19899355 | M40b LOO ×4 + high-confidence cNMF | COMPLETED | 00:08:15 | 8 | 48 G | 1.83 GiB |
| 19899356 | M38 re-run with technical labelling | COMPLETED | 00:02:13 | 4 | 64 G | — |
| 19899358 | M40c technical-gene-free cNMF (full grid) | COMPLETED | 00:25:52 | 8 | 48 G | 1.89 GiB |
| 19899359 | **M40** robustness + **M41** object build | **FAILED** (figures) | 00:10:55 | 4 | 128 G | 22.50 GiB |
| 19899742 | M40 re-run + figures 01–07, 10–11 | **FAILED** (figure 08) | 00:03:20 | 4 | 48 G | 2.63 GiB |
| 19899744 | figures 08/09/12 + **M41 freeze** | COMPLETED | 00:00:56 | 2 | 32 G | 0.19 GiB |
| 19899747 | figure re-render after visual QC | COMPLETED | 00:01:10 | 2 | 32 G | 0.42 GiB |
| 19899750 | figure re-render (robustness panel) | COMPLETED | 00:01:09 | 2 | 32 G | 0.47 GiB |

**Peak memory 29.28 GiB against the 450 G envelope — 6.5%.** No failure was addressed by increasing
memory or walltime.

**Two failures, both root-caused and fixed at source, both preserved above.**
*19899359* — `m41_figures.R` looked for `postint_umap_harmony_1/2` in the exported metadata, but UMAP
coordinates live in the object's **reduction**, not in `meta.data`. Fixed by reading the frozen
Phase 4 embedding table. The object build in the same job had already completed and passed all 8
reload validations, so it was not repeated.
*19899742* — inside a single `mutate()`, a new column named `sig` shadowed the `sig` signature data
frame referenced two lines later. Renamed to `is_sig`.

## 13. Environment

`R_env` **unchanged**: `R_env_PRE_PHASE5.yml` and `R_env_POST_PHASE5.yml` are byte-identical
(asserted in the job script). One isolated addition, `p5_cnmf_env` (python 3.11, cnmf 1.7.1,
numpy/scipy/scikit-learn/anndata/scanpy), which nothing else in the project depends on. Gene sets
downloaded as GMTs to `external/genesets/` with md5s in `MSIGDB_CHECKSUMS.md5`.

## 14. Deliverables

**12 final figures** in `results/phase5/figures/final/` (PDF + PNG) · **16 final tables** in
`results/phase5/tables/final/` · manifest `results/phase5/phase5_manifest.json` (30 sections) ·
`reports/phase5/{PHASE5_FEASIBILITY.md, ROBUSTNESS_REPORT.md, milestones/M36–M41}` ·
`reports/FIGURE_INDEX.tsv` 682 → 706 rows.

## 15. What Phase 5 answers

1. **Do recurrent continuous programs exist despite patient-private discrete states?** **No.**
   0 of 8 recurrent, and 0 recurrent at 11 of 12 ranks on the technical-gene-free universe.
2. **Does patient imbalance determine the solution?** **No** — every program is recovered in the
   balanced run. But patient identity does *confine* each program to one patient.
3. **Are the five Mesenchymal_ECM states one program?** **No** — they map to three, though MPNST_4's
   three ECM states do collapse onto one.
4. **Which programs are recurrent vs private?** 0 / 1 shared-limited (technical) / 7 private. P3
   (ECM) reaches three patients only at a relaxed 10% activity threshold.
5. **What distinguishes malignant ECM-like cells from genuine fibroblasts?** A 92-gene signature,
   consistent in the 2 evaluable patients, with the fibroblast side marking nerve-associated /
   endoneurial stroma.
6. **Can a transparent signature be derived?** **Yes** — `MALIGNANT_ECM_SIGNATURE.tsv`, no classifier.
7. **Do the findings survive the sensitivity analyses?** **Yes**, and the central negative result
   survives the most aggressive one (technical gene removal across the whole K grid).
