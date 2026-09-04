# M37 — Continuous Malignant Program Discovery (cNMF)

**Date** 2026-09-03 · **Status** COMPLETE · **SLURM** 19899333 (env), 19899339 (factorize), 19899353 (consensus), 19899354 (K selection)

**Scientific question.** Do continuous transcriptional programs exist in the 6,434 refined malignant
cells, and at what rank should they be read?

**Method.** cNMF 1.7.1 in the **isolated** `p5_cnmf_env` (python 3.11, numpy/scipy/scikit-learn/
anndata/scanpy). `R_env` was not touched: no package was installed, upgraded or downgraded in it.

**Declared before any factor was inspected.**

```text
cells        malignancy_refined == "Malignant" only (6,434)
expression   RNA raw integer counts (cNMF normalises to TPM internally)
             NOT Harmony, NOT UMAP, NOT PCA scores, NOT SCT residuals, NOT the SCEVAN CNA matrix
genes        detected in >= 0.5% of malignant cells, ^MT- removed  ->  19,663
             ribosomal, cell-cycle, ECM, HLA, interferon, Schwann, neural-crest and
             angiogenesis genes ALL RETAINED in the primary run
K grid       4-15      replicates 100 per K      seed 42      numgenes 2000
consensus    local-density-threshold 0.10 (primary) and 2.00 (unfiltered, for comparison)
```

**Three runs, all declared in advance**: `primary` (6,434 cells), `balanced` (patient-balanced),
`nocc` (95 canonical cell-cycle genes removed). Two more were added later as declared sensitivity
analyses: `loo_<patient>` ×4 and `highconf` (M40b), and `nortech` (M40c).

**Patient-balanced design.** Cap = **651** cells, MPNST_2 being the smallest patient with ≥ 500
malignant cells; patients below the cap contribute all their cells. Result: 2,165 cells with the
dominant-patient fraction cut from **0.573 to 0.301**. Seed 42.

**K-selection rule, declared in `m37_k_selection.py` before it was run.** The largest K on the grid
satisfying: (i) max pairwise cosine between consensus spectra ≤ 0.75 (no duplicated program);
(ii) every program dominant in ≥ 1% of malignant cells (no dead program); (iii) stability ≥ the
median stability of the K values satisfying (i) and (ii). The rule references only measured
properties of the factorization — never a gene, a pathway or a label — so K cannot be chosen to
produce a pleasing biological answer.

| K | silhouette | prediction error | max program cosine | min dominant share | dead programs | no duplicate | no dead |
| ---: | ---: | ---: | ---: | ---: | ---: | :--: | :--: |
| 4 | 0.804 | 9.93e6 | 0.365 | 0.0404 | 0 | yes | yes |
| 5 | 0.834 | 9.65e6 | 0.442 | 0.0228 | 0 | yes | yes |
| 6 | **0.862** | 9.42e6 | 0.517 | 0.0233 | 0 | yes | yes |
| 7 | 0.834 | 9.23e6 | 0.498 | 0.0131 | 0 | yes | yes |
| **8** | **0.843** | 9.09e6 | 0.652 | 0.0131 | 0 | **yes** | **yes** |
| 9 | 0.823 | 8.98e6 | **0.875** | 0.0131 | 0 | **no** | yes |
| 10–15 | 0.78–0.80 | ↓ | 0.62–0.77 | ≤ 0.0005 | 1–3 | mixed | **no** |

Eligible K = {4, 5, 6, 7, 8}; median eligible stability 0.834; **selected K = 8**.

**Resources.**

| JobID | stage | State | Elapsed | CPUs | ReqMem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| 19899333 | build `p5_cnmf_env` | COMPLETED | 00:07:39 | 4 | 16 G | 0.52 GiB |
| 19899339 | prepare + factorize ×3 runs | COMPLETED | 00:31:01 | 8 | 48 G | **1.92 GiB** |
| 19899353 | consensus, 12 K × 2 dt × 3 runs | COMPLETED | 00:21:29 | 4 | 32 G | — |
| 19899354 | K selection + M38 + M39 | COMPLETED | 00:02:37 | 4 | 64 G | 2.19 GiB |

**Warnings.** Peak memory across the whole of M37 is **1.92 GiB** — the factorized matrix is only
6,434 cells × 2,000 over-dispersed genes. The cost is wall-clock (1,200 NMF fits per run), not RAM,
which is why 8 parallel workers were used rather than a larger memory request.

**Failures / fixes.** None.

**Next.** M38 — program annotation, pathways and the comparison against the Phase 4 discrete states.
