# Milestone M18 — Phase 3 Reconstruction and Method Audit

**Phase 3 · MPNST tumour–immune cell–cell communication** · *Generated 2026-09-03*
*Status: **COMPLETE** — proceeded automatically to M19.*

## 1. Deliverables

| Document | Content |
| --- | --- |
| `reports/phase3/PERTURBSEQ_CODE_REUSE_AUDIT.md` | 33 components assessed: 9 ADAPT, 11 REFERENCE ONLY, 13 NOT APPLICABLE |
| `reports/phase3/LOCHNESS_IMPLEMENTATION_AUDIT.md` | Publication vs official MMCA R vs our Python implementation |
| `reports/phase3/PHASE3_METHOD_PLAN.md` | Full Phase 3 method architecture |
| `reports/phase3/environment/` | Pre/post environment exports and the install log |

## 2. Phase 2 state reconstructed from repository evidence

| Property | Value |
| --- | --- |
| Final object | `results/phase2/phase2_final_object.rds`, md5 `153d5f6acc70f9c05aa48cabc4f4ac2d` |
| Cells / features | 19,716 · RNA 31,764 / SCT 29,113 (4 SCT models) |
| CCC field | `annotation_ccc` (17 identities) · detailed: `postint_celltype_level*_refined` |
| Sample field | `sample_id` — simultaneously dataset **and** patient |
| Condition field | **does not exist** |
| Tumour | `MPNST-Tumor` 3,420 cells (17.35%) |
| Not CCC-ready | `Candidate-Malignant-Unresolved` 1,231 · `Uncertain` 720 · `Low-quality-excluded` 438 |

## 3. Perturb-seq audit — headline decisions

**READ-ONLY. Nothing modified; Phase 3 has no runtime dependency on it.**

The valuable component is the **LochNESS core formula** (`local/global − 1`), which matches
the publication exactly — **ADAPT**, with two corrections: add same-sample exclusion (absent
there) and switch *k* to `round(0.5·√N)`.

Two components were rejected on **scientific** rather than technical grounds: energy
distance / MMD (they require a designed control, which MPNST lacks) and ORA / gene-programme
enrichment (Phase 3 §69 prohibits large unrelated pathway screens).

## 4. LochNESS audit — headline findings

Official MMCA (`shendurelab/MMCA` @ `af629c498f421d6f0bfb325e9fd6b2dab0e22837`) is a
**repository of scripts, not an installable R package** — it has no `DESCRIPTION` or `R/`
directory, so `remotes::install_github()` would fail and was **not attempted**. Three scripts
were vendored read-only to `external/MMCA_ref/`.

Six divergences between the official R implementation and our Python port were identified.
The two that matter scientifically:

1. **Same-sample exclusion.** MMCA excludes same-embryo neighbours *by construction* —
   query and reference cell sets are disjoint. **Our Python implementation has no such
   exclusion.** In MPNST, where sample = patient = batch and 78% of nearest neighbours are
   still same-sample after Harmony, omitting it would turn LochNESS into a batch detector.
2. **k.** MMCA uses `round(0.5·√N)`; the Python version uses a fixed configured ~300.

**Decision: implement LochNESS natively in R following MMCA**, with the Python version kept
as an independent cross-check in M24.

## 5. Environment and installation

`conda env export > reports/phase3/environment/R_env_PRE_PHASE3.yml` was taken before any
change. R 4.4.3, Seurat 5.4.0, SeuratObject 5.3.0, harmony 1.2.4, Matrix 1.7.4 —
**all verified unchanged after every install stage.**

Compute nodes have full outbound network access (verified: CRAN, GitHub, Bioconductor, PyPI,
Zenodo all HTTP 200 from `ihc-grid-1-1-1`), so installation runs inside SLURM.

Installed for Phase 3: **liana 0.1.14**, OmnipathR 3.14.0, decoupleR 2.12.0, basilisk 1.18.0,
NMF 0.28, ggraph 2.2.2, ggpubr 1.0.0, svglite 2.2.2, FNN 1.1.4.1, plus a CRAN/Bioconductor
support set. **CellPhoneDB 5.0.1** was installed in a **deliberately isolated** conda
environment `cpdb_env` (Python 3.10) because its numpy/pandas pins would destabilise
`R_env`. Full detail in `reports/phase3/environment/PHASE3_INSTALL_LOG.md`.

## 6. SLURM accounting

| Job | JobID | State | Elapsed | Notes |
| --- | --- | --- | --- | --- |
| Compute-node network probe | 19894719 | COMPLETED | 00:00:02 | all endpoints reachable |
| Dependency install (stage 1–7) | 19894761 | COMPLETED | 00:43:55 | liana OK; CellChat + nichenetr failed on a compile chain |
| Dependency retry 1 | 19895124 | COMPLETED | 00:04:45 | conda-forge binaries fixed 14 packages; `systemfonts` still 1.2.3 < 1.3.0 |
| Dependency retry 2 | 19895459 | see M20 | — | pins `systemfonts ≥ 1.3.2` + nichenetr system-library deps |

**Failures were diagnosed, not brute-forced.** Root cause of the CellChat/nichenetr failures:
source compilation of `systemfonts` and `Deriv` failed, cascading to svglite, doBy, pbkrtest,
car, rstatix and ggpubr; nichenetr additionally needed `shadowtext` and the
`units`/`sf`/`gdtools` system-library chain. The fix was prebuilt conda-forge R 4.4 binaries
with `r-base` pinned, **not** an increase in memory or time and **not** a downgrade of any
Phase 2 package.

## 7. Next

Proceeded automatically to **M19 — sample-aware CCC input preparation**.
