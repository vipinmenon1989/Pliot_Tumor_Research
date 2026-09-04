# M35A — SCEVAN Figure Audit

*Phase 4 · milestone M35A · SCEVAN figure consolidation and evidence visualization · 2026-09-03*

> **Scope.** M35A surfaces, consolidates and improves the visualization of **existing** Phase 4
> SCEVAN evidence. **No analysis was re-run.** SCEVAN was not re-executed, no malignancy call,
> decision rule, amendment, clone assignment, tumour state or CCC result was changed, and
> `results/phase4/phase4_final_object.rds` was neither loaded nor modified.
>
> **SCEVAN infers copy number from gene expression. It is not DNA sequencing and provides no
> DNA-level proof.** Inference is reliable for broad chromosomal, arm-level and large-segment
> events only. No single-gene CNV claim is made in any M35A figure.

---

## 1. What already existed

`results/phase4/tables/SCEVAN_NATIVE_FIGURE_INDEX.tsv` indexes **267 native SCEVAN files**
across 4 patients × 2 runs (primary, sensitivity). By type:

| native file pattern | indexed figure_type | n | biological question it answers |
| --- | --- | ---: | --- |
| `<S>_<run>heatmap.png` | CNA heatmap, all cells, with malignant/normal track | 8 | Which assessed cells does SCEVAN call malignant, and on what inferred CNA pattern? |
| `<S>_<run>onlytumorheatmap.png` | CNA heatmap, tumour cells only | 8 | What is the CNA architecture of the malignant compartment alone? |
| `<S>_<run>heatmap_subclones.png` | CNA heatmap with subclone annotation track | 8 | How does the malignant compartment partition into CNA-defined subclones? |
| `<S>_<run>consensus.png` | consensus clonal CN profile | 8 | What is the segment-level consensus CN profile of the tumour? |
| `<S>_<run>OncoHeat.png`, `OncoHeat2.png` | cytoband-level onco-heatmap by subclone | 16 | Which cytobands are gained/lost in which subclone? |
| `<S>_<run>umap_CNA.png` | UMAP of CNA space, coloured by subclone | 8 | Do the subclones separate in CNA space? |
| `<S>_<run>umap_scRNA.png` | UMAP of expression space, tumour vs normal | 8 | Do malignant and normal cells separate in expression space? |
| `<S>_<run>CloneTree.png` | (indexed as generic SCEVAN output) | 8 | **none — every file is blank** (see §2) |
| `<S>_<run>-DEchr*_subclones.png` | generic SCEVAN output | 151 | subclone-differential expression per segment |
| `<S>_<run>pathwayAnalysis_subclones*.png` | generic SCEVAN output | 44 | per-subclone pathway enrichment |

Machine-readable clonal CN profiles also exist and are **not** images:
`<S>_<run>_Clonal_CN.seg` and `<S>_<run>_subclone<N>_CN.seg` (one per subclone), plus the stored
matrices `<S>_<run>_CNAmtx.RData` (genes × cells relative CNA), `_CNAmtxSubclones.RData`,
`_count_mtx_annot.RData` (gene → chromosome/position) and `_PlotOncoHeat.RData`.

Sixteen Phase 4 final figures (`01`–`16`) plus the M29–M34 working figures already existed.
Of those, only `04_scevan_cna_heatmap`, `05_scevan_clones_by_patient`,
`29_03_clone_composition_by_annotation` and `08/16_cnv_profiles_by_annotation` touched the SCEVAN
CNA evidence, and none of them surfaced a native SCEVAN CNA heatmap.

## 2. Defects found in the existing native output

**All 8 `CloneTree.png` files are blank** (pixel standard deviation exactly 0). This is the
documented ggtree / ggplot2 4.x failure recorded in `PHASE4_HANDOFF.md` §4 defect 3 and
limitation K: `plotCloneTree` calls the removed `ggplot2:::is.waive()`, SCEVAN catches the error
internally, and every clone assignment and CN profile is intact — only the phylogeny plot is
lost. **ggplot2 was not downgraded and clone phylogeny plotting was not forced**, per the M35A
prohibitions. No M35A figure depends on it.

## 3. Audit table

`suitable` = usable as-is in a final presentation. `replacement needed` = a custom or companion
figure was created in M35A.

| existing figure | source file (primary run shown; sensitivity exists for all) | patient | run | biological question | suitable as-is | replacement / companion needed |
| --- | --- | --- | --- | --- | :--: | --- |
| CNA heatmap, subclone-annotated | `results/phase4/scevan/by_sample/MPNST_1/primary/output/MPNST_1_primaryheatmap_subclones.png` | MPNST_1 | primary | CNA-defined subclone architecture | yes — but unlabelled and buried in an index | **surfaced** as `17_scevan_native_cna_MPNST1.pdf` p1 + companion annotation panel |
| CNA heatmap, subclone-annotated | `.../MPNST_2/primary/output/MPNST_2_primaryheatmap_subclones.png` | MPNST_2 | primary | as above | yes | **surfaced** as `18_scevan_native_cna_MPNST2.pdf` p1 + companion panel |
| CNA heatmap, subclone-annotated | `.../MPNST_4/primary/output/MPNST_4_primaryheatmap_subclones.png` | MPNST_4 | primary | as above | yes | **surfaced** as `19_scevan_native_cna_MPNST4.pdf` p1 + companion panel |
| CNA heatmap, subclone-annotated | `.../MPNST_3/primary/output/MPNST_3_primaryheatmap_subclones.png` | MPNST_3 | primary | as above | **no** — partition is an immune artefact | not surfaced as evidence; MPNST_3 is handled as a failure in `25_MPNST3_scevan_failure_qc.pdf` |
| CNA heatmap, all cells + malignant/normal track | `.../<S>_primaryheatmap.png` | 1, 2, 4 | primary | which cells SCEVAN called malignant | yes | **surfaced** on page 2 of `17`/`18`/`19` |
| CNA heatmap, tumour cells only | `.../<S>_primaryonlytumorheatmap.png` | 1, 2, 4 | primary | malignant-compartment CNA architecture | yes | **surfaced** on page 2 of `17`/`18`/`19` |
| consensus clonal CN profile | `.../<S>_primaryconsensus.png` | 1, 2, 4 | primary | segment-level consensus CN | yes | **surfaced** on page 2 of `17`/`18`/`19` |
| clonal / subclonal CN profiles | `.../<S>_primary_Clonal_CN.seg`, `_subclone<N>_CN.seg` | all | primary + sens | per-clone segment CN | n/a (text) | aggregated into `SCEVAN_CNV_SUMMARY.tsv` → `27_broad_cna_recurrence.pdf` |
| cytoband onco-heatmap | `.../<S>_primaryOncoHeat.png`, `OncoHeat2.png` | all | primary + sens | cytoband gain/loss per subclone | partly — cytoband labels are unreadable at page scale and invite gene-level over-reading | **not surfaced.** Broad-event evidence is presented instead at segment level in `27_broad_cna_recurrence.pdf`, which is the resolution the method supports |
| UMAP of CNA space by subclone | `.../<S>_primaryumap_CNA.png` | all | primary + sens | subclone separation in CNA space | yes, but redundant | not surfaced — `05_scevan_clones_by_patient` and `31_03_scevan_clone_umap` already cover clone structure on the project UMAP |
| UMAP of expression space, tumour vs normal | `.../<S>_primaryumap_scRNA.png` | all | primary + sens | malignant/normal separation in expression space | yes, but redundant | not surfaced — `01/02_scevan_malignancy_umap` already shows this on the Phase 2 Harmony UMAP for all cells at once |
| clone phylogeny | `.../<S>_<run>CloneTree.png` | all | both | — | **no — blank file** | **none.** Clone phylogeny plotting is not forced; the clone/state relationship is reported as a cross-tabulation (`CLONE_VS_TUMOR_STATE.tsv`, figure `10`) |
| subclone DE per segment | `-DEchr*_subclones.png` (151) | all | both | subclone-differential expression | no — 151 single-segment panels | not surfaced; out of M35A scope |
| per-subclone pathway analysis | `pathwayAnalysis_subclones*.png` (44) | all | both | subclone pathway enrichment | no | not surfaced; out of M35A scope |
| Phase 4 `04_scevan_cna_heatmap` | `results/phase4/figures/final/04_scevan_cna_heatmap.pdf` | all | primary | project-drawn CNA overview | yes | retained unchanged; complemented, not replaced, by `17`–`20` |
| Phase 4 `29_03_clone_composition_by_annotation` | `results/phase4/figures/final/29_03_clone_composition_by_annotation.pdf` | all | primary | clone composition | partial — all four patients including MPNST_3, no clone size, no refined-malignancy annotation | **replaced for presentation** by `21_scevan_clone_composition_phase2.pdf`; the original is retained as history |
| Phase 4 `08/16_cnv_profiles_by_annotation` | `results/phase4/figures/final/16_cnv_profiles_by_annotation.pdf` | all | primary | CNA profile by Phase 2 annotation | partial — split by annotation only, so it cannot separate malignant from non-malignant fibroblasts | **complemented** by `23_fibroblast_malignant_vs_nonmalignant_cna_profile.pdf` |

**Nothing that already existed was recreated.** The native SCEVAN CNA heatmaps are embedded
verbatim as rasters — the CNA matrix, its cell ordering and its subclone track are SCEVAN's own
output and were not redrawn. Modifying the native plot object itself would have been fragile
(the matrix is not returned in a plot-ready form and the ordering is internal), so the Phase 2 /
clone / malignancy annotation the native plot cannot carry is supplied as an **aligned companion
panel built from the same frozen clone assignments**, exactly as the M35A brief specifies.

## 4. New M35A figures

| figure | what it adds that no existing figure provided |
| --- | --- |
| `17`–`19_scevan_native_cna_MPNST{1,2,4}.pdf` | the native subclone-annotated CNA heatmap, standardised and labelled, beside a companion panel giving each clone's Phase 2 composition, size and refined-malignant fraction |
| `20_scevan_cna_reliable_patients.pdf` | the three reliable patients' CNA/subclone architecture on one page, making the patient-specific structure immediately visible |
| `21_scevan_clone_composition_phase2.pdf` | **the central figure** — clone × Phase 2 identity, faceted by patient, with clone size, refined-malignant fraction and fibroblast-dominance marked |
| `22_fibroblast_cna_burden_by_patient.pdf` | patient-stratified CNA metric distributions with per-patient Cliff's delta instead of a pooled cell-level p-value |
| `23_fibroblast_malignant_vs_nonmalignant_cna_profile.pdf` | genome-wide mean CNA profiles of malignant vs non-malignant fibroblast-labelled cells, from the stored native CNA matrices |
| `24_phase2_to_phase4_malignancy_transition.pdf` | the identity transition for the disputed compartments only, with immune populations deliberately excluded |
| `25_MPNST3_scevan_failure_qc.pdf` | MPNST_3 shown failing — run agreement, immune-lineage clones, immune malignant fractions, depth as context only |
| `26_fibroblast_malignancy_threshold_robustness.pdf` | separates the threshold-stable fibroblast count from the threshold-sensitive cohort fraction |
| `27_broad_cna_recurrence.pdf` | broad recurrent events across the reliable patients at segment resolution, with no gene-level claim |
| `28_phase4_scevan_evidence_summary.pdf` | one publication-style page carrying the whole argument |

## 5. Discrepancy found against the frozen handoff

Every frozen count reproduced exactly (19,716 cells; Fibroblast 5,064 → 4,036/908/120;
Candidate-Malignant-Unresolved 836/395/0; MPNST-Tumor 1,405/2,015; refined
6,434/9,078/3,766/438; clones 7/4/3/8; run agreement 0.9951/0.9663/0.0864/0.9974). **One wording
discrepancy was found and is recorded rather than smoothed over:**

`PHASE4_HANDOFF.md` §12 states that "**Every** MPNST_1 clone mixes
`Candidate-Malignant-Unresolved`, `Fibroblast` and `MPNST-Tumor` together". The frozen clone
assignments show this holds for **6 of 7** MPNST_1 clones. `MPNST_1_clone6` (176 cells) contains
`Candidate-Malignant-Unresolved` = 108 and `MPNST-Tumor` = 68 and **no `Fibroblast` cells at
all**. `MPNST_1_clone7` contains exactly 1. The accurate statement is:

* **all 7** MPNST_1 clones mix `Candidate-Malignant-Unresolved` with `MPNST-Tumor`;
* **6 of 7** additionally contain `Fibroblast` cells;
* the two clones that do not are the two smallest and second-smallest of the mixed set
  (`clone6` n = 176, `clone7` n = 54), together 230 of 1,952 MPNST_1 malignant cells.

This does **not** change any malignancy call, clone assignment, count or conclusion — the
underlying claim (marker annotation split into three identities what copy number sees as one
malignant compartment with seven subclones) stands on all seven clones. Figure 21 plots the
composition as it actually is, and the M35A section appended to the handoff records the
correction without rewriting the earlier text.

## 6. What the figures show, honestly

* **Clone composition is the strongest evidence.** SCEVAN builds subclones from inferred CNA
  profiles without reference to any Phase 2 label, so the Phase 2 composition of a CNA-defined
  clone is not circular. MPNST_2: 4 of 4 clones fibroblast-dominated. MPNST_4: 7 of 8 (clone 8 is
  endothelial, 146/152 — flagged by rule R11, not promoted). MPNST_1: all 7 clones mix the
  disputed identities as described in §5.
* **Malignant-called fibroblasts differ from non-malignant ones in CNA burden, in the two
  patients where the comparison exists.** MPNST_1 median-level separation (mean `cnv_burden`
  0.365 vs 0.194, n = 512 vs 303) and MPNST_2 (0.315 vs 0.194, n = 637 vs 503). **MPNST_4 retains
  only 1 non-malignant fibroblast**, so its contrast is reported as NOT EVALUABLE rather than
  pooled away, and MPNST_3 contributes no malignant fibroblasts because its promotions are
  disabled. The result is therefore consistent where testable, on 2 of 4 patients — stated as
  such, not inflated.
* **The fibroblast count is threshold-stable; the cohort fraction is not.** 4,036 at every tested
  `pop_frac_low` from 0.15 to 0.40, against `MPNST-Tumor` retention moving 3,266 → 1,405 and the
  refined fraction moving 42.07% → 32.63%.
* **MPNST_3 failed, and the figures say so.** Primary↔sensitivity agreement 0.0864 against
  0.9663–0.9974 elsewhere; its three "clones" are T/NK, plasma/pDC/B and plasma-dominated; six of
  nine canonical immune populations are called ~100% malignant. Depth is shown as context and is
  explicitly **not** claimed as the cause.

## 7. Execution record

All figure work ran through `sbatch` on `ihc-grid-1-1-1` (account `ihc`, partition `ihc`).

| JobID | stage | State | Elapsed | AllocCPUS | ReqMem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| 19899301 | figures, first pass | COMPLETED | 00:04:29 | 2 | 24 G | 2.11 GiB |
| 19899307 | figures, label/wrapping fixes | COMPLETED | 00:04:25 | 2 | 24 G | 1.95 GiB |
| 19899308 | figures, summary-layout fix | COMPLETED | 00:04:27 | 2 | 24 G | 1.96 GiB |
| 19899311 | figures, genome-axis and panel-order fix | COMPLETED | 00:04:25 | 2 | 24 G | 2.01 GiB |
| **19899312** | **figures, final** | COMPLETED | 00:04:24 | 2 | 24 G | **2.97 GiB** |
| 19899313 | finalize, first attempt (output discarded, §8) | COMPLETED | 00:00:10 | 1 | 8 G | — |
| **19899315** | **finalize** — checksums, figure index, manifest | COMPLETED | 00:00:09 | 1 | 8 G | — |

Every run exited 0. The four figure reruns are quality iterations found by rendering each PDF and
inspecting it, not failures — the defects fixed were a clipped subtitle, an overlapping percentage
label, a mis-dodged "NOT EVALUABLE" annotation, a per-facet chromosome axis that mislabelled two of
three patients, a patchwork `design` whose area letters did not match the order the plots were
added, and a subclone count that was reporting cells instead of clones. **No plotting problem was
addressed by increasing memory or walltime.** Peak 2.97 GiB against a 24 G request; the 6 GB final
object was never loaded.

## 8. A finalize run that exited 0 and was still wrong

Job 19899313 completed successfully and wrote a corrupted manifest. `jsonlite::toJSON` defaults to
`digits = 4`, so re-serialising `phase4_manifest.json` silently rounded frozen values that were
already in it — run agreement `0.99514117 → 0.9951`, refined fraction `0.32633394 → 0.3263`, every
`elapsed_min`, every `pop_frac`, every CNA reference threshold. Nothing errored.

It was caught by diffing the written manifest against a backup taken before the write. The manifest
was restored from that backup and regenerated under job 19899315 with `digits = NA`.
`m35a_finalize.R` now snapshots the manifest before writing and, afterwards, **re-reads it and
asserts every pre-existing section is identical to the snapshot** — 33 sections verified unchanged,
with only `figures`, `tables` (appended to) and the new `m35a_figure_consolidation` block differing.
`phase4_object` and its md5 `e85ba8486e456917e2483f2773bdbaf3` are asserted before and after.
