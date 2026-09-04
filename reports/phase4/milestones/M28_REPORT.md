# M28 — Phase 4 Reconstruction and SCEVAN Feasibility

**Status: COMPLETE** · 2026-09-03 · continuing automatically to M29

## Tasks (§11)

| # | Task | Outcome |
| --- | --- | --- |
| 1 | reconstruct Phase 2/3 | done — git log, PROJECT/README/PROGRESS/CHANGELOG, both handoffs, both manifests, `MPNST_CCC_MASTER_TABLE.tsv` |
| 2 | identify correct raw-count representation | RNA assay, four per-sample `counts` layers; integer-verified |
| 3 | verify SCEVAN input compatibility | 89.2% of features map to SCEVAN's chr1–22 annotation |
| 4 | inspect gene identifiers | **gene symbols** (0/31,764 are Ensembl IDs) |
| 5 | inspect cell IDs | 19,716; 0 duplicates; round-trip `setequal` TRUE |
| 6 | inspect patient/sample structure | 4 samples; `sample_id` = patient = dataset |
| 7 | inspect normal-reference candidates | 7,448 high-confidence immune; 11,308 disputed excluded |
| 8 | install/validate SCEVAN | **already installed** (1.0.3 + yaGST); validated against the installed API |
| 9 | small technical smoke test | SUCCESS, full pipeline including subclones |

## Input

`results/phase2/phase2_final_object.rds` — md5 `153d5f6acc70f9c05aa48cabc4f4ac2d`, **verified
equal** to the Phase 2 manifest. 19,716 cells, 149 metadata columns, reductions `pca`,
`umap_preintegration`, `postint_harmony`, `postint_umap_harmony`. Read-only.

## Output

```text
results/phase4/scevan/by_sample/MPNST_{1,2,3,4}/counts_raw.rds   sparse integer RNA counts
results/phase4/scevan/m28_feasibility_facts.json                 full structured audit
results/phase4/malignancy/phase4_cell_metadata.tsv.gz            19,716 x 155 (both UMAPs)
results/phase4/malignancy/phase2_harmony_embedding.rds           frozen copy, 19,716 x 30
results/phase4/scevan/smoke_test/                                74 files, technical only
reports/phase4/SCEVAN_FEASIBILITY.md
reports/phase4/MALIGNANCY_DECISION_RULES.md   (written BEFORE any SCEVAN result)
reports/phase4/environment/{R_env_PRE_PHASE4.yml, p4_umap_env_PHASE4.yml, PHASE4_INSTALL_LOG.md}
```

## Scripts

`scripts/phase4/scevan/m28_feasibility.R` · `scripts/shell/phase4/run_m28_feasibility.sh` ·
`scripts/shell/phase4/install_umap_learn.sh` · `scripts/phase4/utils/phase4_plot_utils.R`

## SLURM accounting

| JobID | Purpose | State | Elapsed | ReqMem | MaxRSS |
| --- | --- | --- | --- | --- | --- |
| 19896576 | M28 first attempt | **FAILED** | 00:02:15 | 96 G | 27.60 GiB |
| 19896623 | M28 after output_dir fix | **FAILED** | 00:03:14 | 96 G | 27.60 GiB |
| 19896654 | umap-learn env, attempt 1 | **FAILED** | 00:00:11 | 16 G | 0 |
| 19896672 | umap-learn env | COMPLETED | 00:05:39 | 16 G | 0.48 GiB |
| **19896711** | **M28 final** | **COMPLETED** | **00:05:07** | 96 G | **26.90 GiB** |

**No resource request was raised in response to any failure.** 96 G was set from the Phase 3
measured load and used 28% of the request.

### Failures, root-caused

1. **19896576** — `cannot open the connection` after `"found 363 tumor cells"`. SCEVAN 1.0.3
   hardcodes `path = "./output"` in `getScevanCNV`, `getScevanCNVfinal`, `plotAllClonalCN`,
   `plotAllSubclonalCN`, `plotConsensusCNA` and `analyzeSegm2`, and `plotCNclonal()` does not
   forward `output_dir`. **Fix**: each run sets its own working directory and uses SCEVAN's
   default `output_dir`, so the parameterised and hardcoded paths coincide. Package not patched.
2. **19896623** — `Python module umap was not found`. `subcloneAnalysisPipeline` → `plotTSNE`
   needs Python umap-learn via reticulate, and `plotTSNE` runs **before** the line writing
   subclone labels into `classDf`, so failing it loses the clone assignments. **Fix**:
   umap-learn 0.5.12 in an isolated `p4_umap_env`, reached via `RETICULATE_PYTHON`.
3. **19896654** — `QT_XCB_GL_INTEGRATION: unbound variable`, 11 s. Our own defect: `set -u`
   plus conda's `qt-main` activation hook, triggered by `conda deactivate`. **Fix**: drop `-u`
   from the Phase 4 shell scripts, per conda's own guidance, with the reason in the script.

## Warnings

* SCEVAN's annotation is **chr1–22 only**; X and Y events are not assessable.
* SCEVAN removes cell-cycle genes and **all `HLA-*` genes** before inference, so the CNV
  analysis is structurally blind to the HLA-E / HLA-F loci that Phase 3 highlighted.
* The Phase 2 default assay is `SCT` with four models; it is not used. RNA counts only.

## Scientific decisions

1. **Per-patient SCEVAN is primary** (§14) — `sample_id` = patient = dataset, so pooling first
   would confuse patient-specific CNV architecture with intra-tumour subclones.
2. **Primary run passes `norm_cell = NULL`** — confident-normal detection sees no Phase 2
   label, which is what gives the analysis inferential value (§16).
3. **`FIXED_NORMAL_CELLS = TRUE` is prohibited in this project.** The installed source shows it
   executes `cellType_pred[!cellType_pred %in% norm_cell_names] <- "malignant"`, forcing every
   non-reference cell to malignant — the banned non-immune-equals-tumour inference arriving
   through an argument.
4. **Fibroblast, Pericyte-VSMC, Candidate-Malignant-Unresolved and Uncertain are never fixed
   normal references**, and Endothelial is also excluded so its call stays independent of the
   Phase 3 findings that concern it.
5. **Decision rules committed before results** — `MALIGNANCY_DECISION_RULES.md` was written
   with every threshold fixed a priori, so no rule can be tuned to the answer.
6. **Per-sample counts written to disk** so the four M29 jobs never reload the 6 GB object.
7. **`as.matrix()` never applied to the full dataset** (§13); SCEVAN's own internal densify was
   measured instead and used to size the M29 request.

## Next milestone

**M29 — per-patient SCEVAN**, array 19896712_[1-4], one task per patient, 8 CPUs / 64 G / 12 h
each, `par_cores` matched to `--cpus-per-task`.
