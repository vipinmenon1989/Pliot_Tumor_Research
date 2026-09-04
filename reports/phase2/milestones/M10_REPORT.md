# Milestone M10 — Phase 2 Reconstruction and Integration Planning

**Phase 2 · MPNST single-cell integration**
*Generated: 2026-09-02*
*Status: **COMPLETE — STOPPED, awaiting researcher approval***

> **No Harmony was executed in M10.** No Phase 1 object, reduction, graph or metadata
> column was modified. The only compute performed was a read-only structural inspection
> of the Phase 1 handoff object inside a SLURM allocation.

---

## 1. Executive summary

Phase 1 is frozen and internally consistent at the level of its machine-readable
records. Its combined pre-integration object verifies bit-for-bit against **both**
checksums recorded in Phase 1 provenance, and is the correct and only sensible Phase 2
input.

Three findings require researcher attention before M11:

1. **The Phase 1 `SCT` assay contains four SCTransform models, not one.** Phase 1
   documentation states a single global SCTransform was run. Structurally it was run
   layer-wise, once per dataset. This does **not** invalidate integration, but it
   changes the M14 marker plan (`PrepSCTFindMarkers()` becomes mandatory) and is a real
   discrepancy between the Phase 1 report text and the Phase 1 object. (§4.2)
2. **`percent.mt` is identically zero in all 7,615 MPNST_1 cells**, yet Phase 1 M8
   passed `vars.to.regress = "percent.mt"` to SCTransform unconditionally. For the
   MPNST_1 model this regresses a zero-variance covariate. A numerical-integrity check
   is added to the front of M11. (§4.3)
3. **Batch and biology are perfectly confounded** — dataset = sample = patient =
   presumed technical batch. Phase 2 §8 requires an explicit researcher decision before
   Harmony runs. (§5, and `reports/phase2/HARMONY_VARIABLE_DECISION.md`)

None of the three is judged to invalidate Harmony integration. Finding 1 and finding 2
are reported rather than silently repaired, per Phase 2 §5 and §29.

---

## 2. Repository state at milestone start

```
Branch          : dev  (up to date with origin/dev)
HEAD            : 751a535  "Clean rn and Phase 1 finished"
Working tree    : clean
Project root    : /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
                  (== /autofs/projects-t3/... ; same filesystem, verified with readlink -f)
```

Documents read: `PROJECT.md`, `PROGRESS.md`, `CHANGELOG.md`, `README.md`,
`reports/PHASE1_HANDOFF.md`, `reports/PRE_INTEGRATION_ASSESSMENT.md`,
`reports/INTEGRATION_PREPARATION.md`, `reports/milestones/M0`–`M9_REPORT.md`,
`reports/DATA_AUDIT.md`, and all Phase 1 machine-readable tables and provenance JSONs.

---

## 3. What Phase 1 actually produced

### 3.1 Milestones

M0–M9 all recorded **Completed**; `PROGRESS.md` declares **PHASE 1 = FROZEN**,
**PHASE 2 = NOT STARTED**. `reports/milestones/M9_REPORT.md` records "Unresolved
Issues: **None**".

### 3.2 Per-dataset artefacts (`results/datasets/{ds}/`)

For each of `MPNST_1`, `MPNST_2`, `MPNST_3`, `MPNST_4`: `_raw.rds`, `_filtered.rds`,
`_filtered_specific.rds`, `_normalized.rds`, `_pca.rds`, `_clustered.rds`, plus four
provenance JSONs each. Marker tables live under
`reports/datasets/{ds}/markers/resolution_{r}/`.

Phase 1 outcome per dataset:

| Dataset | Raw cells | Retained cells | Recommended PCs | Recommended res | Clusters | Alt res |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| MPNST_1 | 8,338 | 7,615 | 1–8 | 0.6 | 18 | 0.3 |
| MPNST_2 | 2,830 | 2,284 | 1–6 | 0.3 | 9 | 0.5 |
| MPNST_3 | 3,682 | 2,940 | 1–9 | 0.6 | 13 | 0.3 |
| MPNST_4 | 7,811 | 6,877 | 1–5 | 0.7 | 14 | 0.5 |
| **Total** | **22,661** | **19,716** | | | | |

Doublet detection: `scDblFinder`, run per dataset. QC thresholds were dataset-specific
(Strategy B) — `config/config.yaml: dataset_specific_qc`.

### 3.3 The combined pre-integration handoff (M8)

| Property | Value |
| --- | --- |
| Path | `results/combined/pre_integration/combined_preintegration.rds` |
| Size | 6,041,513,977 bytes |
| Cells | 19,716 |
| Default assay | `SCT` (29,113 features) |
| Secondary assay | `RNA` (31,764 features) |
| Reductions | `pca` (50 dims, assay SCT, loadings + stdev present), `umap_preintegration` (2 dims) |
| Graphs | `SCT_nn`, `SCT_snn` (built on PCs 1–30, k = 15) |
| Metadata columns | 93 |
| Random seed | 42 |

### 3.4 Phase 1 → Phase 2 entry contract

`reports/PHASE1_HANDOFF.md` §4 mandates: `preint_*` columns and Phase 1 coordinates are
immutable; all Phase 2 coordinates and metadata must use a `postint_*` namespace;
Phase 2 must consume `results/phase1_manifest.json` as its machine-readable input.
Phase 2 adopts all three (see `reports/phase2/PHASE2_PLAN.md` §3).

---

## 4. Phase 1 handoff validation (SLURM job 19886411)

Script: `scripts/R/phase2/inspect_phase1_handoff.R` (read-only)
Launcher: `scripts/shell/phase2/run_m10_inspect.sh`
Outputs: `results/phase2/handoff/`

### 4.1 Checksum verification — **PASS**

| Algorithm | Computed (M10) | Recorded in Phase 1 | Match |
| --- | --- | --- | --- |
| md5 | `88a442688f912d882f6c6da01820e329` | `results/combined/pre_integration/prov_analysis.json` | ✅ |
| sha256 | `c3fdce8b61602725977f989d8bbf10030ffe48ef6140372a13d2c59903a8dc66` | `results/phase1_manifest.json` | ✅ |

The Phase 1 handoff object is byte-identical to the object Phase 1 signed off. The two
different-looking checksums in the Phase 1 records are simply md5 and sha256 of the same
file, not a conflict.

### 4.2 FINDING 1 — the SCT assay holds **four** SCTransform models

Observed (`results/phase2/handoff/handoff_assays.tsv`):

| Assay | Class | Default | Features | Layers | SCT models |
| --- | --- | :---: | ---: | --- | ---: |
| `RNA` | `Assay5` | no | 31,764 | `counts.MPNST_1.1; data.MPNST_1.1; scale.data.1; …; counts.MPNST_4.4; data.MPNST_4.4; scale.data.4` (split, 4 sets) | — |
| `SCT` | `SCTAssay` | **yes** | 29,113 | `counts; data; scale.data` (joined) | **4** (`model1`, `model1.1`, `model1.2`, `model1.3`) |

`reports/PRE_INTEGRATION_ASSESSMENT.md` §3.2 and `PHASE1_HANDOFF.md` §3.1 both describe
this as "a single unified `SCTransform` globally on the merged raw counts". What actually
happened is the standard Seurat v5 behaviour: because `merge()` left the `RNA` assay
split into four per-dataset layers, `SCTransform()` ran **once per layer** and produced
four models, with the residual and corrected-count matrices then merged.

**Assessment — does this invalidate integration? No.**
Per-sample SCT residuals feeding a shared PCA, followed by embedding-level batch
correction, is a standard and defensible design (it is essentially Seurat's own
SCT-based integration recipe). The shared PCA, neighbour graph and UMAP that Phase 1
froze are all valid as a *pre-Harmony baseline*, because they were computed once on the
merged residual matrix.

**Consequences that must be carried forward:**

- **M14 (markers): `PrepSCTFindMarkers()` is now mandatory** before `FindAllMarkers()`
  on the `SCT` assay, because differential expression across cells normalised under
  different SCT models is otherwise invalid. Phase 2 §14 asks that Seurat v5 layer
  preparation be inspected explicitly rather than guessed — this is that inspection, and
  the answer is: yes, preparation is required.
- If the `RNA` assay is used for markers instead, `JoinLayers()` is required first — the
  `RNA` counts and data layers are split four ways.
- Residual per-sample normalisation differences persist in the SCT residuals. Harmony
  operates downstream of them, at the embedding level, and does not remove them from
  expression values. This is consistent with Phase 2 §7 and is stated as a limitation.

**No repair is proposed and none was performed.** Re-running a genuinely single-model
SCTransform would mean redoing Phase 1 M8, which Phase 2 §5 forbids without a finding
that integration is invalidated. The recommendation is to correct the Phase 1 *wording*
(a documentation fix), not the Phase 1 *object*.

### 4.3 FINDING 2 — zero-variance `percent.mt` regression in MPNST_1

From `results/phase2/handoff/qc_summary_by_sample.tsv`:

| Sample | nCount_RNA (median) | nFeature_RNA (median) | percent.mt mean ± sd | percent.ribo mean ± sd |
| --- | ---: | ---: | --- | --- |
| MPNST_1 | 7,503 | 3,139 | **0 ± 0** | 5.20 ± 3.42 |
| MPNST_2 | 7,944 | 2,490 | 7.17 ± 2.65 | 15.81 ± 5.84 |
| MPNST_3 | 4,265 | 1,594 | 3.29 ± 1.81 | 17.58 ± 7.92 |
| MPNST_4 | 5,008 | 2,114 | 8.78 ± 3.98 | 15.82 ± 5.41 |

MPNST_1 has no detectable mitochondrial transcripts at all (a known Phase 1 M1/M2
finding — 13 `^MT-` genes exist in the object but MPNST_1 registers 0.00%). Phase 1 M4
explicitly handled this case; Phase 1 M8's `analyze_pre_integration.R` did **not**, and
passed `vars.to.regress = "percent.mt"` unconditionally. For the MPNST_1 model this
regresses a constant.

The Phase 1 PCA, UMAP and mixing diagnostics all completed and look numerically sane, so
this did not produce a catastrophic failure. To close it rather than assume it, **M11
will begin with an explicit numerical-integrity check**: count non-finite values in
`SCT@scale.data` and in `Embeddings(obj, "pca")`, per sample, and abort before Harmony if
any are found. Result recorded in the M11 report either way.

Note also the strong per-sample QC asymmetry above (MPNST_1: 3,139 median genes and
5.2% ribo; MPNST_3: 1,594 median genes and 17.6% ribo). Library depth and complexity
differ substantially between samples and are part of what `sample_id` absorbs.

### 4.4 Confirmed structural facts for M11 design

- `pca` has 50 dimensions available; the frozen baseline used 1–30.
- `DefaultAssay(obj[["pca"]]) == "SCT"`, so `RunHarmony(reduction.use = "pca")` needs no
  `assay.use` argument.
- `sample_id` and `orig.ident` are cell-for-cell identical (`true` in
  `handoff_structure.json`).
- Graph names `SCT_nn` / `SCT_snn` are taken; Phase 2 will use
  `postint_harmony_nn` / `postint_harmony_snn`.
- All 93 metadata columns and their cardinalities are recorded in
  `results/phase2/handoff/handoff_metadata_inventory.tsv`.

---

## 5. Metadata relationships and candidate Harmony variables

Full analysis: **`reports/phase2/HARMONY_VARIABLE_DECISION.md`**. Summary:

| Concept | Field | Levels | Status |
| --- | --- | ---: | --- |
| Dataset | `sample_id` | 4 | present |
| Sample | `sample_id` | 4 | identical to dataset |
| Patient | — | — | **not recorded**; 1 sample = 1 patient |
| Technical batch | — | — | **not recorded** |
| Biological condition | — | — | **does not exist** |
| `orig.ident` | `orig.ident` | 4 | cell-for-cell identical to `sample_id` |

**Dataset = sample = patient = presumed technical batch is a single 4-level variable.**
Correcting technical batch is mathematically indistinguishable from erasing
between-patient tumour biology. Phase 2 §8 therefore requires this to STOP for
researcher review, which is what this milestone does.

Selected (proposed) grouping variable: **`sample_id`**, sole entry in `group.by.vars`.

---

## 6. Environment verification

| Component | Version | Required action |
| --- | --- | --- |
| R | 4.4.3 (`/local/projects-t3/lilab/vmenon/anaconda3/envs/R_env/bin/R`) | none |
| Seurat | 5.4.0 | none |
| SeuratObject | 5.3.0 | none |
| **harmony** | **1.2.4 — present** | none |
| Matrix | 1.7.4 | none |
| ggplot2 | 4.0.1 | none |
| patchwork | 1.3.2 | none |
| sctransform / glmGamPoi | 0.4.3 / 1.18.0 | none |
| presto | 1.0.0 | none (accelerates `FindAllMarkers` in M14) |
| ComplexHeatmap / pheatmap | 2.22.0 / 1.0.13 | none |
| aricode / mclust / cluster | 1.0.3 / 6.1.2 / 2.1.8.1 | none (clustering-stability metrics, M13) |
| `lisi` | **NOT INSTALLED** | see below |
| `kBET` | **NOT INSTALLED** | see below |
| `clustree` | **NOT INSTALLED** | see below |

Full records: `reports/phase2/PHASE2_ENVIRONMENT.tsv`, `reports/phase2/PHASE2_SESSIONINFO.txt`.

**Nothing was installed, upgraded or otherwise changed in the environment.** Per Phase 2
§3, the three missing packages are reported, not installed. Mitigation without them:

- **iLISI / cLISI** → the inverse Simpson index is computed natively from the exact kNN
  graph in M12, reusing the block-wise kNN routine already validated in Phase 1
  `analyze_pre_integration.R`. The deviation from the `lisi` package's Gaussian-kernel
  neighbourhood will be stated in `HARMONY_ASSESSMENT.md`.
- **kBET** → replaced by a per-cluster observed-vs-expected batch composition table plus
  a chi-square-style deviation statistic, reported descriptively.
- **clustree** → the M13 resolution-transition diagram is drawn directly with ggplot2
  from the cluster-membership transition counts.

If the researcher prefers the published implementations, installing `lisi`, `kBET` and
`clustree` is a separate action requiring explicit approval.

---

## 7. Phase 1 defects found (documentation only — nothing repaired)

Reported per Phase 2 §5 and §29. None of these invalidates integration; all are
inconsistencies between Phase 1 prose and Phase 1 machine-readable records, and the
machine-readable records are treated as authoritative.

| # | Location | Documented | Actual (authoritative source) |
| --- | --- | --- | --- |
| D1 | `PHASE1_HANDOFF.md` §2.2 | min_features 500, min_counts 1000, mt 15%, ribo 20% (global) | `config/config.yaml` + `phase1_manifest.json`: min_features 200, min_counts 500, and **dataset-specific** mt (10/15/10/20%) and ribo (20/30/35/30%) |
| D2 | `PHASE1_HANDOFF.md` §2.5 | "Max: 1–30" for all datasets | `reports/PCA_RECOMMENDATIONS.tsv`: max 24 / 25 / 21 / 22 |
| D3 | `PHASE1_HANDOFF.md` §2.6 | `n_resamples = 100` | `config/config.yaml: clustering.n_resamples = 5` |
| D4 | `PHASE1_HANDOFF.md` §2.3 | "scDblFinder (v2.20.2)" | `scDblFinder` 1.20.2 (v2.20.2 is not a real version of that package) |
| D5 | `PRE_INTEGRATION_ASSESSMENT.md` §3.2, `PHASE1_HANDOFF.md` §3.1 | "a single unified SCTransform globally" | four SCT models — Finding 1, §4.2 |
| D6 | `phase1_manifest.json` | `combined_preintegration.feature_count = 31764` | 31,764 is the `RNA` assay; the default `SCT` assay has 29,113 |
| D7 | `phase1_manifest.json` | `shared_pca_configuration.npcs = 30` | `RunPCA(npcs = 50)`; 30 is the number of PCs *used*, not computed |
| D8 | `config/config.yaml: input_rds` | `/local/projects-t3/lilab/vmenon/Pilot_tumor/processed_mpnst.rds` | that directory **no longer exists**; the file is `./processed_mpnst.rds` in the current project root. All `file://` links in the Phase 1 reports point at the same dead root. |
| D9 | `reports/PRE_INTEGRATION_ASSESSMENT.md` §3.4 vs `INTEGRATION_PREPARATION.md` Q1 | "mean same-dataset fraction >98%" (CHANGELOG) | 0.9613 in `neighborhood_mixing_summary.tsv`; the reports themselves say 96.13% — only the CHANGELOG entry overstates it |

D8 is the only one with an operational consequence (a Phase 1 re-run from the original
RDS would fail on a dead path). It does not affect Phase 2, which consumes the combined
object. Recommendation: fix D1–D9 as a documentation-only commit, either now or at the
Phase 2 freeze. **Awaiting researcher instruction; nothing has been changed.**

---

## 8. SLURM accounting — M10

| Field | Value |
| --- | --- |
| JobID | **19886411** |
| Job name | `p2_m10_inspect` |
| Account / Partition / Node | `ihc` / `ihc` / `ihc-grid-1-1-1` |
| Requested CPUs | 4 |
| Requested memory | 96G |
| Requested walltime | 02:00:00 |
| State | **COMPLETED** |
| ExitCode | **0:0** |
| Elapsed | 00:02:17 |
| MaxRSS | 4,959,508 K (**4.73 GiB**) |
| TotalCPU | 00:01:30 |
| Memory efficiency | 4.9% of request — **over-provisioned; corrected for M11** |
| CPU efficiency | ~33% (single-threaded R work) |
| Warnings / errors in log | **none** (`logs/phase2/slurm/m10_inspect_19886411.err` is empty) |

Breakdown from the log: md5 + sha256 of the 6.0 GB file took 72 s; `readRDS()` took 34 s;
all inspection and table writing took 3 s.

---

## 9. Proposed M11 SLURM resources

Evidence-based, from §8 and from Phase 1 M8 benchmarks
(`benchmarks/run_pre_integration_analysis.tsv`: 17.9 GB / 9m43s for the far heavier
SCTransform + PCA + UMAP pass).

| Resource | Proposed | Justification |
| --- | --- | --- |
| CPUs | **8** | Harmony default is `ncores = 1`; extra cores serve BLAS and object serialisation only. 8 matches the Phase 1 M8 precedent. |
| Memory | **64G** | Observed load footprint 4.7 GB; Harmony on a 19,716 × 30 matrix adds <1 GB; the object is then written back out (~6 GB serialised) with the original still resident. 64G is ~10× the observed peak and matches the Phase 1 M8 allocation. Not the 450G maximum. |
| Walltime | **04:00:00** | Load 34 s + Harmony (minutes at this scale) + a second matrix-level Harmony run for convergence diagnostics + save. Phase 1's heaviest comparable job was 9m43s. |
| GPU | none | Harmony is CPU-only. |
| Node | `ihc-grid-1-1-1` | Project standard. |

M12–M17 resources will be re-derived from M11's observed usage rather than assumed.

---

## 10. Files created and modified in M10

### Created

```
scripts/R/phase2/inspect_phase1_handoff.R
scripts/shell/phase2/run_m10_inspect.sh
reports/phase2/PHASE2_PLAN.md
reports/phase2/HARMONY_VARIABLE_DECISION.md
reports/phase2/PHASE2_ENVIRONMENT.tsv
reports/phase2/PHASE2_SESSIONINFO.txt
reports/phase2/milestones/M10_REPORT.md            (this file)
results/phase2/handoff/handoff_structure.json
results/phase2/handoff/handoff_assays.tsv
results/phase2/handoff/handoff_reductions.tsv
results/phase2/handoff/handoff_metadata_inventory.tsv
results/phase2/handoff/candidate_batch_crosstab.tsv
results/phase2/handoff/cells_per_sample.tsv
results/phase2/handoff/qc_summary_by_sample.tsv
results/phase2/handoff/legacy_orig_anno_by_sample_PROVENANCE_ONLY.tsv
logs/phase2/slurm/m10_inspect_19886411.out
logs/phase2/slurm/m10_inspect_19886411.err          (empty)
```

Directory skeleton created: `results/phase2/{harmony,clustering,markers,annotation,composition,handoff}`,
`reports/phase2/{milestones,figures}`, `scripts/R/phase2`, `scripts/shell/phase2`,
`logs/phase2/{slurm,milestones}`, `benchmarks/phase2`.

### Modified

```
PROGRESS.md      (Phase 2 section added; Phase 1 log untouched)
CHANGELOG.md     ([M10] entry appended)
```

### Deliberately NOT modified

`processed_mpnst.rds`, `results/datasets/**`, `results/combined/**`,
`results/phase1_manifest.json`, `reports/FIGURE_INDEX.tsv` (no Phase 2 figures exist
yet), `config/config.yaml`, `workflow/Snakefile`, and every Phase 1 report.

---

## 11. Scientific decisions taken in M10

1. **Input object**: `results/combined/pre_integration/combined_preintegration.rds`,
   verified against both Phase 1 checksums. It is the object Phase 1 designated as the
   handoff and the only one containing all four datasets in a shared space.
2. **Harmony grouping variable**: `sample_id` proposed (not executed) — see §5.
3. **Harmony dimensions**: `1:30` of `pca`, matching the frozen baseline exactly so that
   M12 compares integration rather than dimensionality.
4. **Namespace**: Phase 2 reductions/columns take the `postint_*` prefix required by the
   Phase 1 entry contract, combined with the Phase 2 specification's names
   (`postint_harmony`, `postint_umap_harmony`,
   `postint_harmony_clusters_res_0.1` … `_res_1.0`).
5. **Pre-Harmony baseline**: the existing `pca` and `umap_preintegration` reductions
   *are* the baseline. They are preserved unmodified; no duplicate `pca_pre_harmony`
   copy is created. Deviation from the suggested naming, documented in
   `PHASE2_PLAN.md` §3.1.
6. **No Snakemake in Phase 2**: modular R scripts + explicit `sbatch` scripts.
7. **Mixing metrics computed natively** rather than installing `lisi` / `kBET`.
8. **Legacy `orig.anno` is provenance-only** and is not used to select any parameter.
   Whether it may serve as a reference label for M12 biological-conservation diagnostics
   is deferred to the researcher (§12, item 4). Default if unanswered: not used.

---

## 12. Unresolved questions for the researcher

1. **Approve `sample_id` as the Harmony grouping variable**, accepting that
   batch/patient/tumour biology are inseparable, and choose Option A / B / C in
   `HARMONY_VARIABLE_DECISION.md` §7.1.
2. **Approve `dims.use = 1:30`.**
3. **Approve the `postint_*` naming reconciliation** (`PHASE2_PLAN.md` §3).
4. **May legacy `orig.anno` be used as a reference label for M12 biological-conservation
   diagnostics only** (never for parameter choice or Phase 2 annotation)? Default: no.
5. **Should Phase 1 documentation defects D1–D9 (§7) be corrected**, and if so now or at
   the Phase 2 freeze?
6. **Should `lisi` / `kBET` / `clustree` be installed?** Default: no; native
   implementations used.
7. **Is the four-SCT-model finding (§4.2) acceptable as-is**, with
   `PrepSCTFindMarkers()` in M14, rather than re-running M8 with a joined RNA assay?
   Recommended: accept as-is.

---

## 13. Warnings

- No R warnings or errors were raised by the M10 job; stderr is empty.
- The M10 job was over-provisioned on memory (4.9% efficiency). Corrected in the M11
  proposal.
- `.gitignore` already excludes `results/*`, `logs/`, `benchmarks/`, and `*.tsv`/`*.json`
  outside `reports/`, so no Phase 2 binary or bulk output is at risk of being committed.
  The new Phase 2 scripts and Markdown reports are the only commit-ready additions.

---

## 14. Next proposed milestone

**M11 — Default Harmony Integration.** Blocked on researcher approval of §12 items 1–3.

Planned M11 sequence:
1. Numerical-integrity precheck (Finding 2): non-finite values in `SCT@scale.data` and
   `pca` embeddings, per sample. Abort before Harmony if any are found.
2. Assert Phase 1 reductions `pca` and `umap_preintegration` are present and untouched.
3. `RunHarmony()` exactly as specified in `HARMONY_VARIABLE_DECISION.md` §8.
4. Matrix-level re-run with `return_object = TRUE` for the convergence objective history,
   asserted numerically identical to the primary embedding.
5. Save `results/phase2/harmony/harmony_integrated.rds` (a **new** file; the Phase 1
   object is never overwritten), plus `harmony_parameters.json`,
   `harmony_convergence.tsv`, and md5/sha256 checksums.
6. Inspect log, warnings and `sacct`; write `reports/phase2/milestones/M11_REPORT.md`;
   update `PROGRESS.md` and `CHANGELOG.md`.

---

## 15. STOP

> **STATUS: STOPPED.**
> M10 is complete. No Harmony has been run. Phase 2 will not proceed to M11 until the
> researcher explicitly approves the Harmony grouping variable, the dimensions and the
> naming contract in §12.
