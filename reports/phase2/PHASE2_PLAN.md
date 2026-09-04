# MPNST Phase 2 — Modular Implementation Plan

*Created: 2026-09-02 (Milestone M10)*
*Status: DRAFT — awaiting researcher approval of M10*

---

## 1. Scope

Phase 2 extends the frozen Phase 1 pre-integration baseline into a Harmony-based
integration, post-integration clustering, marker discovery and literature-grounded
cell-type annotation workflow.

**Snakemake is not used in Phase 2.** Every production step is a standalone,
argument-driven R script under `scripts/R/phase2/` invoked by an explicit SLURM
batch script under `scripts/shell/phase2/`. This is a deliberate change from
Phase 1 and is recorded in `CHANGELOG.md`.

---

## 2. Input contract

| Item | Value |
| --- | --- |
| Phase 1 handoff object | `results/combined/pre_integration/combined_preintegration.rds` |
| Machine-readable input manifest | `results/phase1_manifest.json` |
| Cells | 19,716 |
| Constituent datasets | `MPNST_1`, `MPNST_2`, `MPNST_3`, `MPNST_4` |
| Expression assay | `SCT` (single global SCTransform on merged raw counts, `vars.to.regress = percent.mt`, 3,000 variable features) |
| Raw counts assay | `RNA` |
| Pre-Harmony PCA | reduction `pca` (50 PCs computed; 30 used for the Phase 1 neighbor graph and UMAP) |
| Pre-Harmony UMAP | reduction `umap_preintegration` |
| Random seed | 42 (inherited from Phase 1 `config/config.yaml`) |

Phase 1 objects are **read-only** for the whole of Phase 2. No Phase 2 script writes
to `results/datasets/**`, `results/combined/**`, `results/phase1_manifest.json`, or
`processed_mpnst.rds`.

---

## 3. Namespace contract

Phase 1 `reports/PHASE1_HANDOFF.md` §4 mandates a `postint_*` namespace for all
Phase 2 coordinates and metadata. The Phase 2 specification (§10, §13) asks for names
such as `harmony`, `umap_harmony`, `harmony_clusters_res_0.1`. These are reconciled by
prefixing the Phase 2 specification names with `postint_`:

### 3.1 Reductions

| Name | Meaning | Origin |
| --- | --- | --- |
| `pca` | Pre-Harmony shared PCA (the pre-Harmony baseline) | Phase 1 M8 — **preserved unmodified** |
| `umap_preintegration` | Pre-Harmony shared UMAP (the pre-Harmony baseline) | Phase 1 M8 — **preserved unmodified** |
| `postint_harmony` | Harmony-corrected embedding | Phase 2 M11 |
| `postint_umap_harmony` | UMAP computed on `postint_harmony` | Phase 2 M13 |

The Phase 1 `pca` / `umap_preintegration` reductions *are* the pre-Harmony baseline
required by Phase 2 §10; they are neither recomputed nor overwritten, so no
`pca_pre_harmony` duplicate is created. This deviation from the suggested naming is
deliberate, avoids a redundant multi-gigabyte copy of identical coordinates, and is
documented here and in `results/phase2/phase2_manifest.json`.

### 3.2 Graphs

| Name | Meaning |
| --- | --- |
| `SCT_nn`, `SCT_snn` | Phase 1 pre-integration graphs — **preserved unmodified** |
| `postint_harmony_nn`, `postint_harmony_snn` | Phase 2 neighbor graphs built on `postint_harmony` |

### 3.3 Metadata columns

| Name | Meaning | Milestone |
| --- | --- | --- |
| `preint_*` | All Phase 1 columns — **immutable** | M1–M8 |
| `postint_harmony_clusters_res_0.1` … `_res_1.0` | Every evaluated Harmony clustering resolution, all preserved | M13 |
| `postint_harmony_primary_cluster` | Researcher-approved primary resolution | M13/M14 |
| `postint_harmony_alternative_cluster` | Researcher-approved alternative resolution | M13/M14 |
| `postint_celltype_level1` | Broad compartment | M15 |
| `postint_celltype_level2` | Canonical cell type | M15 |
| `postint_celltype_level3` | State / subtype (only where supported) | M15 |
| `postint_annotation_confidence` | `high` / `medium` / `low` / `uncertain` | M15 |
| `postint_celltype_*_v2` | Post-review refinements; the M15 versions are retained | M16 |

---

## 4. Milestone module map

| Milestone | R script | SLURM script | Primary outputs |
| --- | --- | --- | --- |
| M10 | `scripts/R/phase2/inspect_phase1_handoff.R` | `scripts/shell/phase2/run_m10_inspect.sh` | `results/phase2/handoff/*`, `reports/phase2/HARMONY_VARIABLE_DECISION.md`, `reports/phase2/milestones/M10_REPORT.md` |
| M11 | `scripts/R/phase2/run_harmony_integration.R` | `scripts/shell/phase2/run_m11_harmony.sh` | `results/phase2/harmony/harmony_integrated.rds`, `results/phase2/harmony/harmony_parameters.json`, convergence diagnostics |
| M12 | `scripts/R/phase2/evaluate_harmony.R` | `scripts/shell/phase2/run_m12_evaluate.sh` | paired pre/post figures, `results/phase2/harmony/mixing_metrics_*.tsv`, `reports/phase2/HARMONY_ASSESSMENT.md` |
| M13 | `scripts/R/phase2/cluster_sweep_harmony.R` | `scripts/shell/phase2/run_m13_clustering.sh` | `results/phase2/clustering/*`, `reports/phase2/CLUSTERING_ASSESSMENT.md` |
| M14 | `scripts/R/phase2/discover_markers_harmony.R` | `scripts/shell/phase2/run_m14_markers.sh` | `results/phase2/markers/*`, `reports/phase2/MARKER_REPORT.md` |
| M15 | `scripts/R/phase2/annotate_celltypes.R` | `scripts/shell/phase2/run_m15_annotation.sh` | `results/phase2/annotation/ANNOTATION_EVIDENCE.tsv`, `reports/phase2/ANNOTATION_REPORT.md` |
| M16 | `scripts/R/phase2/refine_and_compose.R` | `scripts/shell/phase2/run_m16_composition.sh` | `results/phase2/composition/*` |
| M17 | `scripts/R/phase2/validate_and_freeze.R` | `scripts/shell/phase2/run_m17_freeze.sh` | `results/phase2/phase2_manifest.json`, `reports/phase2/PHASE2_HANDOFF.md` |

Every script:
- parses arguments explicitly (no hard-coded absolute paths beyond the project root),
- sets `set.seed(42)`,
- writes a provenance JSON via the Phase 1 `scripts/R/provenance_utils.R` helpers
  (md5 checksums, git commit, git status, `sessionInfo()`, SLURM JobID),
- exits non-zero on failure so the SLURM job is marked FAILED rather than silently
  producing partial output.

---

## 5. Quantitative integration diagnostics (M12) — implementation note

`lisi` and `kBET` are **not installed** in `R_env`. Per Phase 2 §3, packages will not be
installed without researcher approval. The M12 diagnostics will therefore be computed
natively in R from the k-nearest-neighbour graph, reusing and extending the block-wise
kNN routine already validated in Phase 1 `scripts/R/analyze_pre_integration.R`:

- **same-batch neighbour fraction** (directly comparable to the Phase 1 value of 0.9613),
- **neighbourhood Shannon entropy** (Phase 1 value 0.0966),
- **inverse Simpson index (iLISI-equivalent)** on the batch label — the same estimator
  used by the `lisi` package, computed on exact kNN rather than Gaussian-kernel
  weighted neighbourhoods; the deviation from the published estimator will be stated,
- **batch-composition deviation per cluster** (observed vs. expected dataset fractions),
- **silhouette width** of the batch label and of the biological label in both the
  pre-Harmony and post-Harmony spaces,
- **local structure preservation**: correlation between pre- and post-Harmony
  cell–cell distances, and retention rate of each cell's pre-Harmony k nearest
  neighbours after integration (an over-correction alarm).

If the researcher prefers the published `lisi`/`kBET` implementations, installing them
is a separate, explicitly approved action.

---

## 6. Explicit Phase 2 prohibitions carried into every script

No pseudobulk or condition-level DE, no differential abundance testing, no CNV
inference, no trajectory / velocity / cell–cell communication, no survival or
predictive modelling. `FindAllMarkers()` output is used for cluster characterisation
and annotation only, and every marker table and report carries that statement in its
header.

---

## 7. Open items requiring researcher decision at the M10 gate

1. Approval of `sample_id` as the Harmony grouping variable given that dataset =
   sample = patient = batch are perfectly confounded (see
   `reports/phase2/HARMONY_VARIABLE_DECISION.md`).
2. Approval of Harmony dimensions `1:30` (matching the frozen Phase 1 baseline).
3. Approval of the `postint_*` naming reconciliation in §3 above.
4. Whether the legacy `orig.anno` labels may be used as a *reference label* for
   biological-conservation diagnostics in M12 (never for parameter selection or for
   the Phase 2 annotation itself). Default if no answer: **not used**.
5. Whether any Phase 1 documentation defects listed in `M10_REPORT.md` §7 should be
   corrected before M11.
