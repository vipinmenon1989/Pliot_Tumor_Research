# Milestone M13 — Post-Harmony Neighbours, UMAP and Clustering Sweep

**Phase 2 · MPNST single-cell integration** · *Generated 2026-09-02*
*Status: **COMPLETE** — proceeded automatically to M14 under the continuous-execution authorization.*

> **Recommendation: PRIMARY resolution 1.0 (26 clusters) · ALTERNATIVE resolution 0.7 (21 clusters)**
> Full argument: [`reports/phase2/CLUSTERING_ASSESSMENT.md`](../CLUSTERING_ASSESSMENT.md)

---

## 1. Scripts

**Created:** `scripts/R/phase2/cluster_sweep_harmony.R`, `scripts/shell/phase2/run_m13_clustering.sh`
**Created (reports):** `reports/phase2/CLUSTERING_ASSESSMENT.md`, `reports/phase2/milestones/M13_REPORT.md`

## 2. Command

```bash
sbatch scripts/shell/phase2/run_m13_clustering.sh
#  -> Rscript scripts/R/phase2/cluster_sweep_harmony.R \
#       --input results/phase2/harmony/phase2_harmony_integrated.rds \
#       --out-rds results/phase2/clustering/phase2_harmony_clustered.rds \
#       --out-dir results/phase2/clustering --fig-dir results/phase2/figures/M13 \
#       --dims 30 --k-param 20 --random-seed 42 --expected-cells 19716
```

## 3. SLURM accounting — including three documented failures

| JobID | State | ExitCode | Elapsed | MaxRSS | Outcome |
| --- | --- | --- | --- | ---: | --- |
| 19886683 | FAILED | 1:0 | 00:01:49 | 8.49 GiB | **Guard fired correctly.** A whole-frame metadata digest flagged a change after `FindClusters()`. |
| 19886685 | FAILED | 1:0 | 00:02:20 | 9.88 GiB | `dom_sample[[k]]` subscript error — an unnamed vector indexed by cluster name. |
| 19886687 | FAILED | 1:0 | 00:02:34 | 15.3 GiB | `sprintf("%d", median(...))` — `median()` returns numeric, not integer. |
| **19886690** | **COMPLETED** | **0:0** | **00:09:37** | **12,679,688 K (12.09 GiB)** | Success. |

Requested for all four: 8 CPUs, 128G, 08:00:00 on `ihc-grid-1-1-1`. Final memory
efficiency 9.9%, walltime efficiency 2.0%. **No resource was increased in response to any
failure** — all three were code defects, diagnosed from the logs and fixed at the source.

### Root causes and corrections

1. **19886683.** The preservation check compared a digest of the *entire* metadata frame
   before and after clustering. It reported a change, but a whole-frame digest also reacts
   to attribute-level noise that carries no data. The check was replaced with a
   **per-column** comparison that names the offending column and separates genuine value
   changes from type/attribute changes, plus verbatim restoration of any original column
   Seurat touched. On the successful run this reported
   `original_meta_values_unchanged: TRUE`, `original_meta_types_unchanged: TRUE`,
   `original_metadata_restored_verbatim: TRUE` — i.e. no original value or type had in
   fact changed, and the stricter check now proves it column by column. This is a
   *strengthening* of the guard, not a relaxation.
2. **19886685 / 19886687.** Ordinary R defects (unnamed vector indexing; integer format
   specifier applied to a numeric). Scientifically neutral; no parameter changed.

## 4. Configuration

| Item | Value |
| --- | --- |
| Input | `results/phase2/harmony/phase2_harmony_integrated.rds` (md5 `cf63e84313a91de25fe2e41660f78f06`, verified) |
| Embedding | `postint_harmony` dims 1:30 |
| `FindNeighbors` | `k.param = 20`, graphs `postint_harmony_nn` / `postint_harmony_snn` |
| `RunUMAP` | dims 1:30, `seed.use = 42`, Seurat defaults otherwise → `postint_umap_harmony` |
| `FindClusters` | Louvain (algorithm 1), `random.seed = 42`, resolutions 0.1–1.0 |
| Cells | 19,716 (in = out) |

## 5. Outputs

**Object:** `results/phase2/clustering/phase2_harmony_clustered.rds` — 5.64 GB, md5
`b91f0eecd73e202804ac7c6859db4c46`.

New metadata: `postint_harmony_clusters_res_0.1` … `_res_1.0` (all ten preserved),
`postint_harmony_primary_cluster`, `postint_harmony_alternative_cluster`,
`postint_primary_resolution`, `postint_alternative_resolution`, ten `m12_progscore_*`
columns and `m12_program_argmax`.

**Tables (10)** in `results/phase2/clustering/` and **44 figure files (22 figures ×
PDF + PNG)** in `results/phase2/figures/M13/` — see `CLUSTERING_ASSESSMENT.md` §5.

## 6. Scientific findings

- **13 → 26 clusters** across resolutions 0.1 → 1.0. **No tiny cluster at any
  resolution**; the smallest anywhere is 141 cells (0.7%).
- **Every adjacent-resolution ARI ≥ 0.89** (range 0.892–0.976). Increasing resolution
  splits existing groups rather than reorganising them, so the resolution choice is not
  critical for downstream annotation.
- **Primary = 1.0 (26 clusters)** by a documented composite rule (0.30 stability + 0.20
  compactness + 0.15 compartment coverage + 0.15 program purity + 0.10 low fragmentation +
  0.10 granularity). UMAP appearance is not a criterion.
- **Caveat carried into the assessment:** resolution 1.0 is at the sweep boundary, so its
  stability term is one-sided (only the 0.9 → 1.0 transition), which inflates its
  composite score. It still leads on the two-sided criteria (highest silhouette 0.2078,
  full program coverage, no tiny clusters). **Alternative = 0.7 (21 clusters)** is an
  interior solution validated on both sides and is the conservative fallback.
- **11 of 26 primary clusters are >60% one sample.** Expected: M12 showed Harmony
  deliberately left the Schwann/neural-crest compartment patient-private. These are
  candidate patient-private tumour populations, to be annotated conservatively.
- **M12 caveats carried forward as machine-readable flags:** 3 clusters flagged C1
  (B/plasma > 20%: C2, C21, C22) and 6 flagged C2 (fibroblast > 20%: C0, C3, C5, C7, C15,
  C25). Each must be cross-checked against the non-integrated baseline before annotation.
- Strongest Schwann-lineage/tumour-like candidates: **C9** (713 cells, 79.8% MPNST_1,
  SchwannNC fraction 0.986), **C14** (440, 94.3% MPNST_1, 0.993), **C8** (967, 91.0%
  MPNST_1, 0.799).

## 7. Preservation

All Phase 1 and M11 state verified byte-identical after clustering, both in memory and
after a disk round-trip: `pca`, `umap_preintegration`, `postint_harmony`, the `SCT_nn` /
`SCT_snn` graphs, both assays, the default assay, cell names, and every one of the 93
original metadata columns (value **and** R type). `FindClusters()` overwrites the legacy
`seurat_clusters` column and the active identities as a side effect; both were
snapshotted and restored verbatim, and the round-trip check confirms
`legacy_seurat_clusters_intact: TRUE`.

## 8. Warnings

R warnings on the successful run: **0**. Errors: none. stderr contains only the standard
Seurat `RunUMAP` backend message. All 13 round-trip validation checks returned TRUE,
including `no_nonfinite_umap`.

## 9. Next

Proceeded automatically to **M14 — Marker Discovery** on the primary resolution (1.0,
26 clusters) with the alternative (0.7) computed for robustness.
