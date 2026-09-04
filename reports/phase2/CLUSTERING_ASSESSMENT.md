# Post-Harmony Clustering Assessment

**Phase 2 · Milestone M13 · MPNST single-cell integration**
*Generated: 2026-09-02 · SLURM JobID 19886690*

**Recommendation: PRIMARY resolution 1.0 (26 clusters) · ALTERNATIVE resolution 0.7 (21 clusters)**

---

## 1. Construction

| Step | Configuration |
| --- | --- |
| Input | `results/phase2/clustering/../harmony/phase2_harmony_integrated.rds` (M11, md5 `cf63e843…`) |
| Embedding | `postint_harmony`, dims **1:30** (the M12-accepted Harmony result) |
| Neighbours | `FindNeighbors(reduction = "postint_harmony", dims = 1:30, k.param = 20)` → graphs `postint_harmony_nn`, `postint_harmony_snn` |
| UMAP | `RunUMAP(reduction = "postint_harmony", dims = 1:30, seed.use = 42)`, all other parameters Seurat defaults (n.neighbors 30, min.dist 0.3, cosine) → reduction `postint_umap_harmony` |
| Clustering | `FindClusters(graph.name = "postint_harmony_snn", algorithm = 1 (Louvain), random.seed = 42)` |
| Resolutions | 0.1 … 1.0 in steps of 0.1 — **all ten preserved** as `postint_harmony_clusters_res_0.1` … `_res_1.0` |
| Cells | 19,716 (unchanged) |

`k.param = 20` and the Louvain algorithm match Phase 1 `config/config.yaml`, keeping the
post-integration sweep methodologically comparable to the Phase 1 per-sample sweeps.

**Phase 1 and M11 state preserved.** The `pca`, `umap_preintegration` and
`postint_harmony` embeddings, the `SCT_nn`/`SCT_snn` graphs, the assays and all 93
original metadata columns are byte-identical after clustering. `FindClusters()`
overwrites the legacy `seurat_clusters` column and the active identities as a side
effect; both were snapshotted beforehand and restored verbatim, and a per-column digest
check confirms every original column is unchanged in both value and R type.

---

## 2. Sweep results

| Resolution | Clusters | Min / median / max size | Tiny (<99 cells) | Sample-dominated (>60%) | Mean silhouette | Program coverage | Program purity |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 0.1 | 13 | 287 / 960 / 4210 | 0 | 4 | 0.2012 | 8 | 0.787 |
| 0.2 | 15 | 211 / 960 / 4000 | 0 | 5 | 0.1888 | 9 | 0.757 |
| 0.3 | 17 | 212 / 820 / 3998 | 0 | 6 | 0.1883 | 9 | 0.743 |
| 0.4 | 17 | 212 / 960 / 3990 | 0 | 6 | 0.1907 | 9 | 0.734 |
| 0.5 | 18 | 212 / 730 / 3511 | 0 | 8 | 0.1811 | 9 | 0.729 |
| 0.6 | 20 | 214 / 748 / 3504 | 0 | 9 | 0.1902 | 9 | 0.750 |
| **0.7** | **21** | 212 / 822 / 2369 | 0 | 9 | 0.2024 | 9 | 0.758 |
| 0.8 | 22 | 213 / 670 / 2290 | 0 | 9 | 0.1989 | 9 | 0.766 |
| 0.9 | 25 | 162 / 512 / 2194 | 0 | 10 | 0.2039 | 9 | 0.767 |
| **1.0** | **26** | 141 / 493 / 2192 | 0 | 11 | **0.2078** | 9 | 0.760 |

**No resolution produced a single tiny cluster.** The smallest cluster anywhere in the
sweep is 141 cells (0.7% of the data) at resolution 1.0. Fragmentation is therefore not a
discriminating criterion here, and the sweep contains no degenerate solution.

### Stability (adjacent-resolution ARI)

| Transition | ARI | NMI |
| --- | ---: | ---: |
| 0.1 → 0.2 | 0.908 | 0.930 |
| 0.2 → 0.3 | 0.933 | 0.942 |
| 0.3 → 0.4 | 0.970 | 0.977 |
| 0.4 → 0.5 | 0.946 | 0.963 |
| 0.5 → 0.6 | 0.920 | 0.934 |
| 0.6 → 0.7 | 0.894 | 0.956 |
| 0.7 → 0.8 | 0.976 | 0.976 |
| 0.8 → 0.9 | 0.892 | 0.932 |
| 0.9 → 1.0 | 0.970 | 0.977 |

**Every adjacent transition has ARI ≥ 0.89.** The partition is highly stable across the
entire range: moving from 13 to 26 clusters splits existing groups rather than
reorganising them. The practical consequence is that **the resolution choice is not
critical** — any solution in this sweep describes essentially the same structure at
different granularity.

---

## 3. Selection rule

Deterministic and coded, with UMAP appearance explicitly excluded. Each criterion is
min-max scaled across the ten resolutions and combined with fixed weights:

```
composite = 0.30 · stability (mean adjacent ARI)
          + 0.20 · compactness (mean silhouette in the Harmony space)
          + 0.15 · compartment coverage (distinct canonical programs that are modal in >=1 cluster)
          + 0.15 · program purity (mean per-cluster purity w.r.t. canonical programs)
          + 0.10 · low fragmentation (fewer tiny clusters)
          + 0.10 · granularity (penalty outside a 10-30 cluster band for annotation workability)
```

| Resolution | Composite | Rank |
| ---: | ---: | ---: |
| **1.0** | **0.831** | 1 (PRIMARY) |
| 0.9 | 0.634 | 2 |
| **0.7** | **0.619** | 3 (ALTERNATIVE) |
| 0.8 | 0.608 | 4 |
| 0.4 | 0.578 | 5 |
| 0.3 | 0.552 | 6 |
| 0.2 | 0.445 | 7 |
| 0.1 | 0.406 | 8 |
| 0.5 | 0.374 | 9 |
| 0.6 | 0.372 | 10 |

The alternative is the highest-ranked solution whose cluster count differs from the
primary by at least 2, so that the two offer genuinely different granularity (26 vs 21).

### 3.1 Caveat — the primary sits at the sweep boundary

Resolution 1.0 is the top of the tested range, which has two consequences that must be
stated:

1. **Its stability score is one-sided.** Interior resolutions average the ARI to the
   resolution below *and* above; 1.0 has only the transition from 0.9 (0.970), and 0.1
   only the transition to 0.2 (0.908). This inflates 1.0's stability term relative to
   interior resolutions and is the main reason its composite score leads by a wide
   margin.
2. **The optimum may lie beyond 1.0.** The sweep does not test whether 1.1 or 1.2 would
   continue to add well-populated, well-separated clusters.

Resolution 1.0 nevertheless has the best two-sided evidence available on the other
criteria — the highest silhouette (0.2078), full compartment coverage (9), no tiny
clusters, and a minimum cluster size of 141 — so it is retained as primary. The
alternative, **0.7**, is an interior solution validated on both sides (0.894 and 0.976)
and is the recommended fallback if the researcher prefers coarser, more conservatively
supported granularity. Given ARI ≥ 0.89 throughout, downstream annotation is unlikely to
change materially between them.

---

## 4. Sample composition and the M12 caveats

`sample_id` is simultaneously the dataset, the patient and the Harmony grouping variable,
so cluster composition by sample is also composition by dataset and by patient. There is
no biological-condition field.

At the primary resolution, **11 of 26 clusters are >60% one sample**. This is expected
rather than alarming: M12 established that Harmony deliberately left the
Schwann/neural-crest (presumptive malignant) compartment patient-private while mixing the
shared immune and vascular compartments. Sample-dominated clusters should be read as
candidate patient-private tumour populations, not as integration failures — but each must
be treated conservatively at annotation.

The two M12 caveats were carried forward as explicit machine-readable flags
(`results/phase2/clustering/primary_cluster_m12_caveat_flags.tsv`):

| Flag | Clusters | Meaning |
| --- | ---: | --- |
| **C1** — B/plasma program >20% of cluster | **3** (C2, C21, C22) | M12 found B/plasma cohesion fell 44% under Harmony. These clusters must be cross-checked against the non-integrated baseline before annotation. |
| **C2** — fibroblast program >20% of cluster | **6** (C0, C3, C5, C7, C15, C25) | M12 found fibroblast cohesion fell 44%. Same requirement. |
| Sample-dominated (>60%) | **11** | Candidate patient-private populations; annotate conservatively. |

Notable structure at the primary resolution (canonical program fractions, provisional):

| Cluster | n | Dominant sample | Signal |
| --- | ---: | --- | --- |
| C9 | 713 | MPNST_1 (79.8%) | SchwannNC 0.986 — strongest Schwann-lineage/tumour-like candidate |
| C14 | 440 | MPNST_1 (94.3%) | SchwannNC 0.993 |
| C8 | 967 | MPNST_1 (91.0%) | SchwannNC 0.799 |
| C21 | 287 | MPNST_3 (57.5%) | B/plasma 0.718 |
| C22 | 233 | MPNST_3 (39.9%) | B/plasma 0.403 |
| C5 | 1,313 | MPNST_4 (78.1%) | fibroblast 0.802 |
| C25 | 141 | MPNST_1 (74.5%) | fibroblast 0.638 |
| C0 | 2,192 | MPNST_4 (55.5%) | fibroblast 0.552 |

These are canonical-program sanity-check fractions from M12, **not** annotations.
Annotation is M15 and is evidence-driven from the M14 markers.

---

## 5. Outputs

Object: `results/phase2/clustering/phase2_harmony_clustered.rds` (5.64 GB, md5
`b91f0eecd73e202804ac7c6859db4c46`). New metadata: the ten
`postint_harmony_clusters_res_*` columns, `postint_harmony_primary_cluster`,
`postint_harmony_alternative_cluster`, `postint_primary_resolution`,
`postint_alternative_resolution`, ten `m12_progscore_*` columns and
`m12_program_argmax`. New reduction `postint_umap_harmony`; new graphs
`postint_harmony_nn` / `postint_harmony_snn`.

Tables (`results/phase2/clustering/`): `clustering_comparison_table.tsv`,
`cluster_sizes_all_resolutions.tsv`, `cluster_sizes_primary.tsv`,
`cluster_composition_by_sample_all_resolutions.tsv`, `cluster_composition_primary.tsv`,
`resolution_transition_ari.tsv`, `primary_cluster_m12_caveat_flags.tsv`,
`m13_clustering_record.json`, `prov_m13_clustering.json`, `figure_index_m13.tsv`.

Figures (`results/phase2/figures/M13/`, 22 figures × PDF + PNG): the labelled primary
cluster UMAP, the alternative cluster UMAP, UMAPs by sample and by canonical program, one
UMAP per resolution, a ten-panel sweep grid, cluster size and proportion barplots,
composition-by-sample stacked bars (percentage and count), the adjacent-resolution ARI
stability curve, the primary↔alternative transition matrix, and the four-panel selection
criteria figure.

---

## 6. Limitations

1. The primary resolution sits at the sweep boundary (§3.1); its stability term is
   one-sided and the sweep does not test beyond 1.0.
2. Cluster compactness is modest in absolute terms (silhouette ≈ 0.18–0.21 across the
   whole sweep). Single-cell transcriptomes rarely form well-separated spheres, so this is
   ordinary, but it means silhouette discriminates weakly between resolutions here.
3. The canonical-program coverage and purity criteria depend on the M12 sanity-check
   program labels, which leave 26.3% of cells `Unassigned`.
4. Eleven sample-dominated clusters cannot be distinguished, from clustering alone,
   between genuine patient-private tumour biology and residual uncorrected batch effect —
   the permanent consequence of dataset = patient = batch.
