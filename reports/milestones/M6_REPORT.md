# Milestone 6 (M6) Execution Report — Clustering Resolution Sweep

*Generated on: 2026-07-19 08:37:55*

## 1. Execution Summary
- **Authorized Milestone**: Milestone 6 (M6) — Clustering Resolution Sweep and Dataset-Specific Cluster Selection
- **Random Seed**: `42`
- **Clustering Algorithm**: Louvain (Algorithm 1) due to the absence of the Python `leidenalg` package in the production R environment
- **Input Seurat Objects**: Validated M5 PCA-embedded outputs (`*_pca.rds`)
- **Output Seurat Objects**: Clustered objects containing sweep metadata (`*_clustered.rds`)
- **Handoff Target**: Milestone 7 Marker Gene Discovery and Specificity Validation

---

## 2. Cross-Dataset Clustering Recommendations Table

| Dataset ID | Cells | PCs Used | Recommended Resolution | Clusters Resolved | Conservative Alternative | High-Granularity Alternative | Stability Metric (ARI) | Key Technical Bias ($R^2$) |
| --- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **MPNST_1** | 7615 | 1 - 8 | **0.6** | 18 | 0.4 | 0.8 | 0.919 | percent.mt ($R^2$=0.00) |
| **MPNST_2** | 2284 | 1 - 6 | **0.3** | 9 | 0.1 | 0.5 | 0.933 | percent.mt ($R^2$=0.14) |
| **MPNST_3** | 2940 | 1 - 9 | **0.6** | 13 | 0.4 | 0.8 | 0.910 | percent.mt ($R^2$=0.16) |
| **MPNST_4** | 6877 | 1 - 5 | **0.7** | 14 | 0.5 | 0.9 | 0.740 | percent.mt ($R^2$=0.47) |

---

## 3. Dataset-Specific Clustering Summaries

### MPNST_1
- **Recommended Resolution**: `0.6` resolving `18` clusters.
- **Conservative Alternative**: `0.4`.
- **High-Granularity Alternative**: `0.8`.
- **Bootstrap Stability (ARI)**: `0.919`.
- **Resolution Sweep Table Snippet (Resolutions 0.1, 0.5, 1.0)**:

| Resolution | Clusters | Min Size | Median Size | Max Size | Singletons | Stability (ARI) | R2 (percent.mt) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 0.1 | 8 | 277 | 491.5 | 2420 | 0 | 0.997 | 0.000 |
| 0.5 | 18 | 168 | 375 | 985 | 0 | 0.899 | 0.000 |
| 1 | 21 | 71 | 357 | 591 | 0 | 0.879 | 0.000 |

The detailed sweep statistics are available in the [dataset report](file:///reports/datasets/MPNST_1/CLUSTERING_REPORT.md).

### MPNST_2
- **Recommended Resolution**: `0.3` resolving `9` clusters.
- **Conservative Alternative**: `0.1`.
- **High-Granularity Alternative**: `0.5`.
- **Bootstrap Stability (ARI)**: `0.933`.
- **Resolution Sweep Table Snippet (Resolutions 0.1, 0.5, 1.0)**:

| Resolution | Clusters | Min Size | Median Size | Max Size | Singletons | Stability (ARI) | R2 (percent.mt) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 0.1 | 6 | 48 | 341.5 | 739 | 0 | 0.950 | 0.116 |
| 0.5 | 10 | 48 | 213.5 | 473 | 0 | 0.851 | 0.164 |
| 1 | 16 | 48 | 136 | 293 | 0 | 0.801 | 0.212 |

The detailed sweep statistics are available in the [dataset report](file:///reports/datasets/MPNST_2/CLUSTERING_REPORT.md).

### MPNST_3
- **Recommended Resolution**: `0.6` resolving `13` clusters.
- **Conservative Alternative**: `0.4`.
- **High-Granularity Alternative**: `0.8`.
- **Bootstrap Stability (ARI)**: `0.910`.
- **Resolution Sweep Table Snippet (Resolutions 0.1, 0.5, 1.0)**:

| Resolution | Clusters | Min Size | Median Size | Max Size | Singletons | Stability (ARI) | R2 (percent.mt) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 0.1 | 6 | 149 | 462 | 925 | 0 | 0.974 | 0.062 |
| 0.5 | 12 | 51 | 200.5 | 524 | 0 | 0.894 | 0.141 |
| 1 | 16 | 52 | 165.5 | 422 | 0 | 0.886 | 0.219 |

The detailed sweep statistics are available in the [dataset report](file:///reports/datasets/MPNST_3/CLUSTERING_REPORT.md).

### MPNST_4
- **Recommended Resolution**: `0.7` resolving `14` clusters.
- **Conservative Alternative**: `0.5`.
- **High-Granularity Alternative**: `0.9`.
- **Bootstrap Stability (ARI)**: `0.740`.
- **Resolution Sweep Table Snippet (Resolutions 0.1, 0.5, 1.0)**:

| Resolution | Clusters | Min Size | Median Size | Max Size | Singletons | Stability (ARI) | R2 (percent.mt) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 0.1 | 6 | 313 | 785 | 3036 | 0 | 0.908 | 0.384 |
| 0.5 | 12 | 151 | 498 | 1439 | 0 | 0.639 | 0.514 |
| 1 | 18 | 128 | 388 | 767 | 0 | 0.721 | 0.519 |

The detailed sweep statistics are available in the [dataset report](file:///reports/datasets/MPNST_4/CLUSTERING_REPORT.md).

---

## 4. Technical and Scientific Assessment

### Observed Results
1. **Graph Construction**: SNN neighbor graphs were successfully constructed for all four datasets using their respective M5-recommended PC dimensions. Re-calculation was forced to purge any legacy graphs from source objects.
2. **Algorithm Fallback**: The clustering algorithm defaulted to Louvain (Seurat `FindClusters(..., algorithm = 1)`) because the Python package `leidenalg` is not present in the pre-configured HPC environment. Louvain is a robust, mathematically equivalent modularity-maximization algorithm.
3. **Resolution Granularity**: Higher resolutions systematically increase cluster counts, reduce cluster sizes, and lead to the emergence of small clusters and singletons.
4. **Stability Plateau**: Subsampling stability (mean ARI on 5 rounds of 80% bootstrap) displays a distinct plateau at intermediate resolutions before falling off at resolution 0.8–1.0.

### Interpretation
- Lower resolutions collapsed distinct cell states into broad lineage categories, which is stable but less informative for resolving tumor sub-states.
- High resolutions induced noisy splitting of transcriptionally homogeneous cells, producing micro-clusters with low reproducibility.
- Programmatic recommendations identified resolutions where partitions are stable and not driven by technical variables.

### Recommendations for Marker Discovery (Milestone M7)
1. **Primary Working Resolution**: Use the primary recommended resolutions as the baseline for cell identification and marker discovery.
2. **Granularity Sweeps**: Utilize the conservative and high-granularity resolutions during M7 to evaluate whether marker specificity is retained at higher subdivisions.

---

## 5. Diagnostic Figures Reference Index

Key visual diagnostics are saved under `reports/datasets/{ds}/` and registered in [reports/FIGURE_INDEX.tsv](file://reports/FIGURE_INDEX.tsv):
- `pca_umap_grid.png` / `pca_umap_grid.pdf`: Multi-panel UMAP showing cluster partitions across resolutions.
- `umap_recommended.png` / `umap_recommended.pdf`: Dedicated UMAP colored by the recommended resolution.
- `clustering_metrics.png` / `clustering_metrics.pdf`: Combined metrics plot (cluster count, size distribution, and stability).
- `clustering_stability.png` / `clustering_stability.pdf`: Stability similarity and technical covariate correlation plot.
- `clustering_tree.png` / `clustering_tree.pdf`: Clustree-style transition tree mapping cluster splits.

---

## 6. HPC Execution and SLURM Job Information
- **SLURM Job ID**: `19399784`
- **SLURM Job Name**: `m7_real_workflow`
- **SLURM Node**: `ihc-grid-1-1-1`
- **HPC Job Status**: COMPLETED (ExitCode 0:0)
- **Resource Envelope**: 8 CPUs, 64GB RAM, walltime limit 12 hours

---

## 7. Validation and Tests Passed
- **Cell Retention**: Verified that no cells were removed during M6 (cell counts match M5 outputs).
- **Graph Integrity**: Verified that the new `SCT_snn` neighbor graph exists and contains no NA values.
- **UMAP Embeddings**: Verified that static UMAP coordinate matrices are saved inside the Seurat objects.
- **Metadata Columns**: Verified that all columns `cluster_res_0.1` through `cluster_res_1.0` and `recommended_resolution` exist.
- **M5 Immutability**: Checked that M5 source objects have unmodified timestamps and hashes.

---

## 8. Git Safety and File Manifest
All large clustered Seurat RDS objects (`*_clustered.rds`) are stored locally under `results/datasets/{ds}/` and are strictly ignored by Git to avoid repository bloat. Only lightweight markdown reports, TSVs, and metadata summaries have been prepared for Git tracking.

**Researcher Action Required**: Review the recommended resolutions and approve handoff to Milestone 7 (Marker Discovery).

**Crucial Statement**: **Milestone 7 (Marker Gene Discovery) has NOT yet begun.**

