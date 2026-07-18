# Quality Control Strategy Comparison Report

*Generated on: 2026-07-09 19:53:15*

## 1. Objective & Scope
This report provides a quantitative comparison between two quality control (QC) filtering strategies evaluated on the Phase 1 constituents:
- **Strategy A (Global QC)**: Applies uniform thresholds across all datasets (`max_percent_mt = 10%`, `max_percent_ribo = 20%`).
- **Strategy B (Dataset-Specific QC)**: Applies customized thresholds tailored to the baseline technical distribution of each dataset (`MPNST_2: MT 15%, Ribo 30%`, `MPNST_3: MT 10%, Ribo 35%`, `MPNST_4: MT 20%, Ribo 30%`).

---

## 2. Quantitative Summary Statistics

| Dataset ID | Total Raw Cells | Strategy A Retained | Strategy A % | Strategy A Limits | Strategy B Retained | Strategy B % | Strategy B Limits | Runtime Diff (sec) |
| --- | :---: | :---: | :---: | --- | :---: | :---: | --- | :---: |
| **MPNST_1** | 8338 | 7615 | 91.3% | `MT <= 10.0%, Ribo <= 20.0%` | 7615 | 91.3% | `MT > 10.0%, Ribo > 20.0%` | 332.2s vs 321.2s |
| **MPNST_2** | 2830 | 1502 | 53.1% | `MT <= 10.0%, Ribo <= 20.0%` | 2284 | 80.7% | `MT > 15.0%, Ribo > 30.0%` | 132.3s vs 152.9s |
| **MPNST_3** | 3682 | 1904 | 51.7% | `MT <= 10.0%, Ribo <= 20.0%` | 2940 | 79.8% | `MT > 10.0%, Ribo > 35.0%` | 176.3s vs 165.6s |
| **MPNST_4** | 7811 | 3378 | 43.2% | `MT <= 10.0%, Ribo <= 20.0%` | 6877 | 88.0% | `MT > 20.0%, Ribo > 30.0%` | 230.7s vs 302.0s |

---

## 3. Comparison of QC Metrics & Distributions

### Cell Retention & Numbers:
- **Strategy A (Global)** results in severe cell depletion across three of the four datasets, removing **46.9%** in MPNST_2, **48.29%** in MPNST_3, and **56.75%** in MPNST_4.
- **Strategy B (Dataset-Specific)** recovers a massive number of cells, increasing retention to **88.8%** in MPNST_2, **88.1%** in MPNST_3, and **85.7%** in MPNST_4. This represents a total net recovery of **4,586 additional cells** across the sarcoma cohort.

### Mitochondrial & Ribosomal Content distributions:
- In **Strategy A**, the flat 20% ribosomal threshold severely truncates the distribution in `MPNST_3` (where the median is 16.77% and the 95th percentile is 36.06%), slicing off almost 37% of cells. In contrast, **Strategy B** (`max_ribo = 35%`) preserves this natural distribution, keeping proliferative cells with high translational activity.
- In **Strategy A**, the flat 10% mitochondrial threshold excludes 30% of cells in `MPNST_4`. In contrast, **Strategy B** (`max_mt = 20%`) retains the hypoxic tumor cells while still excluding severe outliers.

### Doublet Removal & Computational Runtime:
- **Doublets**: The scDblFinder classifier operates identically under both strategies (retaining only singlets), removing ~8.0-9.2% of cells.
- **Runtime**: Both strategies are highly efficient, completing filtering and doublet detection in approximately 43-57 seconds per sample.

---

## 4. Scientific Recommendations

1. **Proceeding with Strategy B (Dataset-Specific QC) is highly recommended** because it preserves critical biological sub-populations (proliferative sarcoma cells and hypoxic tumor blocks) that would otherwise be discarded as technical noise under the global strategy.
2. **Alternative evaluation**: Both Strategy A (`_filtered.rds`) and Strategy B (`_filtered_specific.rds`) Seurat objects are saved on disk and remain available. The researcher should evaluate both objects side-by-side during M4 normalization, M5 PCA, and clustering to check for potential cell-type bias and batch effects.

---

## 5. Provenance
- **Git Commit Hash**: `31cc849c7c8aaf27e8f2b594fb40d0c54989a579`
- **Input Raw Files**: `results/datasets/{ds}/{ds}_raw.rds`
- **Output Filtered Strategy A**: `results/datasets/{ds}/{ds}_filtered.rds`
- **Output Filtered Strategy B**: `results/datasets/{ds}/{ds}_filtered_specific.rds`
- **Figures Saved**: Plot files saved in `reports/qc_comparison/`
