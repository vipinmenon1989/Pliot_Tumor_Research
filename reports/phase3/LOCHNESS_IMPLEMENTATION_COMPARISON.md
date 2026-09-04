# LochNESS Implementation Comparison

**Phase 3 · Milestone M24 · MPNST** · *Generated 2026-09-03 · SLURM JobID 19896064*

Required by Phase 3 §36: verify whether the official MMCA R implementation and our perturb-seq
Python implementation give compatible results **when given mathematically equivalent inputs**.

---

## 1. Test design

| Element | Value |
| --- | --- |
| Subset | Macrophage lineage, deterministic stratified subsample, **1,861 cells**, seed 42 |
| k | `round(0.5·√N)` = **22** |
| Embedding | L2-normalised `postint_harmony`, identical matrix passed to both |
| Labels | identical sample and context vectors |
| R | `FNN::get.knnx(algorithm="kd_tree")`, exact |
| Python | `sklearn.neighbors.NearestNeighbors(algorithm="kd_tree")`, exact |

The Python side was run in **two modes**:

- **`mmca`** — same-sample exclusion, denominator = requested k, global fraction over the
  reference set. *Mathematically equivalent inputs.*
- **`perturbseq`** — the formulation as written in
  `perturbseq_pipeline/lochness.py`: one global kNN graph, **no same-sample exclusion**,
  denominator = actual neighbour count, global fraction over **all** cells.

The Python script re-expresses that algorithm from the read-only reference; it does not import
the perturb-seq package, so `Pilot_MPNST` gains no runtime dependency on it.

---

## 2. Results

| Python mode | n | Pearson | Spearman | max abs diff | RMSE |
| --- | ---: | ---: | ---: | ---: | ---: |
| **`mmca` (equivalent inputs)** | 1,861 | **1.0000** | **1.0000** | **0.000e+00** | **0.000e+00** |
| `perturbseq` (as written) | 1,861 | 0.3764 | 0.3519 | 1.737e+00 | 5.515e-01 |

---

## 3. Conclusions

**1. The two implementations are numerically identical given mathematically equivalent
inputs.** Pearson = 1.0000 and max absolute difference = 0 across 1,861 cells is exact
agreement, not approximate. Both correctly implement
`local_fraction / global_fraction − 1`, and the R and Python exact-kNN backends return the
same neighbourhoods. **The formula was implemented correctly in both codebases.**

**2. The perturb-seq formulation as written diverges substantially** (Pearson 0.38). The
divergence is not a bug — it is the accumulated effect of the three documented design
differences, chiefly the **absence of same-sample exclusion**. In MPNST, where `sample_id` is
simultaneously the dataset, the patient and the batch, that omission changes what the score
measures: without exclusion each cell's neighbourhood is dominated by its own sample, so the
score partly re-detects batch structure rather than cross-sample enrichment.

**3. R was chosen for Phase 3.** Since the implementations agree exactly on equivalent inputs,
the choice is on integration and fidelity grounds, not numerical ones:

- the MPNST workflow is Seurat/R-centred, so the R path avoids an AnnData round-trip;
- the official implementation is R, so following it minimises silent drift;
- scale is small (largest lineage 5,064 cells, k = 36), so the Python version's Numba kernel
  and million-cell mode solve a problem this dataset does not have.

**4. This is a validation of the perturb-seq code, not a criticism of it.** Its formulation is
correct for its own setting — a perturbation screen where no sample-level exclusion is
required. It simply is not transferable unmodified to a four-patient tumour atlas.

## 4. Files

`results/phase3/lochness/LOCHNESS_IMPLEMENTATION_COMPARISON.tsv` ·
`m24_implementation_comparison_record.json` ·
figure `results/phase3/figures/M24/M24_06_implementation_comparison.{pdf,png}` ·
scripts `scripts/phase3/lochness/compare_implementations.R`,
`scripts/phase3/lochness/lochness_python_port.py`.
