# Milestone M19 — Sample-Aware CCC Input Preparation

**Phase 3 · MPNST** · *2026-09-03* · Status: **COMPLETE** → proceeded to M20.

## Scripts
`scripts/phase3/ccc/prepare_ccc_inputs.R`, `scripts/shell/phase3/run_m19_ccc_inputs.sh`

## SLURM
| JobID | State | ExitCode | Elapsed | CPUs | Mem | MaxRSS |
| --- | --- | --- | --- | ---: | ---: | ---: |
| **19894939** | COMPLETED | 0:0 | 00:08:24 | 8 | 160G | 28,765,372 K (27.43 GiB) |

No failures. R warnings: 0.

## Expression basis — a documented decision
CCC uses the **RNA assay, `JoinLayers()` + LogNormalize (scale.factor 1e4)**, not SCT.
The SCT assay carries **four SCTransform models**, so SCT values are not on a common footing
across samples, while every CCC framework expects one uniformly log-normalised matrix.
**Harmony coordinates are never used as expression.**

## Minimum-cell policy — stated, not silently chosen
1. Populations flagged not CCC-ready in Phase 2 are excluded outright:
   `Candidate-Malignant-Unresolved` (1,231), `Uncertain` (720), `Low-quality-excluded` (438).
2. A population is EVALUABLE in a sample only with **≥ 10 cells** — CellChat's `min.cells`
   default and the conventional floor for per-group mean-expression estimates.
3. A sample is ANALYSABLE with ≥ 3 evaluable populations.
4. **A population below the minimum is recorded NOT EVALUABLE, never as absence of signalling.**

## Result
**14 populations · 17,327 cells · all 4 samples analysable.** Every population is evaluable in
every sample, and **all 182 directed pairs are evaluable in all four samples** (13 tumour→TME,
13 TME→tumour, 156 TME↔TME). 41 of 42 probe ligand/receptor genes were present (NECTIN2 absent).

Outputs: `results/phase3/ccc/ccc_input_object.rds` (4.13 GB, md5 `aeb02c59df8a42a24074df509c662f7f`);
per-sample CellPhoneDB MTX inputs; `CCC_SAMPLE_CELL_COUNTS.tsv`,
`CCC_POPULATION_EVALUABILITY.tsv`, `CCC_PAIR_EVALUABILITY.tsv`, `CCC_INPUT_EXPRESSION_QC.tsv`;
4 figures × PDF + PNG.
