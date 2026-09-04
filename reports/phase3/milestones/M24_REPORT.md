# Milestone M24 — LochNESS Adaptation

**Phase 3 · MPNST** · *2026-09-03* · Status: **COMPLETE** → proceeded to M25.
Design (written first): [`LOCHNESS_MPNST_DESIGN.md`](../LOCHNESS_MPNST_DESIGN.md) ·
Results: [`LOCHNESS_MPNST_REPORT.md`](../LOCHNESS_MPNST_REPORT.md) ·
Comparison: [`LOCHNESS_IMPLEMENTATION_COMPARISON.md`](../LOCHNESS_IMPLEMENTATION_COMPARISON.md)

## Scripts
`scripts/phase3/lochness/lochness_mpnst.R`, `compare_implementations.R`, `lochness_python_port.py`

## SLURM
**19896064**, COMPLETED, ExitCode 0:0, Elapsed 00:04:14. No failures, 0 R warnings.

## Implementation
Official MMCA formulation: `k = round(0.5·√N)`, L2-normalised `postint_harmony` per lineage,
exact `FNN::get.knnx`, **disjoint query/reference sets (same-sample exclusion)**, global
fraction over the reference set. Two nulls: exact sample-level permutation (6 labelings,
**p-floor 0.167**) and 100 cell-level shuffles.

## Implementation comparison — decisive
| Python mode | Pearson | max abs diff |
| --- | ---: | ---: |
| MMCA-equivalent inputs | **1.0000** | **0.000e+00** |
| perturb-seq formulation as written | 0.3764 | 1.737 |

The two implementations are **numerically identical** given equivalent inputs. The perturb-seq
formulation diverges because it omits same-sample exclusion — which matters here because
sample = patient = batch. R was chosen on integration and fidelity grounds.

## Biological result — a negative
| Lineage | mean lochNESS | Null A p | Harmony-vs-PCA ρ | LOSO ρ |
| --- | ---: | ---: | ---: | ---: |
| Macrophage | −0.163 | 0.667 | 0.419 | 0.497 |
| Fibroblast | +0.104 | 0.833 | 0.023 | 0.347 |
| Endothelial | −0.176 | 0.333 | 0.219 | 0.448 |
| CD8-T | −0.311 | 0.333 | −0.265 | 0.310 |

**No lineage shows receiver-state structure associated with the tumour-derived APP context
beyond chance.** The score is additionally unstable across embeddings and reference samples.
Reported as a legitimate negative. **This does not refute the APP-CD74 LR finding** — LochNESS
is not a communication method.
