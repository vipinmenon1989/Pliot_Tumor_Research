# Milestone M23 — Receiver-Response Analysis (NicheNet)

**Phase 3 · MPNST** · *2026-09-03* · Status: **COMPLETE** → proceeded to M25.
Full report: [`reports/phase3/RECEIVER_RESPONSE_REPORT.md`](../RECEIVER_RESPONSE_REPORT.md)

## Script
`scripts/phase3/receiver_response/run_nichenet.R` · nichenetr 2.2.1.1

## SLURM
| Job | JobID | State | Elapsed |
| --- | --- | --- | --- |
| Zenodo prefetch 1 | 19895781 | COMPLETED | 00:01:53 — 262 MB matrix failed (R `download.file` 60 s default timeout) |
| Zenodo prefetch 2 | 19895842 | COMPLETED | truncated at 112.7 of 262.1 MB |
| Zenodo prefetch 3 | **19895947** | COMPLETED | resumable `curl -C -` with size verification; 262,116,605 bytes, readable 33,354 × 1,226 |
| **M23** | **19896063** | **COMPLETED** | **00:03:13** |

## Design
Gene set of interest = **receiver marker genes** (cluster characterisation, **not** condition
DE). Background = genes in ≥10% of receiver cells present in the prior matrix. Potential
ligands = expressed in ≥10% of tumour cells with a receptor expressed in the receiver.

## Result — a clean split
**Myeloid receivers all converge on CSF1** (Macrophage 0.166, Monocyte 0.123, Dendritic 0.107).
**Lymphoid receivers all converge on IL15** (CD8-T 0.116, NK 0.089, CD4-T 0.084).
**Endothelial**: TGFB1, HMGB1, VEGFA, ANGPT1.
**Fibroblast**: best AUPR **0.022**, rest ≈ 0 — a clear negative: tumour ligands do not explain
the fibroblast programme.

## Divergence reported, not reconciled
**APP is not among the top NicheNet ligands for macrophages**, despite being the strongest LR
finding. The two methods ask different questions, so this is not a contradiction — but it
bounds the APP-CD74 claim to "predicted engagement", not "drives the macrophage state".

160 ligand × receiver receiver-response support flags emitted for M25.
