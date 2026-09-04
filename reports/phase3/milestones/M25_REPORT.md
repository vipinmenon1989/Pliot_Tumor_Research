# Milestone M25 — CCC + LochNESS Integration

**Phase 3 · MPNST** · *2026-09-03* · Status: **COMPLETE** → proceeded to M26.

## Script
`scripts/phase3/concordance/integrate_evidence.R` · SLURM **19896065**, COMPLETED, 00:00:24.

## Design — an evidence matrix, not a composite score
The streams are not commensurable, so they are displayed side by side and **counted**, never
weighted-averaged (Phase 3 §25, §38):

`≥2 LR frameworks` · `≥3/4 patients` · `ligand+receptor expression support` ·
`NicheNet receiver response`

| Evidence streams | Interactions |
| --- | ---: |
| 4 | **325** |
| 3 | 3,406 |
| 2 | 7,365 |
| 1 | 14,518 |
| 0 | 10,872 |

**547 tumour-centric interactions reach ≥3 streams.**

**LochNESS is attached as receiver-LINEAGE-level context only** — it is not interaction-specific
and is never treated as ligand-receptor evidence. Its values (Macrophage −0.163, Fibroblast
+0.104, Endothelial −0.176, CD8-T −0.311) were all non-significant, so it adds no positive
support to any interaction.

## The four-stream core (tumour-centric)
VEGFA→KDR/FLT1/NRP1 (Endothelial) · HLA-E→KLRC1 (NK) · HLA-F→LILRB1/LILRB2 (Monocyte) ·
HLA-A/HLA-E→CD8A/CD8B (CD8-T) · FGF2→FGFR1 (Fibroblast) · COL1A2→ITGA1/2/3/9_ITGB1 (Endothelial) ·
NCAM1→FGFR1, ANGPTL4→CDH5, CD58→CD2.

**APP→CD74 is absent from the four-stream set** — 3 frameworks, 4/4 patients and expression
support, but no NicheNet receiver-response evidence.

Outputs: updated `MPNST_CCC_MASTER_TABLE.tsv`, `CCC_INTEGRATED_EVIDENCE_TOP.tsv`,
3 figures × PDF + PNG.
