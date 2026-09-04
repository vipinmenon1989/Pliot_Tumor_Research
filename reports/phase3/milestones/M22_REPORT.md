# Milestone M22 — Biological Prioritisation

**Phase 3 · MPNST** · *2026-09-03* · Status: **COMPLETE** → proceeded to M23/M24.
Full report: [`reports/phase3/MPNST_INTERACTION_REPORT.md`](../MPNST_INTERACTION_REPORT.md)

## Script
`scripts/phase3/ccc/prioritize_interactions.R` · SLURM **19896062**, COMPLETED, 00:05:42.

## Prioritisation — lexicographic tiers, no invented composite
P1_highest 1,587 · P2_high 1,435 · P3_moderate 163 · P4_moderate_low 1,513 · P5_low 8,697 ·
P6_exploratory 23,091. Tumour-centric: 7,080 supported, **297 tumour→TME** and
**224 TME→tumour** in P1–P3.

## Headline findings
**Tumour → myeloid (3 frameworks, 4/4 patients):** APP→CD74 (5 receivers), CD99→PILRA,
ANXA1→FPR1, THBS1→CD36, HLA-F→LILRB1/2.
**Tumour → lymphoid:** HLA-A/HLA-E→CD8A/CD8B, HLA-E→KLRC1 (inhibitory), BAG6→NCR3 (activating),
CD58→CD2 — a genuinely **mixed** picture.
**Tumour ↔ vasculature/stroma:** VEGFA→KDR/FLT1/NRP1, JAG1→NOTCH3 (pericyte) / NOTCH4
(endothelial), FGF2→FGFR1, SLIT2→ROBO1, COL1A2→integrins.
**TME → tumour:** JAG1/JAG2/DLL4→NOTCH2 from pericytes and endothelium, CRTAM→CADM1 (CD8-T),
TNF→TNFRSF1A (CD4-T), collagens/FN1→ITGAV_ITGB8.

## LochNESS context label emitted for M24
Prioritised tumour ligand **APP**; sample-level split from **tumour cells only**:
MPNST_1 (1.081) and MPNST_4 (0.985) **high**; MPNST_2 (0.682) and MPNST_3 (0.666) low.

## Literature
12 prioritised axes curated with DOIs in `results/phase3/tables/CCC_LITERATURE_EVIDENCE.tsv`.
NOTCH and CSF1R are established in MPNST/NF1; SPP1-CD44 is documented in NF1 tumours;
APP-CD74 is established in other cancers but not characterised in MPNST. **No novelty is
claimed from a failed search.**
