# Milestone 8 Integration Preparation Report
*Generated on: 2026-07-19 23:25:04*

## 1. Integration Justification
Based on the complete spatial segregation observed in the pre-integration baseline, batch correction/integration is scientifically justified and required to perform joint downstream analyses (such as unified cell-type annotation and clustering).

## 2. Integration Design Decisions
1. **Is meaningful separation present?**: Yes, separation is global and complete.
2. **What variables represent technical effects?**: Library preparation chemistry, sequencing batch, and cell capture batch.
3. **What variables represent biological effects?**: Patient-specific tumor biology and cellular composition differences.
4. **What signal could aggressive integration remove?**: Patient-specific cell lineages or tumor-specific stress states (e.g. MPNST_4 stress state).
5. **Should the non-integrated baseline remain a permanent reference?**: Yes, to validate that post-integration alignments do not introduce artificial cell states or over-smooth biological boundaries.
6. **Proposed Phase 2 Methods**: Harmony, Seurat CCA, and Seurat RPCA should be benchmarked.
7. **Proposed Evaluation Metrics**: We recommend evaluating both batch mixing (e.g. LISI, kBET) and biological conservation (cell-type marker retention).
