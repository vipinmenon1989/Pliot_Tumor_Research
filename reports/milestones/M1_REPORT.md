# Milestone 1 (M1) Execution Report — Extended Real-Data Audit
*Generated on: 2026-07-07 15:30:00 EDT*

## 1. Execution Summary
- **Authorized Milestone**: Milestone 1 (M1) - Real-Data Audit
- **Input File audited**: `/local/projects-t3/lilab/vmenon/Pilot_tumor/processed_mpnst.rds`
- **Object Class**: Seurat
- **Dimensions**: 29708 features x 22661 cells
- **Memory Loading Time**: 45.35 seconds
- **File size**: 7.23 GB
- **Memory footprint in R**: 10.08 GB

## 2. Environment Details
- **R version**: R version 4.4.3 (2025-02-28)
- **Seurat version**: 5.1.0
- **Platform**: x86_64-conda-linux-gnu

## 3. Generated Artifacts
- [DATA_AUDIT.md](file://reports/DATA_AUDIT.md)
- [DATASET_STRUCTURE.tsv](file://reports/DATASET_STRUCTURE.tsv)
- [object_inventory.json](file://reports/object_inventory.json)
- [object_inventory.tsv](file://reports/object_inventory.tsv)

## 4. Scientific Findings & Recommendations
1. **Constituent Dataset Identifier**: Recommend using `sample_id` as the primary dataset splitter. Both `sample_id` and `orig.ident` are 100% identical in mappings, but `sample_id` represents standard biological naming conventions.
2. **Reprocessing Feasibility**: **FEASIBLE**. The object contains raw RNA counts in split count layers, allowing complete independent reprocessing from baseline transcripts.
3. **Prior Analysis Inventory**: Stored legacy clustering (e.g. seurat_clusters, harmony_clusters) and reductions (e.g. pca, umap.unintegrated, pca100_harmony, pca100_cca) have been inventoried. They are legacy artifacts and will NOT influence downstream Phase 1 processing.
4. **Feature Statistics**: Zero duplicated genes or duplicated cell names were detected. 13 mitochondrial genes matching `^MT-` and 433 ribosomal genes matching `^RP[SL]` are present.

## 5. Execution Metrics
- **Start Time**: 2026-07-07 14:45:42
- **End Time**: 2026-07-07 14:46:28
- **Elapsed audit time**: 45.35 seconds
- **Peak memory (MaxRSS)**: 9627268K (~9.18 GB)
- **CPU time**: mean_load 54.72%, cpu_time 42.23s

## 6. Verification Status
- Snakemake execution: Successful
- Configuration validation: Successful
- Unit tests: Passed

## 7. Researcher Decisions Required
Please refer to Section 12 of [DATA_AUDIT.md](file://reports/DATA_AUDIT.md) for the complete list of decisions that require approval before proceeding to Milestone 2 (M2).
