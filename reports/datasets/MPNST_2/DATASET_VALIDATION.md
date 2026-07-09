# Dataset Validation Report - MPNST_2
*Generated on: 2026-07-08 19:29:03*

## 1. Object Overview

| Metric | Value |
| --- | --- |
| **Dataset ID** | MPNST_2 |
| **Cells** | 2,830 |
| **Genes** | 29,708 |
| **Dimensions** | 29708 x 2830 |
| **Default Assay** | SCT |
| **Object Size in Memory** | 1216.47 MB |
| **File Path** | `results/datasets/MPNST_2/MPNST_2_raw.rds` |
| **MD5 Checksum** | `cd15a88b7f70a6492a1f05c94725cd8c` |

## 2. Assay Structure

| Assay | Class | Features | Layers | Counts Present? |
| --- | --- | --- | --- | --- |
| RNA | Assay5 | 31,764 | counts.MPNST_2, data.MPNST_2, scale.data | YES |
| mnn.reconstructed | Assay | 2,000 | data | NO |
| SCT | SCTAssay | 29,708 | counts, data, scale.data | YES |

## 3. Metadata Columns

The following metadata columns are preserved in the Seurat object:

- `orig.ident`
- `nCount_RNA`
- `nFeature_RNA`
- `sample_id`
- `percent.mt`
- `orig.anno`
- `unintegrated_clusters`
- `seurat_clusters`
- `pca100_harmony_clusters`
- `pca100_cca_clusters`
- `pca100_mnn_clusters`
- `pca100_rpca_clusters`
- `nCount_SCT`
- `nFeature_SCT`
- `pca100.sct_harmony_clusters`
- `pca100.sct_cca_clusters`
- `percent.ribo`

## 4. Dimensional Reductions (Inventory Only)

The following legacy reductions were present in the original object and have been inventoried but removed from this working object:

- `pca`
- `pca200`
- `pca100`
- `pca50`
- `umap.unintegrated`
- `pca100_cca`
- `umap.pca100_cca`
- `pca100_harmony`
- `umap.pca100_harmony`
- `pca100_mnn`
- `umap.pca100_mnn`
- `pca100_rpca`
- `umap.pca100_rpca`
- `pca100.sct`
- `pca100.sct_harmony`
- `umap.pca100.sct_harmony`
- `pca100.sct_cca`
- `umap.pca100.sct_cca`

## 5. Provenance

- **Input File**: `/local/projects-t3/lilab/vmenon/Pilot_tumor/processed_mpnst.rds`
- **Input Checksum**: `fe839def36e9246ecf4aa0d8df13ba18`
- **Extraction Code**: `scripts/R/extract_dataset.R`
- **Git Commit**: `472429e70d564c230ce08a2fb9d36e52ccadf0f0`
