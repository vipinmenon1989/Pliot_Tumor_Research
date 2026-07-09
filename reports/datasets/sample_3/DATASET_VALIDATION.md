# Dataset Validation Report - sample_3
*Generated on: 2026-07-08 19:25:50*

## 1. Object Overview

| Metric | Value |
| --- | --- |
| **Dataset ID** | sample_3 |
| **Cells** | 150 |
| **Genes** | 120 |
| **Dimensions** | 120 x 150 |
| **Default Assay** | RNA |
| **Object Size in Memory** | 0.27 MB |
| **File Path** | `results/datasets/sample_3/sample_3_raw.rds` |
| **MD5 Checksum** | `92260c22b30e170a03a3b9885d11c140` |

## 2. Assay Structure

| Assay | Class | Features | Layers | Counts Present? |
| --- | --- | --- | --- | --- |
| RNA | Assay5 | 120 | counts | YES |

## 3. Metadata Columns

The following metadata columns are preserved in the Seurat object:

- `orig.ident`
- `nCount_RNA`
- `nFeature_RNA`
- `sample_id`
- `patient_id`
- `condition`
- `cell_type_annotation`
- `percent.mt`
- `percent.ribo`
- `unintegrated_clusters`
- `seurat_clusters`

## 4. Dimensional Reductions (Inventory Only)

The following legacy reductions were present in the original object and have been inventoried but removed from this working object:

No reductions found.

## 5. Provenance

- **Input File**: `data/synthetic/synthetic_mpnst.rds`
- **Input Checksum**: `0d6dd92216a9e3c4a71b28f538e19285`
- **Extraction Code**: `scripts/R/extract_dataset.R`
- **Git Commit**: `472429e70d564c230ce08a2fb9d36e52ccadf0f0`
