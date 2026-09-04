# Milestone M20 — Multi-Method CCC Execution

**Phase 3 · MPNST** · *2026-09-03* · Status: **COMPLETE** → proceeded to M21.

Three independent LR frameworks, each run **per sample**.

## Scripts
`scripts/phase3/ccc/run_liana.R`, `run_cellchat.R`, `run_cellphonedb.py`;
`scripts/shell/phase3/run_m20_{liana,cellchat_one,cellphonedb}.sh`

## SLURM accounting, including failures

| Job | JobID | State | ExitCode | Elapsed | MaxRSS | Note |
| --- | --- | --- | --- | --- | ---: | --- |
| LIANA (attempt 1) | 19895133 | **FAILED** | 1:0 | 00:01:38 | — | liana 0.1.14 calls `GetAssayData(slot=)`, **defunct** in SeuratObject 5.3.0 |
| LIANA | **19895159** | COMPLETED | 0:0 | **00:38:22** | — | |
| CellPhoneDB (attempt 1) | 19895344 | **FAILED** | — | 00:00:10 | — | `download_database()` keyword is `cpdb_version`, not `version` |
| CellPhoneDB | **19895438** | COMPLETED | 0:0 | **00:02:14** | — | |
| CellChat (serial) | 19895564 | **CANCELLED** | — | 00:44 (44% of sample 1) | — | serial ETA ≈ 5 h; deliberately re-run in parallel |
| CellChat MPNST_1 | **19895738** | COMPLETED | 0:0 | 00:44:13 | — | |
| CellChat MPNST_2 | **19895739** | COMPLETED | 0:0 | 00:19:49 | — | |
| CellChat MPNST_3 | **19895740** | COMPLETED | 0:0 | 00:22:51 | — | |
| CellChat MPNST_4 | **19895741** | COMPLETED | 0:0 | 00:24:19 | — | |

**Root causes and fixes**

1. **LIANA.** liana 0.1.14's Seurat code path uses the defunct `GetAssayData(slot=)`.
   **Fix: route through `SingleCellExperiment`**, which uses liana's SCE path. SeuratObject was
   **not** downgraded — that would have destabilised the frozen Phase 2 stack.
2. **CellPhoneDB.** API signature change in 5.0.1. Fixed the keyword.
3. **CellChat.** Not an error — a wall-clock decision. At 44% of sample 1 after 33 minutes the
   serial job projected ≈ 5 hours. The samples are scientifically independent after the shared
   M19 preparation, so it was cancelled and re-run as **four concurrent per-sample jobs**
   (Phase 3 §53), finishing in 44 minutes with identical results.

## Results

| Framework | Version | Rows pooled | Supported | Tumour-involving |
| --- | --- | ---: | ---: | ---: |
| LIANA | 0.1.14 | 192,015 | 36,684 | 6,280 |
| CellChat | 2.2.0.9001 | 17,506 (140 pathways) | 17,506 (significant only) | — |
| CellPhoneDB | 5.0.1 | 1,460,984 | 25,394 | 5,327 |

**Declared partial non-independence:** LIANA's method set includes a CellPhoneDB-style score.

Outputs in `results/phase3/ccc/by_method/` (full tables gzipped; supported tables plain).
