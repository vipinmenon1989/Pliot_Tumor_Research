# Milestone M14 — Marker Discovery

**Phase 2 · MPNST single-cell integration** · *Generated 2026-09-02*
*Status: **COMPLETE** — proceeded automatically to M15.*

> Full detail: [`reports/phase2/MARKER_REPORT.md`](../MARKER_REPORT.md)
> **These are cluster-characterisation markers, not condition-level differential expression.**

## 1. Scripts

Created: `scripts/R/phase2/discover_markers_harmony.R`, `scripts/shell/phase2/run_m14_markers.sh`,
`reports/phase2/MARKER_REPORT.md`, `reports/phase2/milestones/M14_REPORT.md`.

## 2. Command

```bash
sbatch scripts/shell/phase2/run_m14_markers.sh
#  -> Rscript scripts/R/phase2/discover_markers_harmony.R \
#       --input results/phase2/clustering/phase2_harmony_clustered.rds \
#       --out-dir results/phase2/markers --fig-dir results/phase2/figures/M14 \
#       --assay SCT --min-pct 0.25 --logfc 0.25 --random-seed 42 --expected-cells 19716
```

## 3. SLURM

| JobID | State | ExitCode | Elapsed | CPUs | Mem | MaxRSS | Efficiency |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| **19886699** | **COMPLETED** | **0:0** | **00:04:26** | 12 | 250G | **17,193,456 K (16.40 GiB)** | 6.6% memory, 0.6% walltime |

No failures, no retries. **R warnings: 0. R errors: none.** The 250G request was headroom
for `PrepSCTFindMarkers()`, which recorrects counts across four SCT models; actual peak was
16.4 GiB, so later milestones were sized down accordingly.

## 4. Seurat v5 assay/layer handling — inspected, not assumed

| Assay | Class | Default | Layers | SCT models |
| --- | --- | :---: | --- | ---: |
| `RNA` | `Assay5` | no | split four ways (`counts.MPNST_1.1` … `scale.data.4`) | — |
| `SCT` | `SCTAssay` | **yes** | `counts; data; scale.data` (joined) | **4** |

**`PrepSCTFindMarkers()` was required and was run**, recorrecting counts to a common
sequencing depth so cells normalised under different SCT models are comparable. The `RNA`
assay was not used (split layers would need `JoinLayers()` plus renormalisation, discarding
the Phase 1 SCT model). **Harmony coordinates were never used for testing** — Harmony is an
embedding-level method and carries no expression values.

## 5. Parameters

`assay = SCT`, `layer = data`, `test.use = wilcox` (presto-accelerated), `min.pct = 0.25`,
`logfc.threshold = 0.25`, `only.pos = TRUE`, `recorrect_umi = FALSE`, `random.seed = 42`.
These match Phase 1 `config/config.yaml: markers`, keeping M14 comparable to the Phase 1
per-sample runs. Identities: `postint_harmony_primary_cluster` (resolution 1.0, 26
clusters); repeated on the alternative resolution (0.7, 21 clusters) for robustness.

## 6. Results

- **35,437 marker rows** over 26 clusters; **31,774** significant (adj. p < 0.05,
  log2FC ≥ 0.25). **Every cluster returned markers** (min 212 in C15, max 2,901 in C18).
- Key lineage evidence for M15: **C9 retains MPZ** (canonical Schwann myelin protein) with
  GFRA3 and ABCB5; **C8 expresses L1CAM** (reported SCP-like malignant marker); **C14
  expresses SHH** with the Schwann-associated GAL3ST1 and KLK6.
- **C21's markers (LILRA4, CLEC4C, SPIB, GZMB) are decisive for plasmacytoid dendritic
  cells**, overriding the coarse M12 programme score that had called it B/plasma — pDCs
  share IGJ and MZB1 with plasma cells. Marker evidence corrected a programme-score prior.
- **C15 is dominated by mitochondrial transcripts** (a technical signature, matching the
  Phase 1 flag on MPNST_4 cluster C09, `percent.mt` R² = 0.47) and **C23 by ribosomal
  pseudogenes** — neither is a cell type.
- **C12 carries both a hypoxia programme and the reported perineurial signature**
  (SLC2A1/GLUT1, ITGB4, ITGA6, CAV1) — genuinely ambiguous, left `Uncertain` in M15.

## 7. Tables and figures

11 tables in `results/phase2/markers/` (see MARKER_REPORT §4) and 6 figures × PDF + PNG in
`results/phase2/figures/M14/`.

## 8. Limitations

Cluster boundaries inherit the M12 caveats (B/plasma and fibroblast cohesion each fell 44%
under Harmony), so markers for C2/C21/C22 and C0/C3/C5/C7/C15/C25 warrant a cross-check
against the non-integrated baseline. `only.pos = TRUE` means negative markers are not
tabulated. p-values on 19,716 cells are not a measure of biological significance; effect
sizes and detection rates carry the interpretation.

## 9. Next

Proceeded automatically to **M15 — Literature-Grounded Cell-Type Annotation**.
