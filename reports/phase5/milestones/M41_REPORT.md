# M41 — Phase 5 Object Build, Figures and Freeze

**Date** 2026-09-03/04 · **Status** COMPLETE · **SLURM** 19899359 (object), 19899742, 19899744,
19899747, 19899750

**Object.** `results/phase5/phase5_final_object.rds` — 6,058,471,006 bytes, md5
`839e5157bc7c3470ddf86746c2e719e1`, sha256
`fe99ecf51154046145cf21f9c6960664d10205b90736bb13c61cccf409c40f00`. 19,716 cells, 196 metadata
columns (183 Phase 1–4 + 13 Phase 5).

**Preservation.** 10 guards asserted **column by column** against a copy taken before any field was
added, then 8 reload validations after saving. All 18 passed. Phase 4 md5 re-verified unchanged.

**Projection.** Cells outside the malignant compartment carry program scores obtained by NNLS
against the **fixed** spectra and are flagged `program_score_source == "projected"` — they could not
influence a program, and no Phase 4 call was altered.

**Figures.** 12/12 required, PDF + PNG. **Tables.** 16 in `tables/final/`. **Manifest.** 30 sections
with per-figure and per-table md5. **Figure index.** 682 → 706 rows.

**Failures, root-caused and fixed at source.**
*19899359* — the figures script read UMAP coordinates from `meta.data`, but they live in the
object's reduction. Fixed by reading the frozen Phase 4 embedding table; the object build in the
same job had already passed all validations and was not repeated.
*19899742* — a new column named `sig` shadowed the `sig` data frame inside the same `mutate()`.
Renamed `is_sig`.

**Visual QC.** Every figure was rendered and inspected. Three defects were found and fixed:
facet panels ordered `10% / 2% / 5%` and mislabelling 0.025 as "2%"; overlapping subtitles in
figure 08; and a figure-11 verdict panel reading "NOT robust" eight times, which was replaced by the
two-verdict panel described in the robustness report.

**Environment.** `R_env_PRE_PHASE5.yml` and `R_env_POST_PHASE5.yml` byte-identical.

**Next.** Phase 6 (M42–M50).
