# M50 — Phase 6 Object Build, Figures and Freeze

**Date** 2026-09-04 · **Status** COMPLETE · **SLURM** 19899749, 19899752

**Object.** `results/phase6/phase6_final_object.rds` — 6,060,757,871 bytes, md5
`b8c01dd01755dda10b2be4e0f5ef7ef7`, sha256
`bde592d461d22253b646da1426db8545404b9d31b7f1b948e23deed03aa643f6`. 19,716 cells, **241 metadata
columns** (196 Phase 1–5 + 46 Phase 6).

**What went in, and what deliberately did not.** The compact interpretable layer is embedded: the 14
PROGENy pathway activities, the 25 strongest cross-patient-concordant TF activities, and the
clone/program metrics broadcast to each clone's cells. The **full** TF (690 regulons) and Hallmark
(50 sets) matrices stay as separate RDS artefacts referenced by the manifest — pushing hundreds of
columns into object metadata would bloat it for no analytical gain.

**Preservation.** 8 guards asserted column by column against a copy taken before anything was added,
then 8 reload validations. All 16 passed, including that every Phase 5 `program_P*_score` column is
byte-identical and the frozen Phase 4 counts still hold. Phase 4 **and** Phase 5 md5s re-verified
unchanged.

**Figures.** 13/13 required, PDF + PNG. **Tables.** 12 in `tables/final/`. **Manifest.** 31 sections
with per-figure and per-table md5. **Figure index.** 706 → 732 rows.

**Visual QC.** Every figure rendered and inspected. One defect fixed: figure 13 panel B used a shared
colour scale across rows whose values span single digits to several hundred, flattening every row but
one and rendering dark-on-dark text; it now scales within row and prints the raw value.

**Environment.** `R_env_PRE_PHASE6.yml` and `R_env_POST_PHASE6.yml` byte-identical.

**Phase 7 is NOT initiated.**
