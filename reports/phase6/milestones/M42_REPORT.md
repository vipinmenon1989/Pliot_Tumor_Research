# M42 — Phase 6 Reconstruction and Feasibility

**Date** 2026-09-03 · **Status** COMPLETE · **SLURM** 19899748 (00:02:42, 8 CPUs, 128 G, MaxRSS 21.92 GiB)

**Question.** Does the Phase 5 object support the clone/regulatory analyses, and does the frozen
Phase 4 reliability flag independently reproduce the MPNST_3 exclusion?

**Result.** Yes to both. Phase 5 object verified by md5 and sha256 against its manifest; all required
fields present; the frozen `scevan_sample_reliable` flag marks MPNST_3 and only MPNST_3 unreliable.
18 of 19 clones in the reliable patients clear the 20-cell minimum; `MPNST_4_clone8` (5 cells) is
reported NOT EVALUABLE. 178 malignant cells carry no clone label and are excluded from clone
analyses with the number stated.

**Outputs.** `reports/phase6/PHASE6_FEASIBILITY.md` ·
`results/phase6/validation/{m42_feasibility_facts.json, M42_REQUIRED_FIELD_AUDIT.tsv, M42_CLONE_RELIABILITY_BY_PATIENT.tsv, M42_CLONE_SIZES.tsv}` ·
`results/phase6/clone_program/{phase6_cell_metadata.rds, phase6_malignant_lognorm.rds}`

**Failures.** One, in the first submission (19899746): `m43_clone_program.R` grouped by a column
`clone` that exists only after renaming `tumor_clone_phase4`. Root-caused and fixed at source;
M42 itself had already completed with exit 0 in that job.

**Next.** M43 — clone ↔ program coupling.
