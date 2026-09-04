# M36 — Phase 5 Reconstruction and Feasibility

**Date** 2026-09-03 · **Status** COMPLETE · **SLURM** 19899335

**Scientific question.** Can the frozen Phase 4 malignant compartment support a search for
*continuous* transcriptional programs, and which frozen field carries the pre-Phase-4 biological
annotation that M39's fibroblast comparison must be grouped by?

**Inputs.** `results/phase4/phase4_final_object.rds`, opened read-only after verification by md5
**and** sha256.

**Outputs.** `reports/phase5/PHASE5_FEASIBILITY.md` ·
`results/phase5/validation/{m36_feasibility_facts.json, M36_HISTORICAL_ANNOTATION_FIELD_AUDIT.tsv, M36_REQUIRED_FIELD_AUDIT.tsv}` ·
`results/phase5/tables/{M36_MALIGNANT_PATIENT_DISTRIBUTION.tsv, M36_PHASE4_STATE_BY_PATIENT.tsv}` ·
five working artefacts under `results/phase5/programs/`.

**Scripts.** `scripts/phase5/utils/{phase5_common.R, m36_feasibility.R}` ·
`scripts/shell/phase5/run_m36_feasibility.sh`

**Resources.** COMPLETED, Elapsed 00:02:05, AllocCPUS 4, ReqMem 96 G, **MaxRSS 29.28 GiB**.
The 96 G request was set above M32's observed 27.78 GiB peak on the same object, not at the maximum.

**Results.**
* 13 frozen Phase 4 counts re-asserted, all passed.
* The historical annotation field is **`annotation_ccc_phase3`** — *verified* by reproducing the
  frozen fibroblast split 5,064 → 4,036 / 908 / 120, not assumed. `annotation_ccc_refined` was
  excluded by design: with only 908 fibroblasts, all Non-malignant, it visibly encodes the Phase 4
  conclusion.
* Patient imbalance is **17.4×** (MPNST_4 3,685 = 57.3%; MPNST_3 212 = 3.3%) and is confounded with
  depth (MPNST_4 median 5,250 counts vs MPNST_2's 18,901).
* `tumor_clone_phase4` is missing for 178 malignant cells — recorded, not hidden.
* Each Mesenchymal_ECM state is effectively one patient's (ECM-1 = MPNST_1, ECM-4 = MPNST_2,
  ECM-2/-3/-5 = MPNST_4).

**Warnings.** MPNST_3 contributes 212 malignant cells and **0 High-confidence** ones; it can enter
program discovery but cannot carry a per-patient conclusion.

**Failures / fixes.** None.

**Scientific decisions.** Working artefacts were exported so no later milestone reopens the 6 GB
object. The working cell set is malignant ∪ historical Fibroblast (7,462 cells) so M39 needs no
second load. The Phase 4 object md5 was re-checked after M36 and is unchanged.

**Next.** M37 — continuous program discovery.
