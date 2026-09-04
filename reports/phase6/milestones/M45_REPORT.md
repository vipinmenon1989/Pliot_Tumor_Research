# M45 — Regulatory Architecture (TF Activity)

**Date** 2026-09-03/04 · **Status** COMPLETE · **SLURM** 19899749

**Method.** decoupleR 2.12.0 `run_ulm` over **CollecTRI**, on the frozen RNA log-normalised layer —
never Harmony, never UMAP, never the SCEVAN CNA matrix. 6,434 malignant cells; genes detected in
≥ 1% of them; TFs with ≥ 5 targets. **pySCENIC was not introduced**: it would move the frozen stack
for a question decoupleR already answers.

**A dependency defect, worked around without moving a frozen package.**
`decoupleR::get_collectri()` and `OmnipathR::collectri()` both **fail** in the installed stack —
OmnipathR 3.14.0's `unnest_evidences()` errors on the CollecTRI static table. Upgrading OmnipathR to
fix a convenience wrapper would have destabilised a frozen Phase 3/4 dependency. Instead the same
data was taken from **OmniPath's documented REST endpoint**
(`https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&organisms=9606`), cached
with its provenance in `external/networks/`, and converted to the standard source/target/mor form
(mor = −1 only for purely repressive edges). Result: **41,674 edges, 1,201 TFs**.

**Result — the strongest cross-layer confirmation in Phase 6.** Top positive regulators per program,
with `*` marking direction concordant in ≥ 3 evaluable patients **and** |pooled ρ| ≥ 0.20:

| program | positive | negative |
| --- | --- | --- |
| P1 Translation_ribosomal | HDAC5 0.70, PITX3 0.67, ZKSCAN7 0.65* | CREB1 −0.72*, NR1H3 −0.71*, BRCA1 −0.70* |
| P2 Neuronal | GATA3 0.59*, HOXA9 0.54*, RORC 0.53* | ZBTB4 −0.65*, MYC −0.64*, ZBED1 −0.64* |
| P3 **Mesenchymal_ECM** | ZBTB4 0.58*, HMGA2 0.57*, MYC 0.56*, SP1 0.56*, STAT6 0.54* | HDAC3 −0.45*, MBD2 −0.44* |
| P4 Hypoxia_Angio | ID4 0.45*, PGR 0.44*, TBX2 0.43* | ZBTB4 −0.54*, GBX2 −0.52*, MYC −0.50* |
| P5 Translation_ribosomal | PBX1 0.48*, HSF1 0.45*, **HIF1A 0.45***, **ATF4 0.43*** | FOXA3 −0.46*, LHX1 −0.44* |
| P6 Schwann_like | HIVEP2 0.32*, HSF2 0.32*, ATF3 0.31* | HOXA9 −0.27*, NR0B1 −0.24* |
| P7 **Cycling** | **E2F4 0.43***, **MYC 0.34***, **E2F1 0.24*** | OTX2 −0.29*, ARNT −0.28* |
| P8 Myeloid_ambient | ZNF148 0.21*, NKX6-1 0.20* | — |

**E2F4 / E2F1 / MYC leading the Cycling program, and HIF1A / ATF4 / HSF1 leading the
translation-stress program, are textbook and were not imposed** — the phase brief listed AP-1,
TEAD/YAP, STAT/IRF, E2F, MYC, NF-kB and SOX-family as hypotheses only, and the analysis recovered
E2F and MYC without being told to.

**A caution the reader needs.** 305 of ~690 TFs pass the "concordant in ≥ 3 patients" filter for P1.
With four patients, three agreeing in sign is not a stringent bar (≈ 31% by chance alone), so the
filter is weak here and the count should not be read as 305 validated regulators. The **direction
and magnitude** of the named TFs, not the count, is the evidence.

**Outputs.** `PROGRAM_TF_ACTIVITY.tsv` (5,520 rows) ·
`results/phase6/regulatory/m45_tf_activity_matrix.rds` (6,434 × 690).

**Next.** M46 — pathway activity.
