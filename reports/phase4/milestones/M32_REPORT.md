# M32 — Targeted Phase 3 CCC Sensitivity Analysis

**Status: COMPLETE** · 2026-09-03

**The question is sensitivity, not rediscovery**: does the Phase 3 biological architecture
survive a better tumour definition? The whole ~500,000-key interaction space was **not** rerun
blindly (§37) — what was rerun is exactly what a fair tumour-label comparison requires.

## Method reuse — the point of the design (§39, §40)

`scripts/phase3/ccc/run_liana.R`, `run_cellchat.R`, `run_cellphonedb.py` and
`scripts/phase3/concordance/build_concordance.R` were reused **completely unmodified**. Only
the `ccc_label` definition differs. That was possible because Phase 3's own scripts read a
generic `ccc_label` / `ccc_sample` contract, so a Phase 4 input-preparation script
(`prepare_refined_ccc_inputs.R`) was all that was needed.

Identical throughout: liana 0.1.14 (`resource = "Consensus"`, methods natmi/connectome/logfc/
sca/cellphonedb), CellChat 2.2.0.9001 (CellChatDB.human, 3,233 interactions, 290 pathways),
CellPhoneDB 5.0.1 in the isolated `cpdb_env` (1,000 iterations, threshold 0.10),
`min-cells = 10`, seed 42, RNA + JoinLayers + LogNormalize(1e4), per-sample analysis.
**Nothing was reinstalled or upgraded**, so the comparison is a label-sensitivity analysis
rather than a tool-drift confound.

## Refined input

`ccc_input_object_refined.rds` — **15,036 cells, 14 populations**, md5
`99fd641faf0cc1f7e38448406c883f32`. Excluded as not CCC-ready: `Ambiguous-unresolved` (3,766),
`Uncertain` (476), `Low-quality-excluded` (438). Excluded means **NOT EVALUABLE**, not "no
signalling".

Per-patient evaluable populations fell where refinement emptied a population: MPNST_1 14,
MPNST_2 14, MPNST_3 **10**, MPNST_4 **12**. `Fibroblast` retains 1 cell in MPNST_4 and
`Pericyte-VSMC` 4, so both are NOT EVALUABLE there. This is why the change classes distinguish
`Ambiguous` (no longer testable) from `Lost` (testable and no longer supported).

## Per-framework results on the refined labels

| Framework | rows | supported | tumour-involving | runtime |
| --- | ---: | ---: | ---: | ---: |
| LIANA | 164,271 | **31,661** | 7,010 (→TME 4,378, →tumour 3,213) | 33.0 min |
| CellPhoneDB | 1,208,748 | **21,844** | 6,292 | 1.9 min |
| CellChat | — | 5,412 LR + 2,190 pathway rows (MPNST_1) | — | 18–38 min/patient |

## Concordance — Phase 3 versus Phase 4

| | Phase 3 | Phase 4 refined | change |
| --- | ---: | ---: | ---: |
| distinct interaction keys | 502,846 | 494,973 | −1.6% |
| **supported** | **36,486** | **34,893** | **−4.4%** |
| **High concordance** | **5,848** | **5,698** | **−2.6%** |
| Moderate concordance | 499 | 513 | +2.8% |
| Single-method | 24,654 | 23,302 | −5.5% |
| Discordant/ambiguous | 5,485 | 5,380 | −1.9% |
| testable by all three | 1,347 | 1,297 | −3.7% |

**4,036 cells moved into the tumour compartment and 2,015 moved out, and high-concordance
interactions fell by only 2.6%.** That is the answer to the Phase 4 sensitivity question in one
line: the architecture is not an artefact of the Phase 2 tumour definition.

## Change classes (§42) — transparent categories, no invented score

Across 506,605 union keys:

| class | n | reading |
| --- | ---: | --- |
| Not-supported-either | 464,702 | never supported in either phase |
| **Stable** | **18,276** | unchanged |
| Weakened | 8,446 | fewer frameworks or fewer patients |
| Newly-supported | 5,421 | **but 91.8% rest on one patient — weak** |
| Lost | 3,963 | testable and no longer supported |
| **Sender-reassigned** | **3,051** | same LR pair and receiver; sender moved into `MPNST-Tumor` |
| Strengthened | 2,746 | more frameworks or more patients |

## The named axes (§38) — verdicts with magnitudes

`change_class` fires on a single-key difference, so a label alone is misleading: HLA-F→LILRB1
going 72 → 71 is technically "weakened". **Median change across the 23 axes is −3.1%, and not
one axis was lost.**

| behaviour | axes |
| --- | --- |
| **effectively unchanged** (\|Δ\| ≤ 5%) | APP→CD74 −3.1% · ANXA1→FPR1 0.0% · HLA-F→LILRB1 −1.4% · HLA-F→LILRB2 −3.1% · HLA-E→KLRC1 −2.0% · CD58→CD2 −1.4% · JAG1→NOTCH3 0.0% · JAG2→NOTCH2 0.0% · FGF2→FGFR1 −4.3% |
| **strengthened** | **FN1→ITGAV_ITGB8 +33.3%** · CD99→PILRA +14.0% · VEGFA→KDR +11.1% · JAG1→NOTCH2 +7.5% |
| moderately weakened (10–20%) | BAG6→NCR3 −13.0% · VEGFA→NRP1 −13.2% · JAG1→NOTCH4 −11.8% · COL1A2→ITGAV_ITGB8 −10.0% · THBS1→CD36 −8.6% |
| **genuinely weakened** (>20%) | **COL1A1→ITGAV_ITGB8 −44.4%** · SLIT2→ROBO1 −38.9% · DLL4→NOTCH2 −31.6% · COL6A2→ITGAV_ITGB8 −30.0% · VEGFA→FLT1 −28.0% |
| **sender reassigned** | all four ECM→integrin axes · VEGFA→KDR · JAG1→NOTCH3 · DLL4→NOTCH2 |

Axis-by-axis verdicts for §69: **APP–CD74 stable** (−3.1%, still 3 frameworks) ·
**HLA-E–KLRC1 stable** (−2.0%) · **HLA-F–LILRB1/2 stable** (−1.4% / −3.1%) ·
**VEGFA–KDR strengthened, VEGFA–FLT1 and VEGFA–NRP1 weakened** ·
**JAG–NOTCH stable to strengthened** (JAG1→NOTCH3 0.0%, JAG1→NOTCH2 +7.5%, JAG2→NOTCH2 0.0%;
DLL4→NOTCH2 −31.6% with sender reassignment).

## The fibroblast question (§43)

> Were some apparent fibroblast interactions actually tumour-to-tumour interactions caused by
> malignant Mes-NC-like cells being labelled fibroblasts?

**Partly, and the ECM–integrin axis is precisely where it happened.**

Of **9,483** Phase 3 supported interactions involving `Fibroblast`:

```text
6,547 (69.0%)  retained as Fibroblast   -- a real 908-cell fibroblast compartment remains
1,565 (16.5%)  REASSIGNED to the refined tumour compartment
1,371 (14.5%)  lost or no longer testable
```

The signature is unmistakable in the named axes: **all four collagen/FN1 → ITGAV_ITGB8 axes
changed sender**, and they are the most affected axes in the entire set — COL1A1 −44.4%,
COL6A2 −30.0%, COL1A2 −10.0%, while **FN1 rose +33.3%** because the refined tumour compartment
itself expresses FN1. The "ECM → tumour integrin" interaction Phase 3 read as
stroma-to-tumour is substantially **tumour-to-tumour (autocrine)**.

The 16.5% figure is deliberately conservative: it requires the identical LR pair *and* receiver
to reappear with `MPNST-Tumor` as sender. It counts interaction keys, not signalling mass.

## SLURM accounting

| JobID | Stage | State | Elapsed | ReqMem | MaxRSS |
| --- | --- | --- | --- | --- | --- |
| 19897175 | refined CCC input | COMPLETED | 00:07:04 | 96 G | 27.78 GiB |
| 19897176 | LIANA | COMPLETED | 00:34:36 | 64 G | 38.87 GiB |
| 19897177_1 | CellChat MPNST_1 | COMPLETED | 00:39:48 | 48 G | 18.33 GiB |
| 19897177_2 | CellChat MPNST_2 | COMPLETED | 00:20:28 | 48 G | 14.49 GiB |
| 19897177_3 | CellChat MPNST_3 | COMPLETED | 00:19:40 | 48 G | 14.58 GiB |
| 19897177_4 | CellChat MPNST_4 | COMPLETED | 00:23:32 | 48 G | 18.60 GiB |
| 19897178 | CellPhoneDB | COMPLETED | 00:02:31 | 48 G | 1.64 GiB |
| **19897179** | concordance + sensitivity | **FAILED** | 00:03:39 | 48 G | 1.34 GiB |
| 19897616 | sensitivity rerun | COMPLETED | ~00:01 | 48 G | — |

**Failure root-caused, not retried blindly.** 19897179's *concordance* step exited 0 and its
outputs are intact; only my sensitivity script failed, for two reasons: (a) `fib3` was derived
from `t3`, which has no `rk` column — `NULL %in% reasg` returns `logical(0)` and the assignment
aborted; (b) `canon()` strips punctuation, so an upstream resource recording the same complex
as both `"CD8 RECEPTOR"` and `"CD8_RECEPTOR"` collapsed to one key, leaving 4 duplicate keys per
table and a many-to-many join. Fixed by computing `rk` explicitly and by collapsing duplicate
keys while keeping the strongest support and logging how many were affected — rather than
silencing the warning. Only the failed script was rerun; the 96-minute CCC computation was not
repeated. **No resource was increased.**

CellChat was run as four concurrent per-patient jobs from the outset, applying the Phase 3
lesson where a serial run had to be cancelled for projected wall clock.

## Figures

`32_01_phase3_vs_phase4_ccc` · `32_02_axis_sensitivity`, plus the concordance figure suite
`build_concordance.R` writes into `results/phase4/figures/M32/`.

## Tables

`PHASE3_VS_PHASE4_CCC.tsv` (506,605 rows) · `PHASE3_AXIS_SENSITIVITY.tsv` ·
`FIBROBLAST_INTERACTION_REASSIGNMENT.tsv` · `CCC_LABEL_PHASE3_VS_PHASE4.tsv` ·
`CCC_REFINED_SAMPLE_CELL_COUNTS.tsv` · `results/phase4/ccc_refinement/concordance/CCC_CONCORDANCE.tsv`.

## Next milestone

**M33 — tumour-state → TME model**, then **M34 — robustness**, then **M35 — freeze**.
