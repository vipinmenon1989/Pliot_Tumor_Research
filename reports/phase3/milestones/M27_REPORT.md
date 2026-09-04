# Milestone M27 — Phase 3 Final Freeze and Handoff

**Phase 3 · MPNST** · *2026-09-03* · Status: **COMPLETE — PHASE 3 FROZEN**

Handoff: [`reports/phase3/PHASE3_HANDOFF.md`](../PHASE3_HANDOFF.md) ·
Manifest: `results/phase3/phase3_manifest.json`

## 1. Validation performed

| Check | Result |
| --- | --- |
| All 10 milestone reports M18–M27 exist | ✅ |
| All 10 Phase 3 reports exist | ✅ |
| `phase3_manifest.json` written (31 sections) | ✅ |
| Final figure suite | ✅ 36 files |
| Final table suite | ✅ 29 files |
| `reports/FIGURE_INDEX.tsv` updated | ✅ 578 rows, 106 Phase 3 figures registered |
| Phase 2 objects unmodified | ✅ Phase 3 is read-only on Phase 2 |
| Prohibited analyses | ✅ none executed — no CNV, spatial, trajectory, velocity, survival, deep learning, or condition DE |
| No CCC method executed beyond inference | ✅ no Phase 4 spatial work begun |

## 2. Deliverables

- **`results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv`** — the primary deliverable,
  36,486 interactions traceable from methods → patients → receiver response → LochNESS →
  literature → priority tier.
- `results/phase3/phase3_manifest.json` — Phase 2 input and checksums, CCC input object and
  checksum, all five methods with versions and database versions, support rules, concordance
  classes and counts, prioritisation tiers, evidence-stream counts, sample recurrence,
  robustness, LochNESS results and implementation comparison, figures, tables, reports,
  scripts, environment with the recorded igraph change, seeds, every SLURM JobID including
  every failure and its root cause, git state, 11 scientific caveats, and the prohibitions
  respected.
- `results/phase3/figures/final/` (36) and `results/phase3/tables/final/` (29).
- Ten reports: method plan, two audits, LochNESS design, concordance, interactions,
  receiver response, LochNESS results, implementation comparison, robustness, handoff.

## 3. Headline scientific outcome

**547 tumour-centric interactions reach ≥3 independent evidence streams; 325 interactions
reach all 4.** The MPNST tumour compartment is predicted to act principally on **myeloid cells
and the vasculature**: a reproducible myeloid-directed signal set (APP→CD74 across five
receivers, CD99→PILRA, ANXA1→FPR1, HLA-F→LILRB1/2) whose published counterparts elsewhere are
immunosuppressive, with the myeloid *programme* independently best explained by **CSF1**; a
canonical angiogenic axis (VEGFA→KDR/FLT1/NRP1); and a **reciprocal perivascular Notch
circuit** in a tumour type where Notch is already implicated in Schwann-cell transformation.
The lymphoid picture is **mixed, not uniformly suppressive**. The fibroblast compartment shows
abundant LR co-expression but **no receiver-programme evidence**.

## 4. Negative results retained

LochNESS found **no** receiver-state structure associated with the APP context in any lineage
· NicheNet did **not** corroborate APP for macrophages · TGFB1→TGFBR was only weakly supported
· 5,485 interactions are discordant between frameworks · 53.6% of supported interactions rest
on a single patient · the fibroblast receiver-response result is a clear negative.

## 5. Phase 4 readiness

**READY WITH CAVEATS** for spatial validation — see `PHASE3_HANDOFF.md` §26. CNV inference
should precede or accompany it, to resolve the Phase 2 malignant-fraction uncertainty that
propagates into every tumour-centric result.

## 6. STOP

**Phase 3 is complete and frozen. Phase 4 was not begun and requires separate authorization.**
