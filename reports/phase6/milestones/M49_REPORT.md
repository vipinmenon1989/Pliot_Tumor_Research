# M49 — Phase 6 Robustness

**Date** 2026-09-04 · **Status** COMPLETE · **SLURM** 19899749

Five perturbations: clone-size threshold (20 / 50 / 100 cells), dominant-program margin
(0 / 0.05 / 0.10 / 0.15 / 0.20), high-confidence-only restriction, TF and pathway recomputation under
that restriction, and patient-level direction consistency. Leave-one-patient-out is deliberately not
applied to the clone analyses, which are inherently within-patient.

**Headline.** Median η² is **0.0586 at every clone-size threshold**; median effective programs per
clone moves only 1.222 → 1.185 across a fourfold change in the margin; and the High-confidence
restriction reproduces η² at Spearman **0.901**, TF associations at **0.940** and PROGENy at
**0.962**. **Every Phase 6 conclusion survives.**

**Reported rather than smoothed.** 77–81% of strong pooled TF/pathway associations are concordant in
≥ 3 patients — a number that looks better than it is, because with n = 4 three-of-four sign agreement
occurs ≈ 31% of the time by chance. The report says so.

**Outputs.** `PHASE6_ROBUSTNESS.tsv` (14 rows) · `reports/phase6/ROBUSTNESS_REPORT.md`

**Next.** M50 — freeze.
