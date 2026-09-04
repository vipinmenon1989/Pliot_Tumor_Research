# M47 — Broad CNA → Transcriptional Consequences

**Date** 2026-09-04 · **Status** COMPLETE · **SLURM** 19899749

> **THE LIMITATION THAT DEFINES THIS MILESTONE.** SCEVAN infers copy number **from expression**.
> An association between an inferred broad event and the expression of genes on the affected
> chromosome is therefore an **internal transcriptional consistency analysis — NOT independent
> validation** of the copy-number call. Agreement is expected, not confirmatory.

**Design.** Within each reliable patient, cells in clones carrying an inferred broad segment are
compared with cells in clones of the **same patient** that do not carry it, using the mean expression
of **all genes in that interval**. Segment level only; no single-gene claim anywhere. Declared
minima: segment ≥ 10 Mb, ≥ 30 genes in the interval, ≥ 20 cells per clone, and both a carrier and a
non-carrier clone in the same patient.

**Inputs.** 56 broad clonal non-neutral segments in MPNST_1/2/4 and 308 broad subclonal segments
across 19 clone labels, all from the frozen Phase 4 `SCEVAN_CNV_SUMMARY.tsv`. Gene coordinates from
the stored SCEVAN annotation (10,330 genes on chr1–22; 9,910 after the 5% detection filter).
MPNST_3 excluded.

**Result.** 299 clone-level events tested, **263 evaluable**. **Direction matches the inferred event
in 180 of 263 (68%)**. Largest effects are all in MPNST_1's chr8 gains: Cohen's d 3.23, 3.09, 2.97,
2.65 across four adjacent chr8 segments, with chr19 and chr11 losses at d ≈ −1.6 and −1.3.

**68% is the honest number.** It is well above chance and confirms that the clone assignments carry
real, coordinated, chromosome-scale expression structure — but a third of the evaluable events do
*not* move in the inferred direction, which is what one should expect when the events are themselves
expression-derived and when subclonal segment boundaries are approximate.

**Loci.** Segments containing NF1 (chr17) and NF2 (chr22) are tabulated with **both** the permitted
and the prohibited wording printed next to each row, so the distinction cannot be lost downstream:
permitted — *"a broad inferred loss segment of N Mb on chr17 whose interval contains the NF1
locus"*; prohibited — *"NF1 deletion"* or *"NF1-deleted cells"*.

**CNA ↔ program.** Clone-level carrier/non-carrier comparisons were also run against each Phase 5
program within patient, in `CNA_PROGRAM_ASSOCIATION.tsv`. CNA architecture is patient-specific and
is never pooled across patients.

**Outputs.** `BROAD_CNA_EXPRESSION_EFFECTS.tsv` (299 rows) · `CNA_PROGRAM_ASSOCIATION.tsv` ·
`M47_BROAD_SEGMENTS_CONTAINING_LOCI.tsv`

**Next.** M48 — integrated architecture.
