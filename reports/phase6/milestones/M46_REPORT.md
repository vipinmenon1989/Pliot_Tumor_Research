# M46 — Pathway Activity

**Date** 2026-09-04 · **Status** COMPLETE · **SLURM** 19899749

**Method.** Two independent layers, reported side by side and **never merged into a composite score
or a single numeric rank**: **PROGENy** (14 footprint pathways, top-500 weights, decoupleR
`run_mlm`) and **MSigDB v2024.1.Hs Hallmark** (50 sets, mean z-score of the set's genes). Input is
the frozen RNA log-normalised layer.

**Result — the pathway layer independently supports the gene-derived program labels.**

| program | PROGENy (pooled ρ; `*` = concordant in ≥ 3 patients) | Hallmark |
| --- | --- | --- |
| P1 Translation_ribosomal | PI3K −0.71*, p53 −0.69, Hypoxia −0.63*, MAPK −0.40* | OXPHOS 0.13 |
| P2 Neuronal | PI3K 0.51*, MAPK 0.51*, JAK-STAT −0.52, TNFa −0.39* | UV_RESPONSE_DN 0.56*, HEDGEHOG 0.52*, WNT_BETA_CATENIN 0.47* |
| P3 **Mesenchymal_ECM** | **EGFR 0.57***, MAPK −0.53*, Estrogen 0.42*, TNFa 0.36* | COAGULATION 0.43*, OXPHOS 0.42*, ROS 0.41*, **EMT 0.33*** |
| P4 **Hypoxia_Angio** | **Hypoxia 0.47***, PI3K 0.45*, EGFR −0.39* | UV_RESPONSE_DN 0.53*, HEME_METABOLISM 0.44* |
| P5 Translation_ribosomal | p53 0.53*, Hypoxia 0.44*, PI3K 0.37* | **MTORC1 0.58***, **GLYCOLYSIS 0.55***, PI3K_AKT_MTOR 0.50* |
| P6 Schwann_like | EGFR 0.32*, JAK-STAT 0.32*, TNFa 0.28*, VEGF 0.20* | OXPHOS 0.35*, MYC_TARGETS_V1 0.27*, IFN_ALPHA 0.27* |
| P7 **Cycling** | Trail −0.19, NFkB −0.18 | **E2F_TARGETS 0.43***, **MYC_TARGETS_V1 0.44***, DNA_REPAIR 0.33* |
| P8 Myeloid_ambient | JAK-STAT 0.25* | — |

Three of these are worth naming because they were **not** imposed and they cross-check labels
derived from an entirely separate source (the programs' own top genes):

* **P4 was labelled `Hypoxia_Angio` from its genes; the PROGENy Hypoxia footprint is its top
  pathway (0.47*).**
* **P7 was labelled `Cycling` from its genes; Hallmark E2F_TARGETS and MYC_TARGETS_V1 are its top
  sets, matching the E2F4/E2F1/MYC regulators found in M45.**
* **P3 was labelled `Mesenchymal_ECM` from its genes; Hallmark EMT is among its top sets.**

**P5's mTORC1 / glycolysis / PI3K-AKT signature is the informative negative**: it confirms that this
ribosomal-pseudogene factor is a coherent *translation-and-growth* axis rather than pure noise, but
that does not make it a distinct malignant *identity*, and it is still flagged technical-dominated.

**Integration.** `PROGRAM_REGULATORY_EVIDENCE.tsv` lists only cross-patient-concordant rows from all
three layers side by side. **No weighted composite is computed.**

**Outputs.** `PROGRAM_PATHWAY_ACTIVITY.tsv` · `PROGRAM_REGULATORY_EVIDENCE.tsv` ·
`m46_progeny_activity_matrix.rds` (6,434 × 14) · `m46_hallmark_score_matrix.rds` (6,434 × 50)

**Next.** M47 — broad CNA → transcriptional consequences.
