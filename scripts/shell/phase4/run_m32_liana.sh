#!/bin/bash
#SBATCH --job-name=p4_m32liana
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/phase4/M32/m32_liana_%j.out
#SBATCH --error=logs/phase4/M32/m32_liana_%j.err
# Phase 4 · M32 — LIANA on the refined labels.
# scripts/phase3/ccc/run_liana.R is reused UNMODIFIED, with the identical
# resource (Consensus), min-cells (10) and seed (42) Phase 3 used. Only the
# input object's ccc_label differs. That is what makes this a label-sensitivity
# analysis rather than a comparison confounded by tool drift.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase4/M32
echo "=== M32-liana start $(date) | Job ${SLURM_JOB_ID} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
echo "--- versions must match Phase 3 exactly ---"
Rscript -e 'cat(sprintf("liana %s | OmnipathR %s | decoupleR %s | Seurat %s | SeuratObject %s\n", packageVersion("liana"), packageVersion("OmnipathR"), packageVersion("decoupleR"), packageVersion("Seurat"), packageVersion("SeuratObject")))'
Rscript scripts/phase3/ccc/run_liana.R \
  --input results/phase4/ccc_refinement/ccc_input_object_refined.rds \
  --out-dir results/phase4/ccc_refinement/by_method \
  --min-cells 10 --resource Consensus --random-seed 42
rc=$?; echo "=== M32-liana end $(date) | exit ${rc} ==="; exit $rc
