#!/bin/bash
#SBATCH --job-name=p4_m32cc
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --array=1-4
#SBATCH --output=logs/phase4/M32/m32_cellchat_%A_%a.out
#SBATCH --error=logs/phase4/M32/m32_cellchat_%A_%a.err
# Phase 4 · M32 — CellChat on the refined labels, one array task per patient.
# Phase 3 learned this the hard way: the serial CellChat run (19895564) was
# cancelled after projecting ~5 h, and four concurrent per-sample jobs finished
# in 44 min with identical results. Same structure here.
# scripts/phase3/ccc/run_cellchat.R is reused UNMODIFIED.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase4/M32
SAMPLES=(MPNST_1 MPNST_2 MPNST_3 MPNST_4)
S=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}
echo "=== M32-cellchat start $(date) | ${S} | Job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript -e 'cat(sprintf("CellChat %s\n", packageVersion("CellChat")))'
Rscript scripts/phase3/ccc/run_cellchat.R \
  --input results/phase4/ccc_refinement/ccc_input_object_refined.rds \
  --out-dir results/phase4/ccc_refinement/by_method \
  --samples "${S}" --tag "${S}" \
  --min-cells 10 --random-seed 42
rc=$?; echo "=== M32-cellchat end $(date) | ${S} | exit ${rc} ==="; exit $rc
