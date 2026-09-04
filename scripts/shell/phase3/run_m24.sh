#!/bin/bash
#SBATCH --job-name=p3_m24_loch
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=logs/phase3/M24/m24_%j.out
#SBATCH --error=logs/phase3/M24/m24_%j.err
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase3/M24 results/phase3/figures/M24
echo "=== M24 LochNESS start: $(date) | Job ${SLURM_JOB_ID:-NA} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
echo "--- LochNESS (MMCA formulation) ---"
Rscript scripts/phase3/lochness/lochness_mpnst.R \
  --input results/phase3/ccc/ccc_input_object.rds \
  --context results/phase3/lochness/lochness_context_labels.tsv \
  --out-dir results/phase3/lochness --fig-dir results/phase3/figures/M24 \
  --lineages Macrophage,Fibroblast,Endothelial,CD8-T --n-perm 100 --random-seed 42
echo "--- R vs Python implementation comparison ---"
Rscript scripts/phase3/lochness/compare_implementations.R \
  --object results/phase3/ccc/ccc_input_object.rds \
  --context results/phase3/lochness/lochness_context_labels.tsv \
  --out-dir results/phase3/lochness --fig-dir results/phase3/figures/M24 \
  --lineage Macrophage --n-sub 2000
echo "=== End: $(date) ==="
