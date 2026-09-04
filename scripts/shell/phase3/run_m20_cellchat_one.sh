#!/bin/bash
#SBATCH --job-name=p3_cchat1
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/phase3/M20/cellchat_%x_%j.out
#SBATCH --error=logs/phase3/M20/cellchat_%x_%j.err
# Phase 3 / M20 - CellChat for ONE sample. Four of these run concurrently: the samples are
# scientifically independent after the shared M19 input prep, so per-sample parallelism cuts
# wall clock from ~5 h serial to ~1.5 h without changing any result.
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
SAMPLE="${1:?usage: sbatch run_m20_cellchat_one.sh <SAMPLE>}"
mkdir -p logs/phase3/M20 results/phase3/ccc/by_method
echo "=== M20 CellChat [$SAMPLE] start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) | CPUs ${SLURM_CPUS_PER_TASK:-1} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript scripts/phase3/ccc/run_cellchat.R \
  --input results/phase3/ccc/ccc_input_object.rds \
  --out-dir results/phase3/ccc/by_method \
  --samples "$SAMPLE" --tag "$SAMPLE" \
  --min-cells 10 --random-seed 42
echo "=== End [$SAMPLE]: $(date) ==="
