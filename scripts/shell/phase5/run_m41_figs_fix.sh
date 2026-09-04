#!/bin/bash
#SBATCH --job-name=p5_m41d
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/phase5/M41/m41d_%j.out
#SBATCH --error=logs/phase5/M41/m41d_%j.err
# Phase 5 - re-render the figures after visual QC and re-finalize so the
# manifest checksums match the figures actually on disk.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
K=$(cat results/phase5/validation/selected_K.txt)
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
export P5_SLURM_JOBS="19899333,19899335,19899339,19899353,19899354,19899355,19899356,19899358,19899359,19899742,19899744,${SLURM_JOB_ID}"
echo "=== start $(date) | Job ${SLURM_JOB_ID} | K=${K} ==="
Rscript scripts/phase5/utils/m41_figures.R  "${K}" || exit 1
Rscript scripts/phase5/utils/m41_figures2.R "${K}" || exit 1
Rscript scripts/phase5/utils/m41_finalize.R "${K}" || exit 1
echo "=== end $(date) ==="
