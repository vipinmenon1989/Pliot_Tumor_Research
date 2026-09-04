#!/bin/bash
#SBATCH --job-name=p6_figs
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/phase6/m50figs_%j.out
#SBATCH --error=logs/phase6/m50figs_%j.err
# Phase 6 - re-render figures after visual QC and re-finalize so the manifest
# checksums match the figures actually on disk.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
export P6_SLURM_JOBS="19899746,19899748,19899749,${SLURM_JOB_ID}"
echo "=== start $(date) | Job ${SLURM_JOB_ID} ==="
Rscript scripts/phase6/utils/m50_figures2.R || exit 1
Rscript scripts/phase6/utils/m50_finalize.R || exit 1
echo "=== end $(date) ==="
