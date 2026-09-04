#!/bin/bash
#SBATCH --job-name=p5_m36
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --output=logs/phase5/M36/m36_%j.out
#SBATCH --error=logs/phase5/M36/m36_%j.err
# =============================================================================
# Phase 5 - M36 - reconstruction and feasibility.
#
# 96 G: the Phase 4 object is 5.64 GB on disk and prior full-object jobs peaked
# at 27.78 GiB (M32 refined input); this job additionally joins the four
# per-sample RNA layers and slices two 31,764-gene matrices out of them, so the
# request is set above that observed peak rather than at the 450 G maximum.
# 4 CPUs: the work is single-threaded I/O and sparse-matrix slicing.
# =============================================================================
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase5/M36
echo "=== M36 start $(date) | Job ${SLURM_JOB_ID} | node $(hostname) ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript scripts/phase5/utils/m36_feasibility.R; rc=$?
echo "=== M36 end $(date) | rc=${rc} ==="
exit $rc
