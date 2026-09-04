#!/bin/bash
#SBATCH --job-name=p2_m10_inspect
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=logs/phase2/slurm/m10_inspect_%j.out
#SBATCH --error=logs/phase2/slurm/m10_inspect_%j.err

# Phase 2 / M10: read-only structural validation of the Phase 1 handoff object.
# No Harmony, no clustering, no writes to any Phase 1 artifact.

set -euo pipefail

PROJECT_ROOT="/local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST"
cd "${PROJECT_ROOT}"

mkdir -p logs/phase2/slurm results/phase2/handoff

echo "=== SLURM Job Execution Start ==="
echo "Date/Time  : $(date)"
echo "Host       : $(hostname)"
echo "Job ID     : ${SLURM_JOB_ID:-NO_JOB_ID}"
echo "CPUs       : ${SLURM_CPUS_PER_TASK:-1}"
echo "Mem        : 96G"
echo "Walltime   : 02:00:00"
echo "Working Dir: $(pwd)"

source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

echo "R binary   : $(which R)"
echo "Rscript    : $(which Rscript)"

Rscript scripts/R/phase2/inspect_phase1_handoff.R \
  --input results/combined/pre_integration/combined_preintegration.rds \
  --out-dir results/phase2/handoff

echo "=== SLURM Job Execution End ==="
echo "Finished at: $(date)"
