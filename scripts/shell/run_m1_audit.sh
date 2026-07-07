#!/bin/bash
#SBATCH --job-name=m1_real_audit
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/m1_audit_%j.out
#SBATCH --error=logs/slurm/m1_audit_%j.err

# Safe script configurations
set -euo pipefail

# Ensure SLURM log directory exists
mkdir -p logs/slurm

echo "=== SLURM Job Execution Start ==="
echo "Date/Time : $(date)"
echo "Host      : $(hostname)"
echo "Job ID    : ${SLURM_JOB_ID:-NO_JOB_ID}"
echo "CPUs      : ${SLURM_CPUS_PER_TASK:-1}"
echo "Mem       : 128G"
echo "Working Dir: $(pwd)"

# Source conda and activate R_env
echo "Activating environment R_env..."
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

# Force rerun of audit to ensure clean execution on the real dataset
# This is clean because Snakemake will detect input_rds has changed compared to the synthetic data local run.
echo "Executing Snakemake audit..."
snakemake --cores 8 --rerun-incomplete --keep-going

echo "=== SLURM Job Execution End ==="
echo "Finished at: $(date)"
