#!/bin/bash
#SBATCH --job-name=p2_m12_evaluate
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=logs/phase2/slurm/m12_evaluate_%j.out
#SBATCH --error=logs/phase2/slurm/m12_evaluate_%j.err

# ---------------------------------------------------------------------------
# Phase 2 / M12 - Pre/post Harmony integration evaluation (production run).
#
# Compares the frozen Phase 1 'pca' 1:30 space with the M11 'postint_harmony'
# 1:30 space on identical cells, identical parameters and the same seed.
# Harmony is NOT re-run and NOT re-tuned. No clustering, markers or annotation.
#
# Resources sized from M11 accounting (JobID 19886486: 8.38 GiB MaxRSS) plus the
# two n x n distance matrices required for the silhouette widths (~1.55 GB each,
# computed sequentially), not from the 32 CPU / 450G / 72h envelope maximum.
# ---------------------------------------------------------------------------

set -euo pipefail

PROJECT_ROOT="/local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST"
cd "${PROJECT_ROOT}"

mkdir -p logs/phase2/slurm results/phase2/harmony/evaluation reports/phase2/figures/m12

echo "=== SLURM Job Execution Start ==="
echo "Date/Time  : $(date)"
echo "Host       : $(hostname)"
echo "Job ID     : ${SLURM_JOB_ID:-NO_JOB_ID}"
echo "Partition  : ${SLURM_JOB_PARTITION:-NA}"
echo "CPUs       : ${SLURM_CPUS_PER_TASK:-1}"
echo "Mem        : 96G"
echo "Walltime   : 04:00:00"
echo "Working Dir: $(pwd)"

source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

echo "--- Environment record ---"
echo "which R       : $(which R)"
R --version | head -1
echo "which Rscript : $(which Rscript)"
Rscript -e 'cat(sprintf("Seurat %s | SeuratObject %s | harmony %s | Matrix %s | cluster %s | RANN %s | uwot %s | ggplot2 %s\n", packageVersion("Seurat"), packageVersion("SeuratObject"), packageVersion("harmony"), packageVersion("Matrix"), packageVersion("cluster"), packageVersion("RANN"), packageVersion("uwot"), packageVersion("ggplot2")))'
echo "--------------------------"

Rscript scripts/R/phase2/evaluate_harmony.R \
  --input          results/phase2/harmony/phase2_harmony_integrated.rds \
  --out-dir        results/phase2/harmony/evaluation \
  --fig-dir        reports/phase2/figures/m12 \
  --group-by       sample_id \
  --dims           30 \
  --k-primary      15 \
  --k-robust       50 \
  --random-seed    42 \
  --expected-cells 19716

echo "=== SLURM Job Execution End ==="
echo "Finished at: $(date)"
