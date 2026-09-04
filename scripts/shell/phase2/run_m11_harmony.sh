#!/bin/bash
#SBATCH --job-name=p2_m11_harmony
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/phase2/slurm/m11_harmony_%j.out
#SBATCH --error=logs/phase2/slurm/m11_harmony_%j.err

# ---------------------------------------------------------------------------
# Phase 2 / M11 - Default Harmony integration (production run).
#
# Executes the configuration approved at the M10 gate:
#   grouping variable : sample_id
#   input reduction   : pca, dims 1:30 (frozen Phase 1 shared PCA, assay SCT)
#   output reduction  : postint_harmony
#   all Harmony scientific parameters at harmony 1.2.4 defaults
#
# The Phase 1 handoff object is opened read-only. The output is a new file.
# No Snakemake. No clustering, no markers, no annotation, no pre/post assessment.
#
# Resources are sized from M10 SLURM accounting (JobID 19886411: 4.73 GiB MaxRSS,
# 2m17s for checksum + load) plus the Phase 1 M8 benchmark (17.9 GB / 9m43s), not
# from the 32 CPU / 450G / 72h envelope maximum.
# ---------------------------------------------------------------------------

set -euo pipefail

PROJECT_ROOT="/local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST"
cd "${PROJECT_ROOT}"

mkdir -p logs/phase2/slurm results/phase2/harmony

echo "=== SLURM Job Execution Start ==="
echo "Date/Time  : $(date)"
echo "Host       : $(hostname)"
echo "Job ID     : ${SLURM_JOB_ID:-NO_JOB_ID}"
echo "Partition  : ${SLURM_JOB_PARTITION:-NA}"
echo "CPUs       : ${SLURM_CPUS_PER_TASK:-1}"
echo "Mem        : 64G"
echo "Walltime   : 04:00:00"
echo "Working Dir: $(pwd)"

source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

echo "--- Environment record ---"
echo "which R       : $(which R)"
R --version | head -1
echo "which Rscript : $(which Rscript)"
Rscript -e 'cat(sprintf("Seurat %s | SeuratObject %s | harmony %s | Matrix %s\n", packageVersion("Seurat"), packageVersion("SeuratObject"), packageVersion("harmony"), packageVersion("Matrix")))'
echo "--------------------------"

Rscript scripts/R/phase2/run_harmony_integration.R \
  --input          results/combined/pre_integration/combined_preintegration.rds \
  --out-rds        results/phase2/harmony/phase2_harmony_integrated.rds \
  --out-dir        results/phase2/harmony \
  --group-by       sample_id \
  --reduction-use  pca \
  --dims           30 \
  --reduction-save postint_harmony \
  --random-seed    42 \
  --expected-cells 19716

echo "=== SLURM Job Execution End ==="
echo "Finished at: $(date)"
