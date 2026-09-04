#!/bin/bash
#SBATCH --job-name=p4_m28
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=logs/phase4/M28/m28_feasibility_%j.out
#SBATCH --error=logs/phase4/M28/m28_feasibility_%j.err
# Phase 4 · M28 — SCEVAN feasibility + per-sample raw-count extraction.
# Rerun after 19896576: SCEVAN 1.0.3 ignores output_dir in its plot/read-back helpers.
# Extraction outputs from 19896576 are reused (idempotent); only the smoke test re-runs.
# Memory rationale: the Phase 2 object is 6.0 GB on disk with an SCT scale.data
# layer; Phase 3 peaked at 27.4 GiB loading the same object. 96G leaves room for
# the load plus the smoke test's dense subsample. NOT the 450G maximum.
# NOTE: no "set -u". Conda's qt-main activate.d hook dereferences an unset
# QT_XCB_GL_INTEGRATION, so "set -u" aborts the script on conda activate/deactivate
# (root cause of SLURM 19896654). Conda's own guidance is to disable nounset.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase4/M28
echo "=== M28 start $(date) | Job ${SLURM_JOB_ID} | $(hostname) | ${SLURM_CPUS_PER_TASK} cpus ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
# SCEVAN's subclone stage calls umap::umap(method="umap-learn"), i.e. PYTHON
# umap-learn via reticulate. Point reticulate at the isolated p4_umap_env so
# R_env's frozen R stack is never touched (see install_umap_learn.sh).
export RETICULATE_PYTHON=/local/projects-t3/lilab/vmenon/anaconda3/envs/p4_umap_env/bin/python
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
Rscript scripts/phase4/scevan/m28_feasibility.R
rc=$?
echo "=== M28 end $(date) | Rscript exit ${rc} ==="
exit $rc
