#!/bin/bash
#SBATCH --job-name=p4_m29
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=1-4
#SBATCH --output=logs/phase4/M29/m29_%A_%a.out
#SBATCH --error=logs/phase4/M29/m29_%A_%a.err
# Phase 4 · M29 — per-patient SCEVAN. One array task per patient (§14/§19).
#
# Resource rationale (M28 measured, not guessed):
#   largest dense annotated matrix  MPNST_1  30,121 x 7,615 doubles = 1.71 GiB
#   SCEVAN holds smoothed + relative + segmented copies and a parallelDist
#   distance matrix (7,615^2 x 8 B = 0.46 GiB), so 64G is several times the
#   working set. This is NOT the 450G maximum.
#   8 CPUs matches par_cores exactly - no oversubscription (§9).
# NOTE: no "set -u". Conda's qt-main activate.d hook dereferences an unset
# QT_XCB_GL_INTEGRATION, so "set -u" aborts the script on conda activate/deactivate
# (root cause of SLURM 19896654). Conda's own guidance is to disable nounset.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase4/M29
SAMPLES=(MPNST_1 MPNST_2 MPNST_3 MPNST_4)
S=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}
echo "=== M29 start $(date) | Job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | ${S} | $(hostname) | ${SLURM_CPUS_PER_TASK} cpus ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
# SCEVAN's subclone stage calls umap::umap(method="umap-learn"), i.e. PYTHON
# umap-learn via reticulate. Point reticulate at the isolated p4_umap_env so
# R_env's frozen R stack is never touched (see install_umap_learn.sh).
export RETICULATE_PYTHON=/local/projects-t3/lilab/vmenon/anaconda3/envs/p4_umap_env/bin/python
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
Rscript scripts/phase4/scevan/m29_run_scevan.R --sample "${S}" --cores "${SLURM_CPUS_PER_TASK}"
rc=$?
echo "=== M29 end $(date) | ${S} | exit ${rc} ==="
exit $rc
