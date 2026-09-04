#!/bin/bash
#SBATCH --job-name=p4_m30
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/phase4/M30/m30_%j.out
#SBATCH --error=logs/phase4/M30/m30_%j.err
# Phase 4 · M30 — malignancy integration + figures.
# Memory: the CNA matrices are the largest objects (~30k genes x 7.6k cells dense
# for MPNST_1, loaded one sample at a time and freed), plus a 19,716-cell Seurat
# object built from the M28 counts for module scoring. 64 G, not the maximum.
# NOTE: no "set -u" — conda's qt-main activate.d hook dereferences an unset
# QT_XCB_GL_INTEGRATION (root cause of 19896654).
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase4/M30
echo "=== M30 start $(date) | Job ${SLURM_JOB_ID} | $(hostname) | ${SLURM_CPUS_PER_TASK} cpus ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
export RETICULATE_PYTHON=/local/projects-t3/lilab/vmenon/anaconda3/envs/p4_umap_env/bin/python
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo; echo "--- M29 aggregation ---"
Rscript scripts/phase4/scevan/m29_aggregate.R;            rc1=$?
echo; echo "--- M30 malignancy integration ---"
Rscript scripts/phase4/malignancy/m30_integrate_malignancy.R; rc2=$?
echo; echo "--- M30 figures ---"
Rscript scripts/phase4/malignancy/m30_figures.R;          rc3=$?
echo "=== M30 end $(date) | aggregate=${rc1} integrate=${rc2} figures=${rc3} ==="
[ $rc1 -ne 0 ] && exit $rc1
[ $rc2 -ne 0 ] && exit $rc2
exit $rc3
