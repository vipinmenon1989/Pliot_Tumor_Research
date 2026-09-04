#!/bin/bash
#SBATCH --job-name=p5_m41c
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/phase5/M41/m41c_%j.out
#SBATCH --error=logs/phase5/M41/m41c_%j.err
# Phase 5 - figures 08/09/12 and the Phase 5 freeze. M40 and figures 01-07,
# 10-11 already succeeded under 19899742 and are not repeated.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
K=$(cat results/phase5/validation/selected_K.txt)
echo "=== start $(date) | Job ${SLURM_JOB_ID} | K=${K} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
export P5_SLURM_JOBS="19899333,19899335,19899339,19899353,19899354,19899355,19899356,19899358,19899359,19899742,${SLURM_JOB_ID}"
rc=0
Rscript scripts/phase5/utils/m41_figures2.R "${K}"; r=$?; echo "figures2 exit=$r"; [ $r -ne 0 ] && exit $r
Rscript scripts/phase5/utils/m41_finalize.R "${K}"; r=$?; echo "finalize exit=$r"; [ $r -ne 0 ] && exit $r
conda env export > reports/phase5/environment/R_env_POST_PHASE5.yml
if diff -q reports/phase5/environment/R_env_PRE_PHASE5.yml reports/phase5/environment/R_env_POST_PHASE5.yml > /dev/null; then
  echo "R_env UNCHANGED across Phase 5 (PRE == POST)"
else
  echo "WARNING: R_env changed:"; diff reports/phase5/environment/R_env_PRE_PHASE5.yml reports/phase5/environment/R_env_POST_PHASE5.yml
fi
echo "=== end $(date) ==="
