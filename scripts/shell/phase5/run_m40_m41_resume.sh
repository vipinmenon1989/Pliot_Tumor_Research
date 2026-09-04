#!/bin/bash
#SBATCH --job-name=p5_m41b
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/phase5/M41/m41b_%j.out
#SBATCH --error=logs/phase5/M41/m41b_%j.err
# =============================================================================
# Phase 5 - M40 robustness (re-run with the scope-aware verdict) + figures +
# finalize.
#
# The object build is NOT repeated: 19899359 completed it and its 8 reload
# validations passed. Rebuilding a 5.6 GB object that is already verified would
# burn 10 minutes to produce the same bytes. 48 G rather than 128 G because
# nothing here loads the Phase 4 object.
# =============================================================================
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase5/M41
K=$(cat results/phase5/validation/selected_K.txt)
echo "=== M41 resume start $(date) | Job ${SLURM_JOB_ID} | K=${K} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
export P5_SLURM_JOBS="19899333,19899335,19899339,19899353,19899354,19899355,19899356,19899358,19899359,${SLURM_JOB_ID}"
rc=0
run () { echo; echo "=== $* ==="; Rscript "$@"; r=$?; echo "--- exit=${r} ---"
         [ $r -ne 0 ] && { rc=$r; return 1; }; return 0; }
run scripts/phase5/validation/m40_robustness.R "${K}" 0_1 || exit $rc
run scripts/phase5/utils/m41_figures.R         "${K}"     || exit $rc
run scripts/phase5/utils/m41_figures2.R        "${K}"     || exit $rc
run scripts/phase5/utils/m41_finalize.R        "${K}"     || exit $rc
conda env export > reports/phase5/environment/R_env_POST_PHASE5.yml
if diff -q reports/phase5/environment/R_env_PRE_PHASE5.yml \
          reports/phase5/environment/R_env_POST_PHASE5.yml > /dev/null; then
  echo "R_env UNCHANGED across Phase 5 (PRE == POST)"
else
  echo "WARNING: R_env changed across Phase 5:"
  diff reports/phase5/environment/R_env_PRE_PHASE5.yml \
       reports/phase5/environment/R_env_POST_PHASE5.yml
fi
echo "=== M41 resume end $(date) | rc=${rc} ==="
exit $rc
