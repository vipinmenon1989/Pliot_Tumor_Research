#!/bin/bash
#SBATCH --job-name=p5_m40m41
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=logs/phase5/M41/m40m41_%j.out
#SBATCH --error=logs/phase5/M41/m40m41_%j.err
# =============================================================================
# Phase 5 - M40 robustness + M41 object build, figures and freeze.
#
# 128 G: M41 loads the 5.64 GB Phase 4 object (M36 peaked at 29 GiB), keeps a
# metadata copy for the column-by-column preservation comparison, projects all
# 19,716 cells onto the fixed spectra, then reloads the saved object to
# validate it. Evidence-based; the 450 G maximum is not requested.
# =============================================================================
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase5/M41 reports/phase5/environment
K=$(cat results/phase5/validation/selected_K.txt)
echo "=== M40/M41 start $(date) | Job ${SLURM_JOB_ID} | K=${K} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
rc=0
run () { echo; echo "=== $* ==="; Rscript "$@"; r=$?; echo "--- exit=${r} ---"
         [ $r -ne 0 ] && { rc=$r; return 1; }; return 0; }

run scripts/phase5/validation/m40_robustness.R      "${K}" 0_1 || exit $rc
run scripts/phase5/utils/m41_build_final_object.R   "${K}" 0_1 || exit $rc
run scripts/phase5/utils/m41_figures.R              "${K}"     || exit $rc
run scripts/phase5/utils/m41_figures2.R             "${K}"     || exit $rc
export P5_SLURM_JOBS="19899333,19899335,19899339,19899353,19899354,19899355,19899356,19899358,${SLURM_JOB_ID}"
run scripts/phase5/utils/m41_finalize.R             "${K}"     || exit $rc

conda env export > reports/phase5/environment/R_env_POST_PHASE5.yml
if diff -q reports/phase5/environment/R_env_PRE_PHASE5.yml \
          reports/phase5/environment/R_env_POST_PHASE5.yml > /dev/null; then
  echo "R_env UNCHANGED across Phase 5 (PRE == POST)"
else
  echo "WARNING: R_env changed across Phase 5:"
  diff reports/phase5/environment/R_env_PRE_PHASE5.yml \
       reports/phase5/environment/R_env_POST_PHASE5.yml
fi
echo "=== M40/M41 end $(date) | rc=${rc} ==="
exit $rc
