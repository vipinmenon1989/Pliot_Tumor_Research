#!/bin/bash
#SBATCH --job-name=p6_m45
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/phase6/m45_m50_%j.out
#SBATCH --error=logs/phase6/m45_m50_%j.err
# Phase 6 - M45 onward. M42, M43 and M44 already completed under 19899748 and
# their outputs are on disk, so they are not repeated.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase6 reports/phase6/environment
echo "=== Phase 6 (M45-M50) start $(date) | Job ${SLURM_JOB_ID} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
[ -f reports/phase6/environment/R_env_PRE_PHASE6.yml ] || \
  conda env export > reports/phase6/environment/R_env_PRE_PHASE6.yml
rc=0
run () { echo; echo "=== $1 ==="; Rscript "$@"; r=$?; echo "--- $1 exit=${r} ---"
         [ $r -ne 0 ] && { rc=$r; echo "STOPPING: $1 failed"; return 1; }; return 0; }
run scripts/phase6/regulatory/m45_tf_activity.R            || exit $rc
run scripts/phase6/pathways/m46_pathway_activity.R         || exit $rc
run scripts/phase6/cna_expression/m47_cna_expression.R     || exit $rc
run scripts/phase6/utils/m48_integrated_architecture.R     || exit $rc
run scripts/phase6/validation/m49_robustness.R             || exit $rc
run scripts/phase6/utils/m50_build_final_object.R          || exit $rc
run scripts/phase6/utils/m50_figures.R                     || exit $rc
run scripts/phase6/utils/m50_figures2.R                    || exit $rc
export P6_SLURM_JOBS="19899746,19899748,${SLURM_JOB_ID}"
run scripts/phase6/utils/m50_finalize.R                    || exit $rc
conda env export > reports/phase6/environment/R_env_POST_PHASE6.yml
if diff -q reports/phase6/environment/R_env_PRE_PHASE6.yml \
          reports/phase6/environment/R_env_POST_PHASE6.yml > /dev/null; then
  echo "R_env UNCHANGED across Phase 6 (PRE == POST)"
else
  echo "WARNING: R_env changed across Phase 6:"
  diff reports/phase6/environment/R_env_PRE_PHASE6.yml reports/phase6/environment/R_env_POST_PHASE6.yml
fi
echo "=== Phase 6 end $(date) | rc=${rc} ==="
exit $rc
