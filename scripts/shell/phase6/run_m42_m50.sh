#!/bin/bash
#SBATCH --job-name=p6_all
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/phase6/m42_m50_%j.out
#SBATCH --error=logs/phase6/m42_m50_%j.err
# =============================================================================
# Phase 6 - M42 through M50.
#
# 128 G: M42 loads the Phase 5 object (the Phase 4 object peaked at 29 GiB in
# M36 and Phase 5 adds metadata only), and M45/M46 hold a dense
# ~19,000 x 6,434 expression matrix plus its z-scored copy for decoupleR
# (~2 GB each). 8 CPUs for decoupleR. Evidence-based, well under the 450 G max.
# =============================================================================
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase6 reports/phase6/environment
echo "=== Phase 6 start $(date) | Job ${SLURM_JOB_ID} | node $(hostname) ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
conda env export > reports/phase6/environment/R_env_PRE_PHASE6.yml

rc=0
run () {
  echo; echo "=== $1 ==="
  Rscript "$@"; r=$?
  echo "--- $1 exit=${r} ---"
  [ $r -ne 0 ] && { rc=$r; echo "STOPPING: $1 failed"; return 1; }
  return 0
}

run scripts/phase6/utils/m42_feasibility.R                      || exit $rc
run scripts/phase6/clone_program/m43_clone_program.R            || exit $rc
run scripts/phase6/plasticity/m44_within_clone_diversity.R      || exit $rc
run scripts/phase6/regulatory/m45_tf_activity.R                 || exit $rc
run scripts/phase6/pathways/m46_pathway_activity.R              || exit $rc
run scripts/phase6/cna_expression/m47_cna_expression.R          || exit $rc
run scripts/phase6/utils/m48_integrated_architecture.R          || exit $rc
run scripts/phase6/validation/m49_robustness.R                  || exit $rc
run scripts/phase6/utils/m50_build_final_object.R               || exit $rc
run scripts/phase6/utils/m50_figures.R                          || exit $rc
run scripts/phase6/utils/m50_figures2.R                         || exit $rc
export P6_SLURM_JOBS="${SLURM_JOB_ID}"
run scripts/phase6/utils/m50_finalize.R                         || exit $rc

conda env export > reports/phase6/environment/R_env_POST_PHASE6.yml
if diff -q reports/phase6/environment/R_env_PRE_PHASE6.yml \
          reports/phase6/environment/R_env_POST_PHASE6.yml > /dev/null; then
  echo "R_env UNCHANGED across Phase 6 (PRE == POST)"
else
  echo "WARNING: R_env changed across Phase 6:"
  diff reports/phase6/environment/R_env_PRE_PHASE6.yml \
       reports/phase6/environment/R_env_POST_PHASE6.yml
fi
echo "=== Phase 6 end $(date) | rc=${rc} ==="
exit $rc
