#!/bin/bash
#SBATCH --job-name=p4_m35a
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/phase4/M35A/m35a_%j.out
#SBATCH --error=logs/phase4/M35A/m35a_%j.err
# =============================================================================
# Phase 4 - M35A - SCEVAN figure consolidation and evidence visualization.
#
# Resources are deliberately modest and evidence-based. Nothing here re-runs an
# analysis: the inputs are the 5 MB frozen malignancy metadata table, the small
# frozen TSVs, and the stored native SCEVAN PNG/RData outputs. The largest
# single object touched is MPNST_1's stored CNA matrix (9 MB on disk, ~0.5 GB
# in memory), so 24 G is ample. The 6 GB phase4_final_object.rds is NEVER
# loaded - it is not needed and must not be modified.
# 2 CPUs because every step is single-threaded plotting.
# =============================================================================
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase4/M35A results/phase4/figures/final
echo "=== M35A start $(date) | Job ${SLURM_JOB_ID} | node $(hostname) ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
export M35A_SLURM_JOB_ID="${SLURM_JOB_ID}"

rc=0
run () {
  echo; echo "=== $1 ==="
  Rscript "$1"; r=$?
  echo "--- $1 exit=${r} ---"
  [ $r -ne 0 ] && rc=$r
  return 0
}

run scripts/phase4/figures/m35a_native_figures.R
run scripts/phase4/figures/m35a_clone_composition.R
run scripts/phase4/figures/m35a_fibroblast_cna.R
run scripts/phase4/figures/m35a_transition.R
run scripts/phase4/figures/m35a_mpnst3_qc.R
run scripts/phase4/figures/m35a_threshold.R
run scripts/phase4/figures/m35a_broad_cna.R
run scripts/phase4/figures/m35a_summary.R

echo; echo "=== M35A end $(date) | overall rc=${rc} ==="
exit $rc
