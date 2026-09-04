#!/bin/bash
#SBATCH --job-name=p4_m32sens
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=logs/phase4/M32/m32_sensitivity_%j.out
#SBATCH --error=logs/phase4/M32/m32_sensitivity_%j.err
# Rerun ONLY the Phase 3 vs Phase 4 sensitivity comparison. The concordance step
# of 19897179 exited 0 and its outputs are intact; only this script failed.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript scripts/phase4/ccc_refinement/m32_ccc_sensitivity.R
rc=$?; echo "=== exit ${rc} ==="; exit $rc
