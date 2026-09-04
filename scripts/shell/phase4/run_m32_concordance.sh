#!/bin/bash
#SBATCH --job-name=p4_m32conc
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/phase4/M32/m32_concordance_%j.out
#SBATCH --error=logs/phase4/M32/m32_concordance_%j.err
# Phase 4 · M32 — concordance on the refined labels, then the Phase 3 vs Phase 4
# sensitivity comparison. build_concordance.R is reused UNMODIFIED, so the output
# schema is identical to Phase 3's and the comparison is column-for-column fair.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase4/M32 results/phase4/figures/M32
echo "=== M32-concordance start $(date) | Job ${SLURM_JOB_ID} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript scripts/phase3/concordance/build_concordance.R \
  --by-dir results/phase4/ccc_refinement/by_method \
  --out-dir results/phase4/ccc_refinement/concordance \
  --tab-dir results/phase4/tables \
  --fig-dir results/phase4/figures/M32
rc1=$?
echo; echo "--- Phase 3 vs Phase 4 sensitivity ---"
Rscript scripts/phase4/ccc_refinement/m32_ccc_sensitivity.R
rc2=$?
echo "=== M32-concordance end $(date) | concordance=${rc1} sensitivity=${rc2} ==="
[ $rc1 -ne 0 ] && exit $rc1
exit $rc2
