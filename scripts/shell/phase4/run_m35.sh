#!/bin/bash
#SBATCH --job-name=p4_m35
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=logs/phase4/M35/m35_%j.out
#SBATCH --error=logs/phase4/M35/m35_%j.err
# Phase 4 · M35 — build + validate the final object, then assemble and freeze.
# 96 G because it holds the Phase 2 object AND a second copy for the column-by-
# column preservation comparison, then reloads the saved object to validate it.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase4/M35
echo "=== M35 start $(date) | Job ${SLURM_JOB_ID} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
echo; echo "=== build + validate the final Phase 4 object ==="
Rscript scripts/phase4/utils/m35_build_final_object.R; rc1=$?
echo; echo "=== CNA heatmap (§20, §52) ==="
Rscript scripts/phase4/utils/m35_cna_heatmap.R;        rc0=$?
echo; echo "=== assemble figures, tables, manifest, figure index ==="
Rscript scripts/phase4/utils/m35_assemble.R;           rc2=$?
echo "=== M35 end $(date) | build=${rc1} cna_heatmap=${rc0} assemble=${rc2} ==="
[ $rc1 -ne 0 ] && exit $rc1
exit $rc2
