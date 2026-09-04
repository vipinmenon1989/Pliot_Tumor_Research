#!/bin/bash
#SBATCH --job-name=p4_m35a_fin
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/phase4/M35A/m35a_finalize_%j.out
#SBATCH --error=logs/phase4/M35A/m35a_finalize_%j.err
# Phase 4 - M35A finalize: checksums, FIGURE_INDEX.tsv rows, manifest block.
# 8 G and 1 CPU: this reads small TSV/JSON files and md5s ~18 figure files.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase4/M35A
echo "=== M35A finalize start $(date) | Job ${SLURM_JOB_ID} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
export M35A_FINAL_JOBID="${M35A_FINAL_JOBID:-19899312}"
Rscript scripts/phase4/figures/m35a_finalize.R; rc=$?
echo "=== M35A finalize end $(date) | rc=${rc} ==="
exit $rc
