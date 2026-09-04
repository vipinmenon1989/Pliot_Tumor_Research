#!/bin/bash
#SBATCH --job-name=p4_m32prep
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --output=logs/phase4/M32/m32_prepare_%j.out
#SBATCH --error=logs/phase4/M32/m32_prepare_%j.err
# Phase 4 · M32 step 1 — build the refined CCC input object.
# 96 G because this loads the 6 GB Phase 2 object (Phase 3's equivalent step
# peaked near 30 GiB) and joins + renormalises the RNA assay.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase4/M32
echo "=== M32-prepare start $(date) | Job ${SLURM_JOB_ID} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript scripts/phase4/ccc_refinement/prepare_refined_ccc_inputs.R \
  --input results/phase2/phase2_final_object.rds \
  --calls results/phase4/malignancy/PHASE4_MALIGNANCY_CALLS.tsv \
  --out-dir results/phase4/ccc_refinement \
  --cpdb-dir results/phase4/ccc_refinement/cellphonedb_inputs \
  --tab-dir results/phase4/tables \
  --min-cells 10
rc=$?; echo "=== M32-prepare end $(date) | exit ${rc} ==="; exit $rc
