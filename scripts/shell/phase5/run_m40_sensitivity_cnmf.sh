#!/bin/bash
#SBATCH --job-name=p5_m40cnmf
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/phase5/M40/m40cnmf_%j.out
#SBATCH --error=logs/phase5/M40/m40cnmf_%j.err
# =============================================================================
# Phase 5 - M40b - leave-one-patient-out and high-confidence cNMF runs.
#
# Only the selected K is factorized (not the whole grid) and n_iter is 50
# rather than 100, because these runs exist to test whether the PRIMARY
# programs are recoverable, not to select a rank a second time. That is an
# evidence-based reduction in cost, not a shortcut past a needed analysis.
# K is passed in as $1.
# =============================================================================
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase5/M40
K=${1:?"usage: sbatch run_m40_sensitivity_cnmf.sh <K>"}
echo "=== M40b start $(date) | Job ${SLURM_JOB_ID} | K=${K} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh

echo; echo "=== M40a: prepare matrices (R_env) ==="
conda activate R_env
Rscript scripts/phase5/validation/m40_prepare_sensitivity_inputs.R || exit 1
conda deactivate

conda activate p5_cnmf_env
OUTDIR=results/phase5/programs/cnmf
IN=results/phase5/programs/cnmf_input
NITER=50; SEED=42; NUMGENES=2000; WORKERS=8
rc=0
for RUN in loo_MPNST_1 loo_MPNST_2 loo_MPNST_3 loo_MPNST_4 highconf; do
  echo; echo "--- ${RUN} $(date) ---"
  python scripts/phase5/programs/m37_build_anndata.py \
      "${IN}/${RUN}" "${IN}/${RUN}/counts.h5ad" \
      "${IN}/malignant_cell_covariates.tsv" || { rc=1; continue; }
  cnmf prepare --output-dir "${OUTDIR}" --name "${RUN}" \
      -c "${IN}/${RUN}/counts.h5ad" -k ${K} --n-iter ${NITER} --seed ${SEED} \
      --numgenes ${NUMGENES} --total-workers ${WORKERS} || { rc=1; continue; }
  for w in $(seq 0 $((WORKERS-1))); do
    cnmf factorize --output-dir "${OUTDIR}" --name "${RUN}" \
        --worker-index ${w} --total-workers ${WORKERS} &
  done
  wait
  cnmf combine --output-dir "${OUTDIR}" --name "${RUN}" || { rc=1; continue; }
  cnmf consensus --output-dir "${OUTDIR}" --name "${RUN}" --components ${K} \
      --local-density-threshold 0.10 || rc=1
done
echo "=== M40b end $(date) | rc=${rc} ==="
exit $rc
