#!/bin/bash
#SBATCH --job-name=p5_m37
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/phase5/M37/m37_%j.out
#SBATCH --error=logs/phase5/M37/m37_%j.err
# =============================================================================
# Phase 5 - M37 - continuous malignant program discovery with cNMF.
#
# Three runs, all declared before any factor was seen:
#   primary   6,434 malignant cells, full declared gene universe
#   balanced  patient-balanced downsample (sensitivity)
#   nocc      canonical cell-cycle genes removed (sensitivity)
#
# 48 G / 8 CPUs: the factorized matrix is 6,434 cells x 2,000 over-dispersed
# genes - small. The cost is 12 K values x 100 replicate NMF fits, which is
# embarrassingly parallel, so 8 workers are run concurrently per run. This is
# an evidence-based request, not the 450 G maximum.
#
# R_env prepares the matrices; the ISOLATED p5_cnmf_env runs cNMF. Neither
# environment is modified.
# =============================================================================
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase5/M37 results/phase5/programs/cnmf
echo "=== M37 start $(date) | Job ${SLURM_JOB_ID} | node $(hostname) ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh

KGRID="4 5 6 7 8 9 10 11 12 13 14 15"
NITER=100
SEED=42
NUMGENES=2000
WORKERS=8
OUTDIR=results/phase5/programs/cnmf
IN=results/phase5/programs/cnmf_input

echo; echo "=== M37a: prepare matrices (R_env) ==="
conda activate R_env
Rscript scripts/phase5/programs/m37_prepare_input.R; rc=$?
conda deactivate
[ $rc -ne 0 ] && { echo "M37a failed rc=$rc"; exit $rc; }

echo; echo "=== M37b: build AnnData + run cNMF (p5_cnmf_env) ==="
conda activate p5_cnmf_env
python -c "import cnmf,importlib.metadata as m;print('cnmf',m.version('cnmf'))"

for RUN in primary balanced nocc; do
  echo; echo "----------------------------------------------------------------"
  echo "RUN=${RUN}  $(date)"
  echo "----------------------------------------------------------------"
  H5=${IN}/${RUN}/counts.h5ad
  python scripts/phase5/programs/m37_build_anndata.py \
      "${IN}/${RUN}" "${H5}" "${IN}/malignant_cell_covariates.tsv" || exit 1

  cnmf prepare --output-dir "${OUTDIR}" --name "${RUN}" -c "${H5}" \
      -k ${KGRID} --n-iter ${NITER} --seed ${SEED} --numgenes ${NUMGENES} \
      --total-workers ${WORKERS} || exit 1

  for w in $(seq 0 $((WORKERS-1))); do
    cnmf factorize --output-dir "${OUTDIR}" --name "${RUN}" \
        --worker-index ${w} --total-workers ${WORKERS} &
  done
  wait
  echo "factorize done for ${RUN}: $(date)"

  cnmf combine --output-dir "${OUTDIR}" --name "${RUN}" || exit 1
  cnmf k_selection_plot --output-dir "${OUTDIR}" --name "${RUN}" || exit 1
  echo "combine + k_selection_plot done for ${RUN}"
done

echo; echo "=== M37 end $(date) ==="
