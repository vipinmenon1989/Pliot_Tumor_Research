#!/bin/bash
#SBATCH --job-name=p4_m32cpdb
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/phase4/M32/m32_cpdb_%j.out
#SBATCH --error=logs/phase4/M32/m32_cpdb_%j.err
# Phase 4 · M32 — CellPhoneDB on the refined labels.
# Runs in the deliberately ISOLATED cpdb_env, exactly as in Phase 3, and reuses
# scripts/phase3/ccc/run_cellphonedb.py UNMODIFIED with the same database,
# iterations (1000), threshold (0.10) and seed (42).
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase4/M32
echo "=== M32-cpdb start $(date) | Job ${SLURM_JOB_ID} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate cpdb_env
which python && python --version
python -c "import cellphonedb, anndata, scanpy; print('cellphonedb', cellphonedb.__version__)" 2>/dev/null || true
python scripts/phase3/ccc/run_cellphonedb.py \
  --input-dir results/phase4/ccc_refinement/cellphonedb_inputs \
  --out-dir results/phase4/ccc_refinement/by_method \
  --db-dir external/cellphonedb \
  --iterations 1000 --threshold 0.10 --seed 42
rc=$?; echo "=== M32-cpdb end $(date) | exit ${rc} ==="; exit $rc
