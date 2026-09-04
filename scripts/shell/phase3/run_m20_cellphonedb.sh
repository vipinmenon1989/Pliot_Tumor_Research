#!/bin/bash
#SBATCH --job-name=p3_m20_cpdb
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/phase3/M20/cellphonedb_%j.out
#SBATCH --error=logs/phase3/M20/cellphonedb_%j.err
# Phase 3 / M20 - CellPhoneDB, isolated env. Runs in parallel with the R-based methods
# because it is scientifically independent and touches a different conda environment.
set -uo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase3/M20 results/phase3/ccc/by_method external/cellphonedb
echo "=== M20 CellPhoneDB start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) | CPUs ${SLURM_CPUS_PER_TASK:-1} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate cpdb_env
which python && python --version
python scripts/phase3/ccc/run_cellphonedb.py \
  --input-dir results/phase3/ccc/cellphonedb_inputs \
  --out-dir results/phase3/ccc/by_method \
  --db-dir external/cellphonedb \
  --iterations 1000 --threshold 0.10 --seed 42
echo "=== End: $(date) ==="
