#!/bin/bash
#SBATCH --job-name=p5_cnmf_env
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/phase5/install_cnmf_env_%j.out
#SBATCH --error=logs/phase5/install_cnmf_env_%j.err
# =============================================================================
# Phase 5 - build the ISOLATED cNMF environment.
#
# cNMF (Kotliar et al., eLife 2019) is a Python tool. R_env is the frozen
# Phase 1-4 stack and is NOT touched: no package is installed, upgraded or
# downgraded in it. Everything cNMF needs lives in p5_cnmf_env, which nothing
# else in this project depends on.
# =============================================================================
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase5 reports/phase5/environment
echo "=== cNMF env build start $(date) | Job ${SLURM_JOB_ID} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh

ENVNAME=p5_cnmf_env
if conda env list | grep -qE "^${ENVNAME}\s"; then
  echo "${ENVNAME} already exists - not rebuilding"
else
  # python 3.11: cnmf's pinned scikit-learn/scipy stack resolves cleanly there.
  conda create -y -n "${ENVNAME}" -c conda-forge python=3.11 pip || exit 1
fi
conda activate "${ENVNAME}" || exit 1
python -V
pip install --no-input "cnmf" || exit 1

echo; echo "=== resolved versions ==="
python - <<'PY'
import importlib, sys
print("python", sys.version.split()[0])
for m in ["cnmf","numpy","scipy","pandas","sklearn","anndata","scanpy","fastcluster","matplotlib","yaml"]:
    try:
        mod = importlib.import_module(m)
        print(f"{m:14s} {getattr(mod,'__version__','(no __version__)')}")
    except Exception as e:
        print(f"{m:14s} MISSING ({type(e).__name__})")
PY
echo; echo "=== cnmf CLI ==="
which cnmf && cnmf --help 2>&1 | head -20

conda env export > reports/phase5/environment/p5_cnmf_env.yml
echo "=== cNMF env build end $(date) ==="
