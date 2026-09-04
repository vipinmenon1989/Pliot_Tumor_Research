#!/bin/bash
# Submit M32 -> M35 with dependencies. Run after M30 has produced
# results/phase4/malignancy/PHASE4_MALIGNANCY_CALLS.tsv.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
JP=$(sbatch  --parsable scripts/shell/phase4/run_m32_prepare.sh)
JL=$(sbatch  --parsable --dependency=afterok:$JP scripts/shell/phase4/run_m32_liana.sh)
JC=$(sbatch  --parsable --dependency=afterok:$JP scripts/shell/phase4/run_m32_cellchat.sh)
JD=$(sbatch  --parsable --dependency=afterok:$JP scripts/shell/phase4/run_m32_cellphonedb.sh)
JK=$(sbatch  --parsable --dependency=afterok:$JL:$JC:$JD scripts/shell/phase4/run_m32_concordance.sh)
JT=$(sbatch  --parsable --dependency=afterok:$JK scripts/shell/phase4/run_m31_m33_m34.sh)
J35=$(sbatch --parsable --dependency=afterok:$JT scripts/shell/phase4/run_m35.sh)
echo "CHAIN M32prep=$JP M32liana=$JL M32cellchat=$JC M32cpdb=$JD M32conc=$JK M31_33_34=$JT M35=$J35"
