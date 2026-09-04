#!/bin/bash
# Submit the M30 -> M35 chain with SLURM dependencies. Run after M29 completes.
# Each stage starts only if its predecessor exited 0 (afterok), so a failure
# stops the chain instead of propagating bad inputs downstream.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
J30=$(sbatch --parsable scripts/shell/phase4/run_m30.sh)
echo "M30 (M29 aggregation + malignancy integration + figures) = $J30"
JP=$(sbatch  --parsable --dependency=afterok:$J30 scripts/shell/phase4/run_m32_prepare.sh)
echo "M32-prepare (refined CCC input)                          = $JP  (afterok:$J30)"
JL=$(sbatch  --parsable --dependency=afterok:$JP scripts/shell/phase4/run_m32_liana.sh)
echo "M32-LIANA                                                = $JL  (afterok:$JP)"
JC=$(sbatch  --parsable --dependency=afterok:$JP scripts/shell/phase4/run_m32_cellchat.sh)
echo "M32-CellChat (array 1-4, per patient)                    = $JC  (afterok:$JP)"
JD=$(sbatch  --parsable --dependency=afterok:$JP scripts/shell/phase4/run_m32_cellphonedb.sh)
echo "M32-CellPhoneDB (isolated cpdb_env)                      = $JD  (afterok:$JP)"
JK=$(sbatch  --parsable --dependency=afterok:$JL:$JC:$JD scripts/shell/phase4/run_m32_concordance.sh)
echo "M32-concordance + sensitivity                            = $JK  (afterok:$JL:$JC:$JD)"
JT=$(sbatch  --parsable --dependency=afterok:$JK scripts/shell/phase4/run_m31_m33_m34.sh)
echo "M31 states + M33 state->TME + M34 robustness             = $JT  (afterok:$JK)"
J35=$(sbatch --parsable --dependency=afterok:$JT scripts/shell/phase4/run_m35.sh)
echo "M35 final object + validation + freeze                   = $J35 (afterok:$JT)"
echo
echo "CHAIN_JOBIDS M30=$J30 M32prep=$JP M32liana=$JL M32cellchat=$JC M32cpdb=$JD M32conc=$JK M31_33_34=$JT M35=$J35"
