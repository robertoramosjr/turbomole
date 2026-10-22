#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH -N 1 --ntasks-per-node=4
#SBATCH -n 4
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J xtb_hess

# xtb G_thermal (GFN2-xTB --hess, mRRHO at the input geometry; chosen over --bhess by the pilot, 2026-10-22).
# Dirs from prepare_xtb_gtherm_w4.py (mol.xyz + charge.txt).  Validate by TEXT: "G(RRHO) contrib." in hess.out.
# Usage: sbatch --array=1-N%27 job_xtb_hess.sh <joblist_file>   (absolute dirs, one per line)
set -u
JOBLIST=${1:?"Usage: sbatch --array=1-N job_xtb_hess.sh <joblist_file>"}
D=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
[ -n "$D" ] && [ -d "$D" ] || { echo "ERROR: no valid dir on line $SLURM_ARRAY_TASK_ID of $JOBLIST"; exit 1; }
cd "$D"
XTB=$HOME/software/TURBOMOLE/bin/em64t-unknown-linux-gnu/xtb
export OMP_NUM_THREADS=4 OMP_STACKSIZE=2G MKL_NUM_THREADS=4
ulimit -s unlimited
CHG=$(cat charge.txt)
mkdir -p hess
(cd hess && cp ../mol.xyz . && $XTB mol.xyz --gfn 2 --chrg $CHG --hess -P 4 > ../hess.out 2>&1)
ok=0; grep -q "G(RRHO) contrib" hess.out && ok=1
[ $ok = 1 ] && echo XTB_DONE > XTB_DONE
echo "xtb_hess for $D done=$ok"
