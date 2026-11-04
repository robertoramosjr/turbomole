#!/bin/bash
#SBATCH -t 08:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --partition=short
#SBATCH --mem=64G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J diffuse_td

# HSE06-D3(BJ)/def2-TZVPD + COSMO(water) ridft, then fast-response COSMO-TD-DFT escf (same keywords as the production
# water_tddft runs, 30 roots).  Dirs from prepare_diffuse_test.py.  Resumable: ridft is skipped if RIDFT_DONE exists.
# Validate by TEXT: "ridft : all done" / "escf : all done"; markers RIDFT_DONE / ESCF_DONE.
# Usage: sbatch --array=1-N%6 job_diffuse_tddft.sh <joblist_file>   (absolute dirs, one per line)
set -u
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE; source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
JOBLIST=${1:?"Usage: sbatch --array=1-N job_diffuse_tddft.sh <joblist_file>"}
D=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
[ -n "$D" ] && [ -d "$D" ] || { echo "ERROR: no valid dir on line $SLURM_ARRAY_TASK_ID of $JOBLIST"; exit 1; }
cd "$D"
export PARA_ARCH=SMP PARNODES=$SLURM_NTASKS OMP_NUM_THREADS=$SLURM_NTASKS OMP_STACKSIZE=4G OMP_PLACES=cores OMP_PROC_BIND=close
if [ ! -f RIDFT_DONE ]; then
  ridft > ridft.out 2>&1
  if grep -q "ridft : all done" ridft.out && ! grep -qi "ended abnormally\|did not converge" ridft.out; then echo RIDFT_DONE > RIDFT_DONE
  else echo "FAIL ridft in $D"; exit 1; fi
fi
grep "Total energy + OC\|HOMO-LUMO" ridft.out | tail -2
if [ ! -f ESCF_DONE ]; then
  escf > escf.out 2>&1
  # escf prints "all done" even when the Davidson solver did not converge: also require no "No convergence" and an exspectrum file
  grep -q "escf : all done" escf.out && ! grep -qi "ended abnormally\|No convergence" escf.out && [ -s exspectrum ] && echo ESCF_DONE > ESCF_DONE
fi
echo "diffuse TD-DFT $D escf_done=$([ -f ESCF_DONE ] && echo 1 || echo 0)"
