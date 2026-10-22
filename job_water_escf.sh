#!/bin/bash
#SBATCH -t 03:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J water_escf

# COSMO-TDDFT escf (HSE06-D3(BJ)/def2-TZVP, eps 78.4, fast response, 60 roots) on W4 geometries/orbitals.
# Dirs from prepare_water_escf.py (converged COSMO ridft orbitals already in mos, so no ridft here).
# Validate by TEXT: "escf : all done" in escf.out; marker ESCF_DONE.
# Usage: sbatch --array=1-N%27 job_water_escf.sh <joblist_file>   (absolute dirs, one per line)
set -u
ulimit -s unlimited
JOBLIST=${1:?"Usage: sbatch --array=1-N job_water_escf.sh <joblist_file>"}
WORK_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
[ -n "$WORK_DIR" ] && [ -d "$WORK_DIR" ] || { echo "ERROR: no valid dir on line $SLURM_ARRAY_TASK_ID of $JOBLIST"; exit 1; }
cd "$WORK_DIR"
export PARA_ARCH=SMP PARNODES=$SLURM_NTASKS OMP_NUM_THREADS=$SLURM_NTASKS OMP_STACKSIZE=4G OMP_PLACES=cores OMP_PROC_BIND=close
escf > escf.out 2>&1
grep -q "escf : all done" escf.out && ! grep -qi "ended abnormally" escf.out && echo ESCF_DONE > ESCF_DONE
echo "escf for $WORK_DIR done=$([ -f ESCF_DONE ] && echo 1 || echo 0)"
