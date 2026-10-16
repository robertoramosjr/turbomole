#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH -N 1 --ntasks-per-node=4
#SBATCH -n 4
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J cosmo_ridft

# Runs the COSMO(water) single-point ridft for one conformer, on top
# of the already-converged (frozen) HSE06/def2-TZVP geometry --
# setup_cosmo_job.py already built $cosmo into control. ~0.52
# core-hours/conformer measured in a single-conformer trial
# (2026-09-14, job_0044): 31 min cpu-time on 4 cores.

set -u
ulimit -s unlimited

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

JOBLIST=${1:?"Usage: sbatch --array=1-N job_cosmo_ridft.sh <joblist_file>"}
WORK_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")

if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: no valid directory on line $SLURM_ARRAY_TASK_ID of $JOBLIST"
    exit 1
fi

cd "$WORK_DIR"
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS
export OMP_NUM_THREADS=$SLURM_NTASKS
export OMP_STACKSIZE=2G

ridft > ridft_cosmo.out 2>&1
ridft_status=$?
grep -qi "ended abnormally\|did not converge" ridft_cosmo.out && ridft_status=1
echo "ridft (COSMO) exit $ridft_status for $WORK_DIR"
exit $ridft_status
