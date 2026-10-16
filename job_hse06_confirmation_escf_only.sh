#!/bin/bash
#SBATCH -t 03:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J hse06_confirm

# escf-only variant of job_hse06_labeling_retry2.sh, for B2's DFT
# confirmation batch (protonation_and_conformer_investigation.md Sec.
# 9.22): skips aoforce (B1's Hessian, ~27 core-hours/conformer, not
# needed here) -- only ridft (SCF) + escf (TD-DFT), matching the
# ~7.8 core-hour/conformer cost estimate reported before this batch
# was approved. Geometry is the already-converged r2SCAN-3c structure
# (no re-optimization at this stage -- see setup_hse06_tzvp_template.py).

set -u
ulimit -s unlimited

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

JOBLIST=${1:?"Usage: sbatch --array=1-N job_hse06_confirmation_escf_only.sh <joblist_file>"}
WORK_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")

if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: no valid directory on line $SLURM_ARRAY_TASK_ID of $JOBLIST"
    exit 1
fi

cd "$WORK_DIR"
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS
export OMP_NUM_THREADS=$SLURM_NTASKS
export OMP_STACKSIZE=4G
export OMP_PLACES=cores
export OMP_PROC_BIND=close

ridft > ridft.out 2>&1
ridft_status=$?
echo "ridft exit $ridft_status for $WORK_DIR"
grep -qi "ended abnormally\|did not converge" ridft.out && ridft_status=1
if [ $ridft_status -ne 0 ]; then
    exit $ridft_status
fi

escf > escf.out 2>&1
escf_status=$?
grep -qi "ended abnormally" escf.out && escf_status=1
echo "escf exit $escf_status for $WORK_DIR"

exit $escf_status
