#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J noise_floor

# Numerical-noise-floor test for B1's target frequency: re-runs
# ridft+aoforce on an already-converged HSE06/def2-TZVP conformer with
# a tighter grid (m4->m5) and tighter SCF/density convergence
# (scfconv 8->9, denconv 1d-7->1d-8), same geometry, to see how much
# the extracted C=O stretch frequency moves from purely numerical
# settings alone -- not TD-DFT, not needed for this test.

set -u
ulimit -s unlimited

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

JOBLIST=${1:?"Usage: sbatch --array=1-N job_noise_floor_test.sh <joblist_file>"}
WORK_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")

cd "$WORK_DIR"
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS
export OMP_NUM_THREADS=$SLURM_NTASKS
export OMP_STACKSIZE=4G
export OMP_PLACES=cores
export OMP_PROC_BIND=close

ridft > ridft.out 2>&1
ridft_status=$?
grep -qi "ended abnormally\|did not converge" ridft.out && ridft_status=1
echo "ridft exit $ridft_status for $WORK_DIR"
if [ $ridft_status -ne 0 ]; then
    exit $ridft_status
fi

aoforce > aoforce.out 2>&1
aoforce_status=$?
grep -qi "ended abnormally" aoforce.out && aoforce_status=1
echo "aoforce exit $aoforce_status for $WORK_DIR"
exit $aoforce_status
