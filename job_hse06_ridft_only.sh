#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J hse06_ridft

# ridft-ONLY variant (no escf, no aoforce) -- just converges the
# HSE06/def2-TZVP wavefunction on the frozen r2SCAN-3c geometry, the
# minimum needed before setup_population_job.py can compute a Bader
# charge at the SAME level the B1 monoanion_a model (Sec. 9.12) was
# calibrated on (population_jobs reused each labeled conformer's own
# converged HSE06 mos, not r2SCAN-3c's -- see Sec. 9.30 Task 1).
# ~0.27 core-hours/conformer measured previously for this same step.

set -u
ulimit -s unlimited

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

JOBLIST=${1:?"Usage: sbatch --array=1-N job_hse06_ridft_only.sh <joblist_file>"}
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
grep -qi "ended abnormally\|did not converge" ridft.out && ridft_status=1
echo "ridft exit $ridft_status for $WORK_DIR"
exit $ridft_status
