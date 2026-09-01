#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J hse06_label_r2

# Second retry variant of job_hse06_labeling.sh: doubled core count
# (8->16, nodes here have 56 cores, plenty of headroom) and explicit
# OpenMP env vars for Turbomole's SMP binaries (built as *_omp,
# confirmed OpenMP-threaded) -- PARNODES alone gave a good ~7-8x
# speedup on 8 cores in earlier testing on this cluster, but was not
# verified to be exploiting all threads optimally; setting
# OMP_NUM_THREADS explicitly (matching PARNODES, not exceeding it --
# exceeding it would reproduce the CREST/OpenBLAS oversubscription bug
# from earlier in this project) plus a generous OMP_STACKSIZE removes
# any doubt for the large Hessian arrays involved in a 46-atom
# HSE06/def2-TZVP calculation.

set -u
ulimit -s unlimited

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

JOBLIST=${1:?"Usage: sbatch --array=1-N job_hse06_labeling_retry2.sh <joblist_file>"}
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

aoforce > aoforce.out 2>&1
aoforce_status=$?
grep -qi "ended abnormally" aoforce.out && aoforce_status=1
echo "aoforce exit $aoforce_status for $WORK_DIR"

escf > escf.out 2>&1
escf_status=$?
grep -qi "ended abnormally" escf.out && escf_status=1
echo "escf exit $escf_status for $WORK_DIR"

if [ $aoforce_status -ne 0 ] || [ $escf_status -ne 0 ]; then
    exit 1
fi
exit 0
