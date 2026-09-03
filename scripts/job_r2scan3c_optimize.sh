#!/bin/bash
#SBATCH -t 00:30:00
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH -n 8
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J r2scan3c_opt

# Layer A cheap-DFT refinement: runs one Turbomole jobex geometry
# optimization (r2SCAN-3c/def2-mTZVP) per SLURM array task, on the
# directory named by line SLURM_ARRAY_TASK_ID of a job-list file (one
# path per line, as written by list_conformer_jobs.sh).
#
# Stopping criterion is a HARD CYCLE CAP (jobex -c 10), not formal
# gradient convergence -- decided after an interactive test on this
# molecule showed the energy already flat to ~1e-6 Ha (five orders of
# magnitude below the ~1e-3 Ha/kcal-mol scale that matters for
# Boltzmann population weighting) by ~cycle 15-17, while the floppy
# prenyl groups kept the gradient-based convergence test from
# formally triggering even after 35+ cycles. This is a screening-stage
# choice, NOT used for the final production reoptimization, which
# keeps Paper 1's original tight $statpt thresholds untouched (see
# setup_r2scan3c_template.py for the loosened $statpt block used only
# here).
#
# Usage:
#   sbatch --array=1-327%30 job_r2scan3c_optimize.sh <joblist_file>

set -u
ulimit -s unlimited

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

JOBLIST=${1:?"Usage: sbatch --array=1-N job_r2scan3c_optimize.sh <joblist_file>"}
WORK_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")

if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: no valid directory on line $SLURM_ARRAY_TASK_ID of $JOBLIST"
    exit 1
fi

cd "$WORK_DIR"
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS

jobex -c 10 > jobex.out 2>&1
status=$?
echo "jobex exit $status for $WORK_DIR"
exit $status
