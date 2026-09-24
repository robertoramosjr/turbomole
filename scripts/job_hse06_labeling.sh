#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH -n 8
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J hse06_label

# Layer B labeling step: production-level (HSE06-D3(BJ)/def2-TZVP,
# Paper 1 settings) IR (aoforce) + UV-Vis (escf, 60 singlet states) on
# the active-learning labeled subset, single-point on the r2SCAN-3c-
# optimized geometry (no reoptimization at this level -- that is
# reserved for Estagio 6's final top-N). This is the "expensive
# target" data B1/B2's GPR surrogates are trained against.
#
# Same structure as job_pbe0_ir_uvvis.sh: ridft -> aoforce -> escf, and
# the same caveat applies -- Turbomole binaries can exit 0 despite an
# internal failure ("ridft did not converge", "ended abnormally"),
# confirmed once already in this project. Always grep the output text
# too when validating a finished batch, don't trust exit codes alone.
#
# Usage:
#   sbatch --array=1-N%K job_hse06_labeling.sh <joblist_file>

set -u
ulimit -s unlimited

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

JOBLIST=${1:?"Usage: sbatch --array=1-N job_hse06_labeling.sh <joblist_file>"}
WORK_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")

if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: no valid directory on line $SLURM_ARRAY_TASK_ID of $JOBLIST"
    exit 1
fi

cd "$WORK_DIR"
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS

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
