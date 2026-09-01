#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH -n 8
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J pbe0_ir_uvvis

# Layer A spectroscopy step: IR (aoforce) + UV-Vis (escf/TD-DFT) at
# PBE0/def2-SVP (RI-J + RI-K) on r2SCAN-3c-optimized geometries.
# r2SCAN-3c itself cannot do this (confirmed: both aoforce and escf
# abort with "Invalid value of nfun in <mgga_r2>" -- meta-GGA analytic
# Hessian/TD-DFT response isn't implemented in this Turbomole build).
# See setup_pbe0_svp_template.py for why PBE0/def2-SVP was chosen over
# reusing r2SCAN-3c's own def2-mTZVP basis (RI-K's auxiliary basis
# assignment silently fails for def2-mTZVP).
#
# Runs ridft (SCF, needed to reconverge the wavefunction for the new
# functional/basis before any property calc) -> aoforce -> escf in
# sequence per array task. Walltime is generous (4h): a single-molecule
# timing test showed aoforce alone can run well over an hour at this
# level, driven by exact-exchange cost in the Hessian's response
# equations.
#
# Usage:
#   sbatch --array=1-327%N job_pbe0_ir_uvvis.sh <joblist_file>

set -u
ulimit -s unlimited

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

JOBLIST=${1:?"Usage: sbatch --array=1-N job_pbe0_ir_uvvis.sh <joblist_file>"}
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
if [ $ridft_status -ne 0 ]; then
    exit $ridft_status
fi

aoforce > aoforce.out 2>&1
aoforce_status=$?
echo "aoforce exit $aoforce_status for $WORK_DIR"

escf > escf.out 2>&1
escf_status=$?
echo "escf exit $escf_status for $WORK_DIR"

if [ $aoforce_status -ne 0 ] || [ $escf_status -ne 0 ]; then
    exit 1
fi
exit 0
