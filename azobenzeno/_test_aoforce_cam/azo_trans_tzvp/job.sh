#!/bin/bash
#SBATCH --job-name=aof_cam_test
#SBATCH --time=04:00:00
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

# Smoke/timing test: aoforce CAM-B3LYP-D3(BJ)/def2-TZVP + $senex, TM 8.0.1,
# on the (non-stationary) HSE06 trans geometry. Imaginary modes are expected;
# the goal is to prove the run completes and to measure wall-time.
# Usage: sbatch job.sh   (from this directory)

set -e
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE801   # TM 8.0.1 (NOT ~/software/TURBOMOLE = 7.8.1)
export PATH=$TURBODIR/scripts:$PATH
export PARA_ARCH=SMP
export PARNODES=$SLURM_CPUS_PER_TASK
source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/bin/$(sysname):$PATH
which ridft aoforce | grep -q TURBOMOLE801 || { echo "wrong TURBODIR"; exit 1; }
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd $SLURM_SUBMIT_DIR
ridft > ridft.out
grep -q "ridft : all done" ridft.out
aoforce > aoforce.out
grep -q "force : all done" aoforce.out
echo "TEST_OK"
