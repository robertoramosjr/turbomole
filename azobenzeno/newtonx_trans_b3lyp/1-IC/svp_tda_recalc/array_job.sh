#!/bin/bash
#SBATCH --job-name=svp_tda_100
#SBATCH --time=00:30:00
#SBATCH --output=%a_job.out
#SBATCH --error=%a_job.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -e
ulimit -s unlimited

export TURBODIR=$HOME/software/TURBOMOLE
export PATH=$TURBODIR/scripts:$PATH
export PARA_ARCH=SMP
export PARNODES=$SLURM_CPUS_PER_TASK
source $TURBODIR/Config_turbo_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

POINT_IDS=($(ls -d P*/ | sed 's#P##; s#/##' | sort -n))
POINT_ID=${POINT_IDS[$SLURM_ARRAY_TASK_ID]}
cd "P${POINT_ID}"

ridft > ridft.out
escf > escf.out
