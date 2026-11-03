#!/bin/bash
#SBATCH --job-name=nx_ic_1000
#SBATCH --time=20:00:00
#SBATCH --output=I%a/job.out
#SBATCH --error=I%a/job.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

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

export NX=/home/rramos/software/newton-x-package-26a/bin

cd $SLURM_SUBMIT_DIR/I${SLURM_ARRAY_TASK_ID}
$NX/initcond.pl > initcond.log 2>&1
