#!/bin/bash
#SBATCH --job-name=scan_s1_b3lyp_resume
#SBATCH --time=24:00:00
#SBATCH --output=job_resume.out
#SBATCH --error=job_resume.err
#SBATCH --partition=medium
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

cd $SLURM_SUBMIT_DIR
bash concatenated_scan_resume.sh
