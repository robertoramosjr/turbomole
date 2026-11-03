#!/bin/bash
#SBATCH --job-name=nx_traj_test
#SBATCH --time=01:30:00
#SBATCH --output=%x.out
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE
export PATH=$TURBODIR/scripts:$PATH
export PARA_ARCH=SMP
export PARNODES=$SLURM_CPUS_PER_TASK
source $TURBODIR/Config_turbo_env
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NX=$HOME/software/newton-x-package-26a/bin
export CIOVERLAP=$HOME/software/newton-x-package-26a/CIOVERLAP/bin
cd $SLURM_SUBMIT_DIR/TRAJ1
s=$(date +%s)
$NX/nx_moldyn --clean > md.log 2>&1
echo "nx_moldyn $(( $(date +%s)-s )) s" > ../timing.txt
