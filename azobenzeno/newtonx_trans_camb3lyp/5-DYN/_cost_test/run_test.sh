#!/bin/bash
#SBATCH --job-name=nx_cost_test
#SBATCH --time=02:00:00
#SBATCH --output=%x_%a.out
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --array=0-1
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
D=(tzvp_tda svp_tda)
cd $SLURM_SUBMIT_DIR/${D[$SLURM_ARRAY_TASK_ID]}
for step in ridft escf egrad; do
  s=$(date +%s); $step > $step.out 2>&1; e=$(date +%s)
  echo "$step $((e-s)) s" >> timing.txt
done
