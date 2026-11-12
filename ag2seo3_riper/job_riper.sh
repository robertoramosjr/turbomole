#!/bin/bash
# riper (TURBOMOLE 8.0.1, SMP/OpenMP) on harpia. Usage: sbatch [-p par16|par48a] --cpus-per-task=N --mem=<900M*N> job_riper.sh <rundir> <sp|opt>
# RAM is 900 MB/CPU (16 CPUs -> 14400M). par8 is busy with CP2K array 789; par16 allows one 16-CPU job at a time.
#SBATCH -J riper
#SBATCH -p par16
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=14400M
#SBATCH -t 1-00:00:00
#SBATCH -o riper_%j.out
#SBATCH -e riper_%j.err
#SBATCH --mail-user=robeerto.aguiar.ramos@gmail.com --mail-type=ALL
set -eu
DIR=${1:?rundir}; MODE=${2:?sp|opt}
module load modmat/turbomole/8.0.1
ulimit -s unlimited
export PARA_ARCH=SMP OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_STACKSIZE=512M
export PATH=$TURBODIR/bin/$(sysname):$PATH
cd "$DIR"
grep -n -E "functional|xctype" control          # stale $xctype from copied controls changes the functional
echo "node=$(hostname) threads=$OMP_NUM_THREADS mode=$MODE"
if [ "$MODE" = opt ]; then jobex -riper -c 200 > jobex.out 2>&1; else riper > riper.out 2>&1; fi
grep -a -E "ENERGY|FREE ENERGY|converged|ended normally|ERROR" riper.out job.last 2>/dev/null | tail -8
