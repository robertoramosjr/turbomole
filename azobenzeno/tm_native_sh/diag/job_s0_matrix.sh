#!/bin/bash
#SBATCH -J s0_matrix
#SBATCH -t 01:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o /home/rramos/work_turbomole/azobenzeno/tm_native_sh/diag/slurm_s0matrix_%A_%a.out
# Diagnostic (plan 2026-10-19, Task 3): S0 continuation variants. Usage: sbatch --array=1-13 job_s0_matrix.sh <tasks: "<dir> <ncycles> <omp78|serial78|omp801>">
set -u
ulimit -s unlimited
read DIR NC MODE < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$1")
if [ "$MODE" = "omp801" ]; then export TURBODIR=$HOME/software/TURBOMOLE801; else export TURBODIR=$HOME/software/TURBOMOLE; fi
if [ "$MODE" = "serial78" ]; then N=1; else N=16; fi
export PARNODES=$N OMP_NUM_THREADS=$N MKL_NUM_THREADS=$N
[ "$MODE" = "serial78" ] && unset PARA_ARCH || export PARA_ARCH=SMP
source $TURBODIR/Config_turbo_env
[ "$MODE" = "serial78" ] && unset PARA_ARCH
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
cd "$DIR" || exit 1
echo "MODE=$MODE TURBODIR=$TURBODIR PARA_ARCH=${PARA_ARCH:-unset} host=$(hostname) $(date)" > run_info.txt
jobex -ri -md -c $NC > jobex.out 2>&1
echo "done $(date) frames=$(grep -c '^t=' mdlog.1)" >> run_info.txt
