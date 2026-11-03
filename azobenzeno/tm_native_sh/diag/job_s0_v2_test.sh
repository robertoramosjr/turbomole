#!/bin/bash
#SBATCH -J s0v2_test
#SBATCH -t 01:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o /home/rramos/work_turbomole/azobenzeno/tm_native_sh/diag/slurm_s0v2_%A.out
# Diagnostic (plan 2026-10-16, Task 3): S0 continuation of t3/t9 with generator v2 (no stale $actual step). Dirs = "$@"
set -u
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE PARA_ARCH=SMP PARNODES=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
for d in "$@"; do ( cd "$d" && echo "$(date) start" > run_info.txt && jobex -ri -md -c 25 > jobex.out 2>&1; echo "$(date) frames=$(grep -c '^t=' mdlog.1)" >> run_info.txt ); done
