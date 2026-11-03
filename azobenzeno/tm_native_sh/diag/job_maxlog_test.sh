#!/bin/bash
#SBATCH -J maxlog_test
#SBATCH -t 03:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o /home/rramos/work_turbomole/azobenzeno/tm_native_sh/diag/slurm_maxlog_%A.out
# Diagnostic (plan 2026-10-16, Task 2): does the S1 SH run die when the SH logs roll over at maxlog? Dir = $1
set -u
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE PARA_ARCH=SMP PARNODES=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
cd "$1" || exit 1
echo "TM78 host=$(hostname) $(date)" > run_info.txt
jobex -ri -ex -md -c 40 > jobex.out 2>&1
echo "done $(date) frames=$(grep -c '^t=' mdlog.1)" >> run_info.txt
