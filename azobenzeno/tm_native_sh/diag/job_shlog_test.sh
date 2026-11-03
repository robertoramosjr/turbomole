#!/bin/bash
#SBATCH -J shlog_test
#SBATCH -t 03:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o /home/rramos/work_turbomole/azobenzeno/tm_native_sh/diag/slurm_shlog_%A_%a.out
# Diagnostic (plan 2026-10-19, Task 2): do the SH log groups, pre-seeded with a small maxlog, roll over and truncate the control?
# Usage: sbatch --array=1-3 job_shlog_test.sh <tasks: "<dir> <ncycles> <78|801>" per line>
set -u
ulimit -s unlimited
read DIR NC TM < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$1")
if [ "$TM" = "801" ]; then export TURBODIR=$HOME/software/TURBOMOLE801; else export TURBODIR=$HOME/software/TURBOMOLE; fi
export PARA_ARCH=SMP PARNODES=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
cd "$DIR" || exit 1
echo "TM=$TM host=$(hostname) $(date)" > run_info.txt
jobex -ri -ex -md -c $NC > jobex.out 2>&1
echo "done $(date) frames=$(grep -c '^t=' mdlog.1) control_lines=$(wc -l < control)" >> run_info.txt
