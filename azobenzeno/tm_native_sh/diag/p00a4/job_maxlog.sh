#!/bin/bash
#SBATCH -J p00a4_maxlog
#SBATCH -t 00:20:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o /home/rramos/work_turbomole/azobenzeno/tm_native_sh/diag/p00a4/slurm_%A_%a.out
# Usage: sbatch job_maxlog.sh
set -u
ulimit -s unlimited
cd /home/rramos/work_turbomole/azobenzeno/tm_native_sh/diag/p00a4
DIR=snap_maxlog1000
export TURBODIR=$HOME/software/TURBOMOLE
export PARA_ARCH=SMP PARNODES=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
cd $DIR || exit 1
echo "TM=7.8.1 host=$(hostname) $(date)" > run_info.txt
n0=$(grep -c '^t=' mdlog.1)
l0=$(wc -l < control)
case $DIR in
  snap_*)
    ridft > ridft.out 2>&1
    grep -q "ridft : all done" ridft.out || { echo "ridft FAILED $(date)" >> run_info.txt; exit 1; }
    echo egrad > nextstep
    echo "control_lines_after_ridft=$(wc -l < control)" >> run_info.txt ;;
esac
jobex -ri -ex -md > jobex.out 2>&1
echo "done $(date) frames_before=$n0 frames_after=$(grep -c '^t=' mdlog.1) control_lines_before=$l0 control_lines_after=$(wc -l < control)" >> run_info.txt
