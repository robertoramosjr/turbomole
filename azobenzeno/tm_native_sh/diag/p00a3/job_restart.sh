#!/bin/bash
#SBATCH -J p00a3_restart
#SBATCH -t 00:20:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o /home/rramos/work_turbomole/azobenzeno/tm_native_sh/diag/p00a3/slurm_%A_%a.out
# P00a3: is the control truncation (p00a2) a property of the 256-snapshot or of the restart route (nextstep=ridft -> dscf-like step + egrad)?
# Usage: sbatch --array=1-3 job_restart.sh
#   1=young_ridft (127-frame traj, nextstep=ridft = same route as p00a2)
#   2=snap_full_manual   (256 snapshot, ridft run by hand, nextstep=egrad)
#   3=snap_pruned_manual (same, $nac/$excitation pruned to 2 cycles)
set -u
ulimit -s unlimited
cd /home/rramos/work_turbomole/azobenzeno/tm_native_sh/diag/p00a3
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" <<< $'young_ridft\nsnap_full_manual\nsnap_pruned_manual')
export TURBODIR=$HOME/software/TURBOMOLE
export PARA_ARCH=SMP PARNODES=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
cd $DIR || exit 1
echo "TM=7.8.1 host=$(hostname) $(date)" > run_info.txt
n0=$(grep -c '^t=' mdlog.1)
l0=$(wc -l < control)
case $DIR in
  *_manual)
    ridft > ridft.out 2>&1
    grep -q "ridft : all done" ridft.out || { echo "ridft FAILED $(date)" >> run_info.txt; exit 1; }
    echo egrad > nextstep
    echo "control_lines_after_ridft=$(wc -l < control)" >> run_info.txt ;;
esac
jobex -ri -ex -md > jobex.out 2>&1
echo "done $(date) frames_before=$n0 frames_after=$(grep -c '^t=' mdlog.1) control_lines_before=$l0 control_lines_after=$(wc -l < control)" >> run_info.txt
