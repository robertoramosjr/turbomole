#!/bin/bash
#SBATCH -J p00a2_restart
#SBATCH -t 00:20:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o /home/rramos/work_turbomole/azobenzeno/tm_native_sh/diag/p00a2/slurm_%A_%a.out
# P00a2: restart TRAJ2 snapshot (cycle 256) with a COMPLETE control. Cycles limited by $nsteps (4) in mdmaster.
# Usage: sbatch --array=1-4 job_restart.sh   (1=full_78 2=pruned_78 3=full_801 4=pruned_801)
set -u
ulimit -s unlimited
cd /home/rramos/work_turbomole/azobenzeno/tm_native_sh/diag/p00a2
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" <<< $'full_78\npruned_78\nfull_801\npruned_801')
case $DIR in *_801) TMDIR=$HOME/software/TURBOMOLE801; V=8.0.1;; *) TMDIR=$HOME/software/TURBOMOLE; V=7.8.1;; esac
export TURBODIR=$TMDIR
export PARA_ARCH=SMP PARNODES=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
cd $DIR || exit 1
echo "TM=$V host=$(hostname) $(date)" > run_info.txt
n0=$(grep -c '^t=' mdlog.1)
jobex -ri -ex -md > jobex.out 2>&1
echo "done $(date) frames_before=$n0 frames_after=$(grep -c '^t=' mdlog.1) control_lines=$(wc -l < control)" >> run_info.txt
