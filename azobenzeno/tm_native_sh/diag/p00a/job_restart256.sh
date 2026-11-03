#!/bin/bash
#SBATCH -J p00a_restart
#SBATCH -t 00:40:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o /home/rramos/work_turbomole/azobenzeno/tm_native_sh/diag/p00a/slurm_restart_%A_%a.out
# P00a test 1/2: restart from TRAJ2 snapshot at cycle 256. Ciclos limitados por $nsteps no mdmaster (4), NÃO por jobex -c.
# Usage: sbatch --array=1-2 job_restart256.sh   (1=full, 2=pruned)
set -u
ulimit -s unlimited
cd /home/rramos/work_turbomole/azobenzeno/tm_native_sh/diag/p00a
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" <<< $'full\npruned')
export TURBODIR=$HOME/software/TURBOMOLE
export PARA_ARCH=SMP PARNODES=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
cd $DIR || exit 1
echo "TM=7.8 host=$(hostname) $(date)" > run_info.txt
n0=$(grep -c '^t=' mdlog.1)
jobex -ri -ex -md -c 20 > jobex.out 2>&1
echo "done $(date) frames_before=$n0 frames_after=$(grep -c '^t=' mdlog.1) control_lines=$(wc -l < control)" >> run_info.txt
