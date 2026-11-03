#!/bin/bash
#SBATCH -J tm_s0_cont
#SBATCH -t 04:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --array=1-10
#SBATCH -o slurm_s0cont_%A_%a.out
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE; source $TURBODIR/Config_turbo_env
export PARA_ARCH=SMP PARNODES=16 OMP_NUM_THREADS=16
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
D=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SLURM_SUBMIT_DIR/s0_cont/tasks.txt)
cd $D && jobex -ri -md -c 305 > jobex.out 2>&1
python3 /home/rramos/work_turbomole/scripts/classify_tm_md.py mdlog.1 > result.txt 2>&1
