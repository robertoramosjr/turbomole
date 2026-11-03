#!/bin/bash
#SBATCH --job-name=tm_sh_native
#SBATCH --time=12:00:00
#SBATCH --output=slurm_%A_%a.out
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --array=0-9
ulimit -s unlimited
read idx dir tm < <(awk -v i=$SLURM_ARRAY_TASK_ID '$1==i' $SLURM_SUBMIT_DIR/tasks.txt)
if [ "$tm" = "801" ]; then export TURBODIR=$HOME/software/TURBOMOLE801; else export TURBODIR=$HOME/software/TURBOMOLE; fi
export PARA_ARCH=SMP PARNODES=$SLURM_CPUS_PER_TASK OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
cd $SLURM_SUBMIT_DIR/$dir
echo "TM=$tm TURBODIR=$TURBODIR host=$(hostname) $(grep -m1 'model name' /proc/cpuinfo)" > run_info.txt
start=$(date +%s)
jobex -ri -ex -md -c 100 > jobex.out 2>&1
rc=$?
echo "rc=$rc wall=$(( $(date +%s)-start )) s" >> run_info.txt
exit $rc
