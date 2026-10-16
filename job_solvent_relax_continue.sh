#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH -n 8
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH -o slurm_cont_%A_%a.out
#SBATCH -J solv_relax_c
# Continues the r2SCAN-3c/COSMO relaxation of a conformer that hit the cycle cap (removes the failure flags and
# runs jobex -ri for MAX more cycles). Usage: sbatch --array=1-N job_solvent_relax_continue.sh <joblist> [max_cycles]
set -u
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE; source $TURBODIR/Config_turbo_env
export PARA_ARCH=SMP PARNODES=$SLURM_NTASKS
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
JOBLIST=${1:?joblist}; CYC=${2:-40}
read SRC OUT < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
cd "$OUT" || exit 1
rm -f GEO_OPT_FAILED not.converged
jobex -ri -c "$CYC" >> jobex.out 2>&1
echo "jobex exit $? for $OUT"
