#!/bin/bash
#SBATCH -t 02:00:00
#SBATCH -N 1 --ntasks-per-node=4
#SBATCH -n 4
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J ref_reopt

# r2SCAN-3c/def2-mTZVP + COSMO(water) re-optimisation of a pKa reference species (jobex -ri).
# Dirs from prepare_ref_water.py reopt.  Resumable: GEO_OPT_FAILED is cleared and the run continues.
# Validate by TEXT: GEO_OPT_CONVERGED.   Usage: sbatch --array=1-N job_ref_reopt.sh <joblist_file> [cycles]
set -u
source /opt/gridunesp/dist/miniconda/24.4.0/libmamba/etc/profile.d/conda.sh; conda activate paper2_funnel
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE; source $TURBODIR/Config_turbo_env
export PARA_ARCH=SMP PARNODES=$SLURM_NTASKS OMP_NUM_THREADS=$SLURM_NTASKS
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
JOBLIST=${1:?joblist}; CYC=${2:-300}
D=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
[ -n "$D" ] && [ -d "$D" ] || { echo "ERROR: no valid dir on line $SLURM_ARRAY_TASK_ID of $JOBLIST"; exit 1; }
cd "$D"
[ -f GEO_OPT_CONVERGED ] && { echo "already converged $D"; exit 0; }
if [ -f GEO_OPT_FAILED ]; then rm -f GEO_OPT_FAILED not.converged; cat jobex.out >> jobex.prev.out 2>/dev/null; fi
jobex -ri -c "$CYC" > jobex.out 2>&1
echo "jobex for $D conv=$([ -f GEO_OPT_CONVERGED ] && echo 1 || echo 0)"
