#!/bin/bash
#SBATCH -t 02:00:00
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH -n 8
#SBATCH --partition=short
#SBATCH -o slurm_dscfpop_%A_%a.out
#SBATCH -e slurm_dscfpop_%A_%a.error
#SBATCH -J dscf_pop

# Refaz so o `dscf -proper` (populacao Mulliken/SEN) em dirs ja convergidos,
# garantindo $pop no control antes (dscf/ridft o consomem).
# Uso: sbatch --array=1-N job_dscf_pop.sh <joblist_file>

set -u
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE
export PATH=$TURBODIR/scripts:$PATH
source "$TURBODIR/Config_turbo_env"
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS

WORK_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${1:?joblist}")
[ -d "$WORK_DIR" ] || { echo "ERROR: dir invalido"; exit 1; }
cd "$WORK_DIR"
grep -q '^\$pop' control || sed -i '/^\$end/i \$pop mulliken loewdin nbo paboon' control
dscf -proper > pop.out 2>&1
grep -qi "ended abnormally" pop.out && echo "WARNING: dscf -proper abnormal, ver pop.out"
grep -ci "mulliken\|shared electron" pop.out
