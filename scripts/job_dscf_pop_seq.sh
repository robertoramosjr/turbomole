#!/bin/bash
#SBATCH -J dscf_pop
#SBATCH --partition=par48-x
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -t 02:00:00
#SBATCH -o slurm_dscfpop_%j.out

# Sequential version of job_dscf_pop.sh for coaraci (no job arrays there):
# runs `dscf -proper` (Mulliken/Loewdin/SEN population) on every directory of a
# joblist, one after the other. Validates by text, never by exit code.
# Usage: sbatch job_dscf_pop_seq.sh <joblist_file>   (one absolute dir per line)
set -u
ulimit -s unlimited
export PARA_ARCH=SMP PARNODES=${SLURM_CPUS_PER_TASK:-16} OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
source "$HOME/software/turbomole_env.sh"
export TURBOTMPDIR="${TMPDIR:-/tmp}/tm_pop_${SLURM_JOB_ID:-local}_$$"; mkdir -p "$TURBOTMPDIR"
trap 'rm -rf "$TURBOTMPDIR"' EXIT
while read -r W; do
    [ -d "$W" ] || { echo "ERROR: $W not a dir"; continue; }
    cd "$W" || continue
    grep -q '^\$pop' control || sed -i '/^\$end/i \$pop mulliken loewdin nbo paboon' control
    actual -r > /dev/null 2>&1
    t0=$(date +%s); dscf -proper > pop.out 2>&1
    echo "$(basename "$W"): $(( $(date +%s) - t0 )) s, ended normally: $(grep -ci 'ended normally' pop.out), abnormal: $(grep -ci 'ended abnormally' pop.out), mulliken rows: $(grep -ci 'mulliken\|shared electron' pop.out)"
done < "${1:?joblist}"
