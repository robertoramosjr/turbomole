#!/bin/bash
#SBATCH -J mc_stab
#SBATCH --partition=par48-x
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH -o slurm_mc_stab_%j.out

# Wavefunction stability tests (ridft from the converged MOs, then escf $scfinstab
# singlet|triplet|urpa with $soes all 1), prepared by setup_stability_test.sh.
# All dirs run in parallel on one node: 48 / N cores each, private TURBOTMPDIR.
# Usage: sbatch job_stability_test.sh <dir1> [<dir2> ...]
# Validation is by TEXT (exit code is always 0); the eigenvalues are NOT parsed here:
# read the printed lines in escf.out, because the format was not verified yet.
set -u
[ $# -ge 1 ] || { echo "Usage: sbatch $0 <dir1> [<dir2> ...]"; exit 1; }
N=$#; NC=$(( ${SLURM_CPUS_PER_TASK:-48} / N ))
source "${TURBOMOLE_ENV:-$HOME/software/turbomole_env.sh}"
export PARA_ARCH=SMP PARNODES=$NC OMP_NUM_THREADS=$NC

run_one() {
    local d=$1 first=$2
    cd "$d" || return 1
    export TURBOTMPDIR="${TMPDIR:-/tmp}/tm_stab_${SLURM_JOB_ID:-local}_$(basename "$d")_$$"
    mkdir -p "$TURBOTMPDIR"
    rm -f STAB_DONE STAB_FAILED
    echo "$(date '+%F %T') ridft start ($d, $NC cores)"
    taskset -c "$first-$((first+NC-1))" ridft > ridft_ref.out 2>&1
    if ! grep -qi "ridft ended normally" ridft_ref.out || grep -qi "ended abnormally\|did not converge" ridft_ref.out; then
        echo "ridft FAILED" > STAB_FAILED; rm -rf "$TURBOTMPDIR"; return 1
    fi
    echo "$(date '+%F %T') escf start"
    taskset -c "$first-$((first+NC-1))" escf > escf_stab.out 2>&1
    if grep -qi "escf ended normally" escf_stab.out && ! grep -qi "ended abnormally" escf_stab.out; then
        date > STAB_DONE
    else
        echo "escf FAILED" > STAB_FAILED
    fi
    rm -rf "$TURBOTMPDIR"
    echo "$(date '+%F %T') finished $d: $(ls STAB_DONE STAB_FAILED 2>/dev/null)"
}

i=0
for d in "$@"; do
    run_one "$(readlink -f "$d")" $((i*NC)) > "$(readlink -f "$d")/stab.log" 2>&1 &
    i=$((i+1))
done
wait
for d in "$@"; do echo "== $d: $(ls "$d"/STAB_DONE "$d"/STAB_FAILED 2>/dev/null | xargs -n1 basename | tr '\n' ' ')"; done
