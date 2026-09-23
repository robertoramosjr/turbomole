#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --partition=medium
#SBATCH -o slurm_pbe0cont_%j.out
#SBATCH -e slurm_pbe0cont_%j.error
#SBATCH -J pbe0_continue

# Self-chaining continuation of a single pbe0_m4_d3 jobex that is too
# expensive (RI-K) to finish geometry optimization within one 3-day
# `medium` SLURM allocation -- see job_metal_complex_pipeline.sh for
# the sibling hse06/pbe pipeline this is adapted from (same
# aoforce/escf/pop tail, just split off from jobex so a walltime kill
# mid-optimization doesn't waste the rest of the pipeline on an
# unconverged geometry).
#
# Each invocation: cd into WORK_DIR, run `jobex -c 300` (resumes
# automatically from whatever coord/gradient/mos/hessapprox history is
# already there -- a prior walltime kill leaves this state intact, no
# manual cleanup needed between rounds). If jobex reaches
# GEO_OPT_CONVERGED, run the rest of the property pipeline once and
# stop chaining. Otherwise (walltime kill or jobex's own cycle cap),
# resubmit this same script for the same WORK_DIR with one fewer
# round left in the budget, so it keeps going across multiple 3-day
# allocations without needing anyone to babysit it.
#
# ROUNDS_LEFT is a hard safety cap on this self-resubmission chain
# (default 20 => up to ~60 days of medium-partition wall time for one
# system) so a genuine stuck/broken case does not resubmit forever
# unnoticed on a shared cluster.
#
# Usage:
#   sbatch job_pbe0_continue.sh <WORK_DIR> [ROUNDS_LEFT]

set -u
ulimit -s unlimited

export TURBODIR=$HOME/software/TURBOMOLE
export PATH=$TURBODIR/scripts:$PATH
source "$TURBODIR/Config_turbo_env"
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS

WORK_DIR=${1:?"Usage: sbatch job_pbe0_continue.sh <WORK_DIR> [ROUNDS_LEFT]"}
ROUNDS_LEFT=${2:-20}
THIS_SCRIPT="$(readlink -f "$0")"

if [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: $WORK_DIR is not a directory"
    exit 1
fi

cd "$WORK_DIR"

TIMING_LOG="$WORK_DIR/timing.log"
if [ ! -f "$TIMING_LOG" ]; then
    echo "stage,start_epoch,end_epoch,elapsed_seconds,exit_status" > "$TIMING_LOG"
fi

run_stage() {
    local stage="$1"
    local cmd="$2"
    echo "--- $stage ---"
    local t0 t1
    t0=$(date +%s)
    eval "$cmd"
    local status=$?
    t1=$(date +%s)
    echo "$stage,$t0,$t1,$((t1 - t0)),$status" >> "$TIMING_LOG"
    echo "$stage: ${status} exit, $((t1 - t0))s (node=$(hostname), cores=$PARNODES)"
    return $status
}

echo "=== $(date) job_pbe0_continue in $WORK_DIR (round budget: $ROUNDS_LEFT, PARNODES=$PARNODES) ==="

if [ -f GEO_OPT_CONVERGED ]; then
    echo "GEO_OPT_CONVERGED already present -- skipping jobex, going straight to the property pipeline"
else
    run_stage "jobex" "jobex -c 300 > jobex.out 2>&1"
    if grep -qi "abnormally" jobex.out; then
        echo "WARNING: jobex reported an abnormal exit, see jobex.out"
    fi
    if [ ! -f coord ]; then
        echo "FATAL: jobex left no usable coord, aborting for $WORK_DIR (not resubmitting)"
        exit 1
    fi
fi

if [ ! -f GEO_OPT_CONVERGED ]; then
    echo "jobex not yet converged in $WORK_DIR after this round."
    if [ "$ROUNDS_LEFT" -le 1 ]; then
        echo "GIVING UP: round budget exhausted -- resubmit manually with a fresh budget if more time is warranted."
        exit 0
    fi
    NEXT=$((ROUNDS_LEFT - 1))
    echo "Resubmitting for another round ($NEXT round(s) left after this one) ..."
    sbatch "$THIS_SCRIPT" "$WORK_DIR" "$NEXT"
    exit 0
fi

echo "jobex: geometry converged -- running the rest of the property pipeline"

# snapshot/restore em vez de sed no $pointval (o sed apagava ate EOF e perdia $end)
cp control control.pre_proper

run_stage "ridft_proper" "ridft -proper > proper.out 2>&1"
if grep -qi "ended abnormally" proper.out; then
    echo "FATAL: ridft -proper ended abnormally, see proper.out"
    exit 1
fi

BADER_BIN=$(command -v bader || echo ~/local_bin/bader)
run_stage "bader" "\"$BADER_BIN\" td.cub > bader_run.out 2>&1"
tail -5 ACF.dat

echo "--- restoring pre-proper control (drops the leftover \$pointval grid block) ---"
mv control.pre_proper control

run_stage "ridft_tight" "ridft > ridft_tight.out 2>&1"
if grep -qi "ended abnormally\|did not converge" ridft_tight.out; then
    echo "FATAL: ridft (tight) failed to converge, see ridft_tight.out"
    exit 1
fi

run_stage "escf" "escf > escf.out 2>&1"
if grep -qi "ended abnormally" escf.out; then
    echo "WARNING: escf ended abnormally, see escf.out"
fi

run_stage "aoforce" "aoforce > aoforce.out 2>&1"
if grep -qi "ended abnormally" aoforce.out; then
    echo "WARNING: aoforce ended abnormally, see aoforce.out"
fi
grep -i "imaginary" aoforce.out

# dscf/ridft consomem $pop do control -- reinserir logo antes do dscf -proper
if ! grep -q '^\$pop' control; then
    sed -i '/^\$end/i \$pop mulliken loewdin nbo paboon' control
fi

run_stage "dscf_proper" "dscf -proper > pop.out 2>&1"
if grep -qi "ended abnormally" pop.out; then
    echo "WARNING: dscf -proper ended abnormally, see pop.out"
fi

echo "=== $(date) pbe0 pipeline finished in $WORK_DIR -- no further resubmission ==="
echo "--- timing summary ($WORK_DIR) ---"
column -s, -t "$TIMING_LOG"
