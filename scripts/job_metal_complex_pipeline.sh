#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --partition=medium
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J metal_complex_pipeline

# Full pipeline_extracao_dft.md Estagio 1 (production), one SLURM array
# task per directory listed in <joblist_file> (one path per line):
#
#   jobex (full geometry optimization, no cycle cap -- production
#          $statpt thresholds already baked into each control by
#          setup_metal_complex_template.py)
#   -> ridft -proper (density cube for Bader)
#   -> bader (Bader charges)
#   -> remove leftover $pointval from control (otherwise every later
#      dscf/ridft -proper call recomputes the full cube grid -- known
#      slowdown, already hit once for azo_trans)
#   -> ridft (tight reconverge before TD-DFT)
#   -> escf (TD-DFT, singlet states 1-N from the control's $soes)
#   -> aoforce (Hessian/IR)
#   -> dscf -proper (Mulliken/SEN population)
#
# Turbomole exit codes are not trustworthy on their own (a binary can
# print "did not converge" / "ended abnormally" and still return 0) --
# grep each .out for failure text in addition to checking $?.
#
# Per-stage wall-clock time is logged to <WORK_DIR>/timing.log (CSV:
# stage,start_epoch,end_epoch,elapsed_seconds,exit_status) so the three
# functionals (hse06/pbe/pbe0, identical everything else except the
# functional itself) can be compared as a benchmark once all 9 finish.
#
# Usage:
#   sbatch --array=1-N job_metal_complex_pipeline.sh <joblist_file>

set -u
ulimit -s unlimited

export TURBODIR=$HOME/software/TURBOMOLE
export PATH=$TURBODIR/scripts:$PATH
source "$TURBODIR/Config_turbo_env"
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS

JOBLIST=${1:?"Usage: sbatch --array=1-N job_metal_complex_pipeline.sh <joblist_file>"}
WORK_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")

if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: no valid directory on line $SLURM_ARRAY_TASK_ID of $JOBLIST"
    exit 1
fi

cd "$WORK_DIR"

TIMING_LOG="$WORK_DIR/timing.log"
if [ ! -f "$TIMING_LOG" ]; then
    echo "stage,start_epoch,end_epoch,elapsed_seconds,exit_status" > "$TIMING_LOG"
fi

# Runs one pipeline stage, timing it and appending a row to timing.log.
# $1 = stage name (also used as the label in stdout)
# $2 = the shell command to run for this stage (as a single string)
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

echo "=== $(date) starting pipeline in $WORK_DIR (PARNODES=$PARNODES) ==="

run_stage "jobex" "jobex > jobex.out 2>&1"
jobex_status=$?
if grep -qi "OPTIMIZATION DID NOT CONVERGE\|abnormally" jobex.out; then
    echo "WARNING: jobex reported non-convergence/abnormal exit, see jobex.out"
fi
if [ ! -f coord ]; then
    echo "FATAL: jobex left no usable coord, aborting pipeline for $WORK_DIR"
    exit 1
fi
if [ -f GEO_OPT_CONVERGED ]; then
    echo "jobex: geometry converged"
else
    echo "WARNING: jobex did not reach GEO_OPT_CONVERGED (walltime cap or non-convergence) -- proceeding with the last available geometry anyway"
fi

run_stage "ridft_proper" "ridft -proper > proper.out 2>&1"
if grep -qi "ended abnormally" proper.out; then
    echo "FATAL: ridft -proper ended abnormally, see proper.out"
    exit 1
fi

BADER_BIN=$(command -v bader || echo ~/local_bin/bader)
run_stage "bader" "\"$BADER_BIN\" td.cub > bader_run.out 2>&1"
tail -5 ACF.dat

if grep -q '^\$pointval' control; then
    echo "--- removing leftover \$pointval from control ---"
    sed -i '/^\$pointval/,/^origin/d' control
fi

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

run_stage "dscf_proper" "dscf -proper > pop.out 2>&1"
if grep -qi "ended abnormally" pop.out; then
    echo "WARNING: dscf -proper ended abnormally, see pop.out"
fi

echo "=== $(date) pipeline finished in $WORK_DIR ==="
echo "--- timing summary ($WORK_DIR) ---"
column -s, -t "$TIMING_LOG"
