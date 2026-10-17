#!/bin/bash
# Idempotent, restartable stage runner for ONE metal-complex directory.
# Designed to be launched (in the background) by job_metal_complex_trio.sh,
# which runs several of these side by side on a single whole node.
#
# Stages (each is skipped if already validated, so re-running after a
# walltime kill simply continues where it stopped):
#
#   jobex        -> done when GEO_OPT_CONVERGED exists
#   ridft_proper -> proper.out "ended normally", newer than GEO_OPT_CONVERGED
#   bader        -> ACF.dat newer than GEO_OPT_CONVERGED (skipped with a
#                   warning if no bader binary is installed)
#   ridft_tight  -> ridft_tight.out "ended normally", newer than GEO_OPT_CONVERGED
#   escf         -> escf.out "ended normally", newer than GEO_OPT_CONVERGED
#   numforce     -> vibspectrum newer than GEO_OPT_CONVERGED and
#                   numforce/aoforce.out "ended normally"
#                   (NumForce instead of aoforce: aoforce has no usable
#                   checkpoint and did not fit 3 d; NumForce restarts
#                   automatically from numforce/KraftWerk)
#   dscf_pop     -> pop.out "ended normally" + Mulliken rows, newer than
#                   GEO_OPT_CONVERGED
#
# Exit codes: 0 = all stages validated, 2 = stopped early but
# restartable (e.g. jobex hit its cycle cap), 1 = FATAL (needs a human;
# a .stages/FAILED marker is written so the trio job stops resubmitting).
#
# Turbomole exit codes are not trustworthy: every stage is validated by
# the text of its output, never by $? alone.
#
# Usage: run_metal_complex_stages.sh <work_dir> <ncores>

set -u
ulimit -s unlimited

WORK_DIR=$(readlink -f "${1:?Usage: run_metal_complex_stages.sh <work_dir> <ncores>}")
NCORES=${2:?Usage: run_metal_complex_stages.sh <work_dir> <ncores>}
JOBEX_CYCLES=${JOBEX_CYCLES:-1500}
NUMFORCE_OPTS=${NUMFORCE_OPTS:--central}

cd "$WORK_DIR" || exit 1
mkdir -p .stages
rm -f .stages/ALL_DONE

export PARA_ARCH=SMP
export PARNODES=$NCORES
export OMP_NUM_THREADS=$NCORES
source "${TURBOMOLE_ENV:-$HOME/software/turbomole_env.sh}"

# Private scratch per directory: turbomole_env.sh sets TURBOTMPDIR=/tmp, and
# ridft keeps diff_densmat/olddensmat there under FIXED names -- several runs
# on the same node read each other's files (seen 2026-10-14: Cu UHF ridft
# crashed with "forrtl severe (67) ... /tmp/diff_densmat" written by Ni).
export TURBOTMPDIR="${TMPDIR:-/tmp}/tm_${SLURM_JOB_ID:-local}_$(basename "$(dirname "$WORK_DIR")")_$$"
mkdir -p "$TURBOTMPDIR"
trap 'rm -rf "$TURBOTMPDIR"' EXIT

TIMING_LOG="$WORK_DIR/timing.log"
log() { echo "[$(date '+%F %T')] [$(basename "$(dirname "$WORK_DIR")")] $*"; }

fatal() {
    log "FATAL: $*"
    echo "$(date '+%F %T') $*" > .stages/FAILED
    exit 1
}

run_stage() {
    local stage="$1" cmd="$2" t0 t1 status
    log "--- $stage (PARNODES=$PARNODES) ---"
    t0=$(date +%s)
    eval "$cmd"
    status=$?
    t1=$(date +%s)
    echo "$stage,$t0,$t1,$((t1 - t0)),$status" >> "$TIMING_LOG"
    log "$stage: exit $status, $((t1 - t0)) s, node=$(hostname)"
    return $status
}

# file exists, is newer than the converged geometry and reports success
fresh_ok() {
    local file="$1" pattern="${2:-ended normally}"
    [ -f GEO_OPT_CONVERGED ] && [ -f "$file" ] && [ "$file" -nt GEO_OPT_CONVERGED ] \
        && grep -qi "$pattern" "$file"
}

# --- recovery from a previous kill -----------------------------------------
# A kill during ridft -proper leaves the $pointval grid block in control.
if [ -f control.pre_proper ]; then
    log "restoring control.pre_proper left by an interrupted ridft -proper"
    mv control.pre_proper control
fi
# A kill mid-program leaves "$actual step <x>", which blocks every later run.
actual -r > /dev/null 2>&1
grep -q '^\$end' control || fatal "control has no \$end (truncated?)"

# --- jobex -------------------------------------------------------------------
if [ ! -f GEO_OPT_CONVERGED ]; then
    rm -f GEO_OPT_RUNNING GEO_OPT_FAILED not.converged
    run_stage "jobex" "jobex -c $JOBEX_CYCLES > jobex.out 2>&1"
    if [ ! -f GEO_OPT_CONVERGED ]; then
        log "jobex did not converge in this round (last: $(tail -2 energy | head -1 | awk '{print "cycle "$1", E="$2}')) -- restartable"
        exit 2
    fi
fi
log "jobex: GEO_OPT_CONVERGED"

# --- ridft -proper + bader -----------------------------------------------------
if ! fresh_ok proper.out; then
    cp control control.pre_proper
    run_stage "ridft_proper" "ridft -proper > proper.out 2>&1"
    grep -qi "ended normally" proper.out || fatal "ridft -proper failed, see proper.out"
fi
if [ ! -f ACF.dat ] || [ ! ACF.dat -nt GEO_OPT_CONVERGED ]; then
    BADER_BIN=$(command -v bader || { [ -x "$HOME/local_bin/bader" ] && echo "$HOME/local_bin/bader"; } || true)
    if [ -n "$BADER_BIN" ] && [ -f td.cub ]; then
        run_stage "bader" "\"$BADER_BIN\" td.cub > bader_run.out 2>&1"
    else
        log "WARNING: bader skipped (binary or td.cub missing); td.cub kept for later"
    fi
fi
[ -f control.pre_proper ] && mv control.pre_proper control

# --- ridft (tight) + escf ------------------------------------------------------
if ! fresh_ok ridft_tight.out; then
    run_stage "ridft_tight" "ridft > ridft_tight.out 2>&1"
    fresh_ok ridft_tight.out || fatal "ridft (tight) failed, see ridft_tight.out"
fi
if ! fresh_ok escf.out; then
    run_stage "escf" "escf > escf.out 2>&1"
    fresh_ok escf.out || log "WARNING: escf did not end normally, see escf.out (continuing)"
    actual -r > /dev/null 2>&1
fi

# --- NumForce (numerical Hessian, restartable) --------------------------------
if ! { [ -f vibspectrum ] && [ vibspectrum -nt GEO_OPT_CONVERGED ] \
       && grep -qi "ended normally" numforce/aoforce.out 2>/dev/null; }; then
    nf_opts="$NUMFORCE_OPTS"
    if [ -d numforce/KraftWerk ]; then
        # reference gradient already checked in an earlier round: skip it (-c)
        nf_opts="$nf_opts -c"
        log "NumForce: restarting from numforce/KraftWerk ($(ls numforce/KraftWerk 2>/dev/null | wc -l) entries)"
        [ -f numforce.out ] && mv numforce.out "numforce.out.$(date +%s)"
    fi
    # NumForce runs PARNODES independent serial displacement jobs in SMP mode
    run_stage "numforce" "NumForce $nf_opts > numforce.out 2>&1"
    if ! { [ -f vibspectrum ] && [ vibspectrum -nt GEO_OPT_CONVERGED ]; }; then
        grep -qi "error\|abnormal" numforce.out && fatal "NumForce reported an error, see numforce.out"
        log "NumForce not finished in this round -- restartable"
        exit 2
    fi
fi
log "NumForce done; imaginary modes: $(awk '/^ *[0-9]+ +a/ && $3 < 0 {n++} END {print n+0}' vibspectrum)"

# --- dscf -proper population ---------------------------------------------------
if ! fresh_ok pop.out || ! grep -qi "mulliken" pop.out; then
    grep -q '^\$pop' control || sed -i '/^\$end/i \$pop mulliken loewdin nbo paboon' control
    run_stage "dscf_pop" "dscf -proper > pop.out 2>&1"
    fresh_ok pop.out || fatal "dscf -proper failed, see pop.out"
fi

touch .stages/ALL_DONE
log "=== all stages validated ==="
column -s, -t "$TIMING_LOG" | tail -10
exit 0
