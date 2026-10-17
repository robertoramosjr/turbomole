#!/bin/bash
#SBATCH -J mc_trio
#SBATCH --partition=par48-x
#SBATCH -t 3-00:00:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH -o slurm_mc_trio_%j.out

# Runs several metal-complex directories IN PARALLEL on ONE whole node
# (coaraci allocates whole nodes, so Slurm arrays waste cores). Each
# pending directory gets CORES_TOTAL / N_pending cores (3 x 16 = 48) and its own
# run_metal_complex_stages.sh, pinned to a disjoint core range.
#
# Self-chaining across the 3-day walltime: at start the job queues its
# own successor with --dependency=afterany:<this job>, so a walltime
# kill still continues. The successor is cancelled when every directory
# is validated (.stages/ALL_DONE), when any directory FAILED, or when
# MAX_ROUNDS is reached.
#
# Usage (from anywhere):
#   sbatch job_metal_complex_trio.sh <dir1> [<dir2> <dir3> ...]
# Optional environment (pass with sbatch --export=ALL,VAR=...):
#   MAX_ROUNDS (default 6), ROUND (internal), NUMFORCE_OPTS, JOBEX_CYCLES

set -u

SCRIPT_DIR=${SCRIPT_DIR:-$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")}
# under sbatch BASH_SOURCE points to the spooled copy: prefer the submit dir copy
[ -f "$SCRIPT_DIR/run_metal_complex_stages.sh" ] || SCRIPT_DIR=$HOME/work_turbomole/metal_complexes/scripts/scripts
RUNNER="$SCRIPT_DIR/run_metal_complex_stages.sh"
SELF="$SCRIPT_DIR/job_metal_complex_trio.sh"

ROUND=${ROUND:-1}
MAX_ROUNDS=${MAX_ROUNDS:-6}
CORES_TOTAL=${SLURM_CPUS_PER_TASK:-48}

[ $# -ge 1 ] || { echo "Usage: sbatch $0 <dir1> [<dir2> ...]"; exit 1; }
DIRS=()
for d in "$@"; do DIRS+=("$(readlink -f "$d")"); done

echo "=== $(date) mc_trio round $ROUND/$MAX_ROUNDS, job ${SLURM_JOB_ID:-none} on $(hostname) ==="
echo "dirs: ${DIRS[*]}"

# --- decide whether there is anything left to do -----------------------------
pending=()
for d in "${DIRS[@]}"; do
    if [ -f "$d/.stages/FAILED" ]; then
        echo "STOP: $d has .stages/FAILED: $(cat "$d/.stages/FAILED")"
        echo "Fix it, rm the FAILED marker and resubmit manually."
        exit 1
    fi
    [ -f "$d/.stages/ALL_DONE" ] || pending+=("$d")
done
if [ ${#pending[@]} -eq 0 ]; then
    echo "All directories already validated, nothing to do."
    exit 0
fi

# cores are split among the directories still pending (3 -> 16 each; once
# only one is left it gets the whole node)
CORES_PER_DIR=$((CORES_TOTAL / ${#pending[@]}))
echo "pending: ${#pending[@]} -> $CORES_PER_DIR cores each"

# --- queue successor so a walltime kill does not end the chain ----------------
NEXT_JOB=""
if [ -n "${SLURM_JOB_ID:-}" ] && [ "$ROUND" -lt "$MAX_ROUNDS" ]; then
    NEXT_JOB=$(cd "$SLURM_SUBMIT_DIR" && sbatch --parsable \
        --dependency=afterany:"$SLURM_JOB_ID" \
        --export=ALL,ROUND=$((ROUND + 1)),MAX_ROUNDS="$MAX_ROUNDS",SCRIPT_DIR="$SCRIPT_DIR" \
        "$SELF" "${DIRS[@]}")
    echo "successor queued: $NEXT_JOB (round $((ROUND + 1)))"
else
    echo "WARNING: no successor (round $ROUND = MAX_ROUNDS or not under Slurm)"
fi

# --- launch one runner per directory, pinned to disjoint cores -----------------
pids=()
i=0
for d in "${pending[@]}"; do
    first=$((i * CORES_PER_DIR)); last=$((first + CORES_PER_DIR - 1))
    echo "start: $d on cores $first-$last"
    taskset -c "$first-$last" bash "$RUNNER" "$d" "$CORES_PER_DIR" \
        > "$d/stages_round${ROUND}_${SLURM_JOB_ID:-local}.log" 2>&1 &
    pids+=($!)
    i=$((i + 1))
done

for p in "${pids[@]}"; do wait "$p"; done

# --- summary and chain control ------------------------------------------------
all_done=1; failed=0
for d in "${DIRS[@]}"; do
    if [ -f "$d/.stages/ALL_DONE" ]; then s="ALL_DONE"
    elif [ -f "$d/.stages/FAILED" ]; then s="FAILED"; failed=1; all_done=0
    else s="pending"; all_done=0; fi
    echo "$(date '+%F %T') $s  $d"
done

if [ -n "$NEXT_JOB" ] && { [ $all_done -eq 1 ] || [ $failed -eq 1 ]; }; then
    scancel "$NEXT_JOB" && echo "successor $NEXT_JOB cancelled (all_done=$all_done failed=$failed)"
fi
echo "=== $(date) mc_trio round $ROUND finished ==="
