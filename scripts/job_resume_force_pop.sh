#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --partition=medium
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J resume_force_pop

# Resumes ONLY the aoforce + dscf -proper tail of the metal-complex
# pipeline, for directories where jobex/ridft-proper/bader/ridft_tight/
# escf already succeeded but aoforce got walltime-killed mid-Hessian
# (job_metal_complex_pipeline.sh's own 3-day cap). Avoids redundantly
# re-running the already-successful earlier stages.

set -u
ulimit -s unlimited

export TURBODIR=$HOME/software/TURBOMOLE
export PATH=$TURBODIR/scripts:$PATH
source "$TURBODIR/Config_turbo_env"
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS

JOBLIST=${1:?"Usage: sbatch --array=1-N job_resume_force_pop.sh <joblist_file>"}
WORK_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")

if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: no valid directory on line $SLURM_ARRAY_TASK_ID of $JOBLIST"
    exit 1
fi

cd "$WORK_DIR"

TIMING_LOG="$WORK_DIR/timing.log"
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

echo "=== $(date) resuming aoforce/dscf-proper in $WORK_DIR (PARNODES=$PARNODES) ==="

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

echo "=== $(date) finished in $WORK_DIR ==="
column -s, -t "$TIMING_LOG"
