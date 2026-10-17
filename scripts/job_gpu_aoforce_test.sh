#!/bin/bash
#SBATCH -J tm_gpu_test
#SBATCH --partition=gpu-x
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH -o slurm_gpu_test_%j.out

# TURBOMOLE 8.0.1 GPU test on coaraci gpu-x.
#  Part 1 (correctness): <small_dir> ridft + aoforce on CPU, GPU fp32, GPU fp64;
#          compare vibspectrum afterwards.
#  Part 2 (performance): aoforce on GPU in <big_dir> for the rest of the
#          walltime, logging progress every PROGRESS_EVERY seconds so the
#          total time can be extrapolated.
# Usage: sbatch -t <time> job_gpu_aoforce_test.sh <small_dir> <big_dir>

set -u
ulimit -s unlimited
SMALL=$(readlink -f "${1:?small_dir}"); BIG=$(readlink -f "${2:?big_dir}")
PROGRESS_EVERY=${PROGRESS_EVERY:-300}

export PARA_ARCH=SMP PARNODES=${SLURM_CPUS_PER_TASK:-16} OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
source "$HOME/software/turbomole_gpu_env.sh"
export TURBOTMPDIR="${TMPDIR:-/tmp}/tm_gpu_${SLURM_JOB_ID:-local}_$$"; mkdir -p "$TURBOTMPDIR"
trap 'rm -rf "$TURBOTMPDIR"' EXIT

echo "=== $(date) node $(hostname) ==="
nvidia-smi --query-gpu=index,name,driver_version,memory.total --format=csv
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"

run_case() {   # <label> <gpu: on|off> <precision: fp32|fp64|->
    local label=$1 gpu=$2 prec=$3 W="$SMALL/$1" t0 t1 t2
    rm -rf "$W"; mkdir -p "$W"; cp -p "$SMALL"/{control,coord,basis,auxbasis,mos} "$W"/; cd "$W"
    if [ "$gpu" = on ]; then
        export TM_GPGPU=on
        [ "$prec" = fp64 ] && sed -i '/^\$end/i \$gpgpu fp64' control
    else
        unset TM_GPGPU
    fi
    t0=$(date +%s); ridft > ridft.out 2>&1; t1=$(date +%s)
    aoforce > aoforce.out 2>&1; t2=$(date +%s)
    printf "%-6s ridft %5d s (%s)  aoforce %5d s (%s)  gpu_lines=%s\n" "$label" $((t1-t0)) \
        "$(grep -c 'ended normally' ridft.out)" $((t2-t1)) "$(grep -c 'ended normally' aoforce.out)" \
        "$(grep -c 'GPU initialized' ridft.out aoforce.out | awk -F: '{s+=$2} END{print s}')"
    grep "total energy" ridft.out | tail -1
}

echo "--- Part 1: correctness on $SMALL ---"
run_case cpu  off -
run_case gpu32 on fp32
run_case gpu64 on fp64

echo "--- Part 2: aoforce GPU (fp32) on $BIG ---"
cd "$BIG"
export TM_GPGPU=on
actual -r > /dev/null 2>&1
sed -i '/^\$restart/d' control
t0=$(date +%s)
aoforce > aoforce.out 2>&1 &
pid=$!
while kill -0 $pid 2>/dev/null; do
    sleep "$PROGRESS_EVERY"
    echo "$(date '+%F %T') +$(( $(date +%s) - t0 )) s | $(grep -E 'CPKS|iteration|Hessian|ended' aoforce.out | tail -1 | tr -s ' ')" >> progress.log
done
wait $pid
echo "aoforce finished after $(( $(date +%s) - t0 )) s: $(grep -c 'ended normally' aoforce.out) 'ended normally'"
