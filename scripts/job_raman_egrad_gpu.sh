#!/bin/bash
#SBATCH -J tm_raman_eg
#SBATCH --partition=gpu-x
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH -t 12:00:00
#SBATCH -o slurm_raman_egrad_%j.out

# Raman on GPU, TURBOMOLE analytic route: egrad ($scfinstab polly) -> intense.
# Needs a COPY of a finished aoforce dir (hessian, vib_normal_modes, vibspectrum, mos...).
# Usage: sbatch job_raman_egrad_gpu.sh <dir>
set -u
ulimit -s unlimited
W=$(readlink -f "${1:?dir}"); cd "$W" || exit 1
export PARA_ARCH=SMP PARNODES=${SLURM_CPUS_PER_TASK:-16} OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
source "$HOME/software/turbomole_gpu_env.sh"
export TURBOTMPDIR="${TMPDIR:-/tmp}/tm_ramg_${SLURM_JOB_ID:-local}_$$"; mkdir -p "$TURBOTMPDIR"
trap 'rm -rf "$TURBOTMPDIR"' EXIT
grep -q '^\$scfinstab polly' control || sed -i '/^\$end/i \$scfinstab polly' control
actual -r > /dev/null 2>&1; sed -i '/^\$restart/d' control
t0=$(date +%s); egrad > egrad.out 2>&1; t1=$(date +%s)
echo "egrad $((t1-t0)) s, all done: $(grep -c 'all done' egrad.out)"
grep -q 'all done' egrad.out || { echo "egrad FAILED"; exit 1; }
intense > intense.out 2>&1
echo "intense: $(grep -c 'ended normally' intense.out) ended normally"
sed -n '/^\$raman spectrum/,/^\$/p' control > raman_spectrum.dat
