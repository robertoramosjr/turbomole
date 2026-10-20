#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH -n 8
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH -o slurm_raman_%A_%a.out
#SBATCH -e slurm_raman_%A_%a.error
#SBATCH -J raman_polly

# Raman, analytic TURBOMOLE route on a COPY of a finished PBE0/def2-SVP aoforce dir (prepare_raman_jobs.py):
# egrad ($scfinstab polly) -> intense -> raman_spectrum.dat. Resumable: a dir with raman_spectrum.dat is skipped.
# Validate by TEXT: egrad.out "all done", intense.out "ended normally", non-empty raman_spectrum.dat.
# If egrad/intense complain about $rik, report the exact error; do not change the exchange scheme here.
# Usage: sbatch --array=1-N%T job_raman_polly.sh <joblist.txt>   (absolute paths on the shared FS)
set -u
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE; source $TURBODIR/Config_turbo_env
export PARA_ARCH=SMP PARNODES=$SLURM_NTASKS OMP_NUM_THREADS=$SLURM_NTASKS
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH

W=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${1:?joblist}")
[ -d "$W" ] || { echo "ERROR: no dir on line $SLURM_ARRAY_TASK_ID of $1"; exit 1; }
cd "$W" || exit 1
[ -s raman_spectrum.dat ] && { echo "already done $W"; exit 0; }
grep -q '^\$scfinstab polly' control || sed -i '/^\$end/i \$scfinstab polly' control
actual -r > /dev/null 2>&1; sed -i '/^\$restart/d' control
t0=$(date +%s); egrad > egrad.out 2>&1; t1=$(date +%s)
echo "egrad $((t1 - t0)) s, all done: $(grep -c 'all done' egrad.out)"
grep -q 'all done' egrad.out || { echo "egrad FAILED in $W"; exit 1; }
intense > intense.out 2>&1
echo "intense: $(grep -c 'ended normally' intense.out) ended normally (stderr banner may be in .error)"
sed -n '/^\$raman spectrum/,/^\$/p' control > raman_spectrum.dat
[ -s raman_spectrum.dat ] && echo "raman_spectrum.dat $(wc -l < raman_spectrum.dat) lines" || { echo "intense produced no \$raman spectrum in $W"; exit 1; }
