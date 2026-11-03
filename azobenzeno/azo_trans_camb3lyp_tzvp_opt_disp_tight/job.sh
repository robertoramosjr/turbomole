#!/bin/bash
#SBATCH --job-name=azo_trans_tight
#SBATCH --time=1-00:00:00
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

# Test 2 (imaginary-mode diagnosis): same CAM-B3LYP-D3(BJ)/def2-TZVP + $senex, TM 8.0.1, but
# standard grid 4 (manual l.3888: multiple grids not recommended for hybrids/response),
# $scfconv 10, $denconv 1d-10, $forceconv 7. jobex re-relaxes on the fine grid, then aoforce only.
# Each stage validated by TEXT; markers in .stages/. Re-run sbatch to resume.
# Usage: sbatch job.sh   (from this directory)
set -e
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE801   # TM 8.0.1 (NOT ~/software/TURBOMOLE = 7.8.1)
export PATH=$TURBODIR/scripts:$PATH
export PARA_ARCH=SMP
export PARNODES=$SLURM_CPUS_PER_TASK
source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/bin/$(sysname):$PATH
which ridft aoforce jobex | grep -c TURBOMOLE801 | grep -q 3 || { echo "wrong TURBODIR"; exit 1; }
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd $SLURM_SUBMIT_DIR
mkdir -p .stages

stage() {  # stage <name> <ok-pattern> <logfile> <cmd...>
  local n=$1 pat=$2 log=$3; shift 3
  [ -f .stages/$n ] && { echo "skip $n"; return; }
  local t0=$(date +%s)
  "$@" > $log 2>> job.err || true
  grep -q "$pat" $log || { echo "FAIL $n (no '$pat' in $log)"; exit 1; }
  grep -q "ended abnormally\|did not converge" $log && { echo "FAIL $n (abnormal)"; exit 1; }
  echo "$n ok $(($(date +%s)-t0)) s" | tee .stages/$n
}

[ -f .stages/opt ] || {
  rm -f GEO_OPT_FAILED not.converged
  jobex -c 800 > jobex.out 2>> job.err || true
  [ -f GEO_OPT_CONVERGED ] || { echo "FAIL opt (GEO_OPT_FAILED)"; exit 1; }
  echo "opt ok" > .stages/opt
}
stage ridft "ridft : all done" ridft_tight.out ridft
stage aoforce "force : all done" aoforce.out aoforce
if grep -q '\*\*\*\*' vibspectrum; then echo "FAIL aoforce: overflow in vibspectrum"; exit 1; fi
awk '$1 ~ /^[0-9]+$/ && $3+0 != 0 {print $3}' vibspectrum | sort -g | head -3 > lowest_modes.txt
if grep -q '^-' lowest_modes.txt; then echo "WARN imaginary mode: $(head -1 lowest_modes.txt)"; else echo "no imaginary mode"; fi
echo ALL_OK
