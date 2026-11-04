#!/bin/bash
#SBATCH --job-name=lig
#SBATCH --time=2-00:00:00
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --partition=par48a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=36000M
# CAM-B3LYP-D3(BJ)/def2-TZVP + $senex, TM 8.0.1 (harpia). Passo 2 do plano 2026-11-03:
# jobex -ri -c 800 -> ridft -> rdgrad -> NumForce (Hessiana de referência; aoforce dá -23 cm-1 espúrio no trans,
# ver cofre _log/azobenzeno 2026-11-04) -> freeh (escala 1.0, sigma default). Cada estágio validado por TEXTO;
# marcadores em .stages/ (re-sbatch retoma). Uso: sbatch job.sh na pasta do ligante (control/coord/basis/mos).
set -e
ulimit -s unlimited
module load modmat/turbomole/8.0.1
export PARA_ARCH=SMP PARNODES=$SLURM_CPUS_PER_TASK OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_PLACES=cores OMP_PROC_BIND=close
cd $SLURM_SUBMIT_DIR
mkdir -p .stages
grep -q '^$scfiterlimit 1500' control || { echo "FAIL: scfiterlimit != 1500"; exit 1; }
stage() {  # stage <nome> <padrão-ok> <log> <cmd...>
  local n=$1 pat=$2 log=$3; shift 3
  [ -f .stages/$n ] && { echo "skip $n"; return; }
  local t0=$(date +%s)
  "$@" > $log 2>> job.err || true
  grep -q "$pat" $log || { echo "FAIL $n (sem '$pat' em $log)"; exit 1; }
  grep -q "ended abnormally\|did not converge\|cannot be satisfied" $log && { echo "FAIL $n (anormal)"; exit 1; }
  echo "$n ok $(($(date +%s)-t0)) s" | tee .stages/$n
}
[ -f .stages/opt ] || {
  actual -r
  rm -f GEO_OPT_FAILED not.converged
  jobex -ri -c 800 > jobex.out 2>> job.err || true
  [ -f GEO_OPT_CONVERGED ] || { echo "FAIL opt (GEO_OPT_FAILED)"; exit 1; }
  echo "opt ok" > .stages/opt
}
stage ridft "ridft : all done" ridft_final.out ridft
stage rdgrad "rdgrad : all done" rdgrad.out rdgrad
[ -f .stages/numforce ] || {
  NumForce -ri -central -d 0.02 -thrgrd 5.0d-4 > numforce.out 2>> job.err || true
  [ -f vibspectrum ] && grep -q "END OF NumForce" numforce.out || { echo FAIL numforce; exit 1; }
  echo ok > .stages/numforce
}
if grep -q '\*\*\*\*' vibspectrum; then echo "FAIL overflow vibspectrum"; exit 1; fi
awk '$1 ~ /^[0-9]+$/ && $3+0 != 0 {print $3}' vibspectrum | sort -g | head -3 > lowest_modes.txt
if grep -q '^-' lowest_modes.txt; then echo "WARN imaginary mode: $(head -1 lowest_modes.txt)"; else echo "no imaginary mode"; fi
[ -f .stages/freeh ] || {
  mkdir -p freeh && cp control coord hessian vibspectrum vib_normal_modes freeh/
  (cd freeh && printf '\n1.0\n\nq\n' | freeh > freeh.out 2> freeh.err)
  grep -qi "zpe=" freeh/freeh.out || { echo FAIL freeh; exit 1; }
  echo ok > .stages/freeh
}
echo ALL_OK
