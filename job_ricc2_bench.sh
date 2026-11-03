#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --partition=medium
#SBATCH --mem=64G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J ricc2_bench

# COSMO-ADC(2) and COSMO-CC2 (post-SCF reaction field, LR) vertical excitations on top of an HF/COSMO reference,
# plus the same-basis HSE06 COSMO-TD-DFT in the sibling dir.  Dirs from prepare_ricc2_bench.py (hf/ and dft/).
# Stages are resumable by markers: HF_DONE, ADC2_DONE, CC2_DONE, ESCF_DONE; each validated by TEXT.
# Usage: sbatch --array=1-N%4 job_ricc2_bench.sh <joblist_file>   (joblist = absolute .../hf dirs)
set -u
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE; source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
JOBLIST=${1:?"Usage: sbatch --array=1-N job_ricc2_bench.sh <joblist_file>"}
HF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
[ -n "$HF" ] && [ -d "$HF" ] || { echo "ERROR: no valid dir on line $SLURM_ARRAY_TASK_ID of $JOBLIST"; exit 1; }
BASE=$(dirname "$HF"); DFT=$BASE/dft
export PARA_ARCH=SMP PARNODES=$SLURM_NTASKS OMP_NUM_THREADS=$SLURM_NTASKS OMP_STACKSIZE=4G OMP_PLACES=cores OMP_PROC_BIND=close

# 1) same-basis HSE06 COSMO-TD-DFT
cd "$DFT"
if [ ! -f ESCF_DONE ]; then
  [ -f RIDFT_DONE ] || { ridft > ridft.out 2>&1; grep -q "ridft : all done" ridft.out && ! grep -qi "ended abnormally\|did not converge" ridft.out && echo RIDFT_DONE > RIDFT_DONE; }
  [ -f RIDFT_DONE ] && { escf > escf.out 2>&1; grep -q "escf : all done" escf.out && ! grep -qi "ended abnormally" escf.out && echo ESCF_DONE > ESCF_DONE; }
fi

# 2) HF/COSMO reference
cd "$HF"
if [ ! -f HF_DONE ]; then
  ridft > ridft.out 2>&1
  if grep -q "ridft : all done" ridft.out && ! grep -qi "ended abnormally\|did not converge" ridft.out; then echo HF_DONE > HF_DONE
  else echo "FAIL HF ridft in $HF"; exit 1; fi
fi

# 3) ricc2 for each model on a private copy of the converged reference (ricc2 rewrites files and cosmo charges)
for MODEL in "adc(2)" "cc2"; do
  TAG=$(echo "$MODEL" | tr -d '()' | tr 'a-z' 'A-Z')
  W=$BASE/ricc2_$TAG
  [ -f "$W/${TAG}_DONE" ] && continue
  rm -rf "$W"; mkdir -p "$W"; cp "$HF"/{control,coord,basis,mos,auxbasis,out.ccf} "$HF"/*.tmp "$W"/ 2>/dev/null
  cd "$W"
  sed -i "s/^\$freeze/\$ricc2\n $MODEL\n maxiter=100\n\$freeze/" control
  ricc2 > ricc2.out 2>&1
  if grep -qi "ricc2 : all done" ricc2.out && ! grep -qi "ended abnormally" ricc2.out; then echo done > "${TAG}_DONE"; else echo "FAIL ricc2 $MODEL in $W"; fi
done
echo "ricc2 bench $BASE: dft=$([ -f $DFT/ESCF_DONE ] && echo 1 || echo 0) adc2=$([ -f $BASE/ricc2_ADC2/ADC2_DONE ] && echo 1 || echo 0) cc2=$([ -f $BASE/ricc2_CC2/CC2_DONE ] && echo 1 || echo 0)"
