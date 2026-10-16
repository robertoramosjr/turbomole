#!/bin/bash
#SBATCH -t 01:30:00
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH -n 8
#SBATCH --partition=short
#SBATCH --mem=24G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J relax_rdgrad

# One HSE06/def2-TZVP gradient (rdgrad) at the FROZEN r2SCAN-3c geometry
# of an already-converged labeled conformer, on a COPY (production
# labeled_jobs untouched). Together with the aoforce Hessian already in
# that job, estimate_relaxation_energy.py turns it into the Newton-step
# relaxation energy dE = -1/2 g^T H^+ g (Sec. 9.36) -- a cheap test of
# whether not re-optimizing the artepilin C conformers at HSE06 biases
# the pKa cycle. joblist lines: <source labeled dir>|<scratch copy dir>
#
# Usage: sbatch --array=1-N job_rdgrad_relax_test.sh <joblist>

set -u
ulimit -s unlimited
CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${1:?Usage: sbatch --array=1-N job_rdgrad_relax_test.sh <joblist>}")
SRC=${LINE%%|*}; DST=${LINE##*|}
[ -d "$SRC" ] || { echo "ERROR: no source dir $SRC"; exit 1; }
rm -rf "$DST"; mkdir -p "$DST"
for f in control coord basis auxbasis mos hessian; do cp "$SRC/$f" "$DST/"; done
cd "$DST"
sed -i '/^\$last step/d;/^\$actual step/d' control

export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS
export OMP_NUM_THREADS=$SLURM_NTASKS
export OMP_STACKSIZE=4G

rdgrad > rdgrad.out 2>&1
if ! grep -q "rdgrad ended normally" rdgrad.out || grep -qi "ended abnormally" rdgrad.out || [ ! -s gradient ]; then
    echo "FAIL: rdgrad in $DST"; tail -n 5 rdgrad.out; exit 1
fi
grep -m1 "cpu-time" rdgrad.out
echo "GRADIENT OK: $DST"
