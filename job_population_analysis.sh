#!/bin/bash
#SBATCH -t 00:30:00
#SBATCH -N 1 --ntasks-per-node=4
#SBATCH -n 4
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J pop_analysis

# Population-analysis step for B1 feature expansion: bond order/SEN
# (dscf -proper, no $pointval -- fast) then a Bader density cube
# (ridft -proper, $pointval injected only after dscf finishes -- see
# setup_population_job.py's docstring for why the two must not share a
# control with $pointval already present) then the external `bader`
# binary. Reuses the already-converged `mos` from the HSE06 labeling
# jobs -- no new SCF-from-scratch, no new Hessian.
#
# Real per-conformer timing (measured in a 1-conformer trial): dscf
# ~16s, ridft cube ~3min (does not parallelize well -- keep core count
# modest), bader ~25s. 4 cores is enough; more would not meaningfully
# speed up the cube step.
#
# Usage:
#   sbatch --array=1-N%K job_population_analysis.sh <joblist_file>

set -u
ulimit -s unlimited

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

JOBLIST=${1:?"Usage: sbatch --array=1-N job_population_analysis.sh <joblist_file>"}
WORK_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")

if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: no valid directory on line $SLURM_ARRAY_TASK_ID of $JOBLIST"
    exit 1
fi

cd "$WORK_DIR"
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS
export OMP_NUM_THREADS=$SLURM_NTASKS

dscf -proper < /dev/null > pop.out 2>&1
dscf_status=$?
grep -qi "ended abnormally" pop.out && dscf_status=1
echo "dscf -proper exit $dscf_status for $WORK_DIR"
if [ $dscf_status -ne 0 ]; then
    exit $dscf_status
fi

python3 -c "
with open('control') as f:
    text = f.read()
with open('pointval_block.txt') as f:
    block = f.read()
text = text.replace('\$end', block + '\$end', 1)
with open('control', 'w') as f:
    f.write(text)
"

ridft -proper < /dev/null > proper.out 2>&1
ridft_status=$?
grep -qi "ended abnormally" proper.out && ridft_status=1
echo "ridft -proper exit $ridft_status for $WORK_DIR"
if [ $ridft_status -ne 0 ] || [ ! -f td.cub ]; then
    exit 1
fi

BADER_BIN=$(command -v bader || echo ~/local_bin/bader)
"$BADER_BIN" td.cub > bader_run.out 2>&1
bader_status=$?
echo "bader exit $bader_status for $WORK_DIR"
if [ ! -f ACF.dat ]; then
    bader_status=1
fi
exit $bader_status
