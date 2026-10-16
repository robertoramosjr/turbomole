#!/bin/bash
#SBATCH -t 03:00:00
#SBATCH -N 1 --ntasks-per-node=4
#SBATCH -n 4
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J pka_refcmpd

# Full pKa reference-compound pipeline for ONE species directory built
# by setup_hse06_opt_template.py (protonation_and_conformer_investigation.md
# Sec. 9.34): jobex (cycle cap 800 -- user directive for optimization
# stages) -> aoforce (Hessian) -> freeh (thermal G correction, run on a
# COPY via run_freeh_thermal_correction.py) -> COSMO(water) setup
# (setup_cosmo_job.py) -> COSMO single-point ridft. Each stage's text
# output is grepped ("ended abnormally"), never just the exit code.
# Small molecules (<=15 atoms): measured for acetic acid ~0.15 core-h
# for jobex, ~0.07 for aoforce, ~0.01 for the COSMO single point.
#
# Usage: sbatch --array=1-N job_reference_compound_pipeline.sh <joblist>
#   (joblist = absolute species dirs, one per line, on the SHARED
#   filesystem -- never /tmp, compute nodes cannot see it)

set -u
ulimit -s unlimited

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

JOBLIST=${1:?"Usage: sbatch --array=1-N job_reference_compound_pipeline.sh <joblist>"}
WORK_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: no valid directory on line $SLURM_ARRAY_TASK_ID of $JOBLIST"
    exit 1
fi

SCRIPTS=/home/rramos/work_turbomole/scripts
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS
export OMP_NUM_THREADS=$SLURM_NTASKS
export OMP_STACKSIZE=2G

cd "$WORK_DIR"

[ -f GEO_OPT_CONVERGED ] || jobex -c 800 > jobex.out 2>&1   # resumable: skip jobex if already converged
if [ ! -f GEO_OPT_CONVERGED ]; then
    echo "FAIL: no GEO_OPT_CONVERGED in $WORK_DIR"; exit 1
fi
echo "jobex converged: $(grep -c 'OPTIMIZATION CYCLE' jobex.out) cycles"

aoforce > aoforce.out 2>&1
if ! grep -q "aoforce ended normally" aoforce.out || grep -qi "ended abnormally" aoforce.out; then
    echo "FAIL: aoforce in $WORK_DIR"; exit 1
fi
# aoforce can print "ended normally" with a numerically corrupt Hessian (seen 2026-10-15, ferulate acid:
# elements ~5e2 vs <1 normally, frequencies printed as '**********', freeh then unparsable). Catch it here.
if grep -q '\*\*\*\*\*\*\*\*\*\*' vibspectrum hessian; then
    echo "FAIL: aoforce produced overflowed ('*****') Hessian/frequencies in $WORK_DIR -- rerun aoforce from scratch"; exit 1
fi
NIMAG=$(awk '$1 ~ /^[0-9]+$/ && $3+0 < -10 {c++} END{print c+0}' vibspectrum)
echo "imaginary modes (< -10 cm-1): $NIMAG"

# freeh on a copy, via the wrapper (needs parent-dir/job-name layout)
PARENT=$(dirname "$WORK_DIR"); NAME=$(basename "$WORK_DIR")
echo "$NAME" > freeh_joblist.txt
python3 "$SCRIPTS/run_freeh_thermal_correction.py" \
    --job-list freeh_joblist.txt --labeled-jobs "$PARENT" \
    --output freeh_thermal_correction.csv > freeh.log 2>&1
[ -s freeh_thermal_correction.csv ] && [ "$(wc -l < freeh_thermal_correction.csv)" -ge 2 ] \
    || { echo "FAIL: freeh in $WORK_DIR"; cat freeh.log; exit 1; }

COSMO_DIR="${WORK_DIR}_cosmo"
rm -rf "$COSMO_DIR"
python3 "$SCRIPTS/setup_cosmo_job.py" --source-dir "$WORK_DIR" \
    --output-dir "$COSMO_DIR" --epsilon 78.4 > "$WORK_DIR/cosmo_setup.log" 2>&1
if ! grep -q '^\$cosmo$' "$COSMO_DIR/control"; then
    echo "FAIL: no \$cosmo group written in $COSMO_DIR (silent cosmoprep failure)"; exit 1
fi

cd "$COSMO_DIR"
ridft > ridft_cosmo.out 2>&1
if ! grep -q "ridft ended normally" ridft_cosmo.out || grep -qi "ended abnormally\|did not converge" ridft_cosmo.out; then
    echo "FAIL: COSMO ridft in $COSMO_DIR"; exit 1
fi
grep "Total energy + OC" ridft_cosmo.out
echo "PIPELINE OK: $WORK_DIR"
