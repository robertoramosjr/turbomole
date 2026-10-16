#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J cosmo_tddft

# COSMO-TDDFT, HSE06/def2-TZVP, eps 78.4, 60 roots on frozen r2SCAN-3c geometries (fast/nonequilibrium
# response through the COSMO refractive index). Dirs prepared by paper_2/scripts/prepare_eps_cosmo_tddft.py;
# ridft (COSMO SCF) then escf. Paper 2 plans 2026-10-16 (E3 of epsilon, H4 of anion band).

set -u
# COSMO-TDDFT (eps 78.4, fast/nonequilibrium response via refind) escf on frozen r2SCAN-3c geometries, HSE06/def2-TZVP, 60 roots.
# Dirs prepared by paper_2/scripts/prepare_eps_cosmo_tddft.py; ridft (COSMO SCF) then escf. Plans 2026-10-16 (E3/H4).
ulimit -s unlimited

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

JOBLIST=${1:?"Usage: sbatch --array=1-N job_cosmo_tddft_escf.sh <joblist_file>"}
WORK_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")

if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: no valid directory on line $SLURM_ARRAY_TASK_ID of $JOBLIST"
    exit 1
fi

cd "$WORK_DIR"
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS
export OMP_NUM_THREADS=$SLURM_NTASKS
export OMP_STACKSIZE=4G
export OMP_PLACES=cores
export OMP_PROC_BIND=close

ridft > ridft_cosmo.out 2>&1
ridft_status=$?
echo "ridft exit $ridft_status for $WORK_DIR"
grep -qi "ended abnormally\|did not converge" ridft_cosmo.out && ridft_status=1
if [ $ridft_status -ne 0 ]; then
    exit $ridft_status
fi

escf > escf.out 2>&1
escf_status=$?
grep -qi "ended abnormally" escf.out && escf_status=1
echo "escf exit $escf_status for $WORK_DIR"

exit $escf_status
