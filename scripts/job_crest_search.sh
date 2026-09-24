#!/bin/bash
#SBATCH -t 04:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --mail-user=robeerto.aguiar.ramos@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o crest_search.out
#SBATCH -e crest_search.error
#SBATCH -J crest_search

# Runs a CREST iMTD-GC / GFN2-xTB conformational search in a target
# directory (a protonation-state seed folder containing 'coord' and
# 'CHARGE.txt', as produced by build_protonation_state.py). One state
# per submission -- submit once per protonation-state folder.
#
# Usage:
#   sbatch ~/work_turbomole/scripts/job_crest_search.sh /path/to/artepillin_c_<state>_d3_hse06 [ewin_kcal]
#
# Walltime is set generously (4h) relative to the ~15-60 min estimate
# from a direct single-thread MD throughput benchmark on this cluster
# (17 s/ps GFN2-xTB, 45-46 atoms) -- CREST may need a second MTD
# iteration if new low-lying minima are found, which the benchmark
# alone does not account for.

set -e
ulimit -s unlimited

module purge
module load miniconda/24.4.0-libmamba

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

WORK_DIR=${1:?"Usage: sbatch job_crest_search.sh <state_dir> [ewin_kcal]"}
EWIN=${2:-6}

cd "$WORK_DIR"

if [ ! -f coord ] || [ ! -f CHARGE.txt ]; then
    echo "ERROR: '$WORK_DIR' must contain 'coord' and 'CHARGE.txt' "
    echo "       (run build_protonation_state.py first)."
    exit 1
fi

CHARGE=$(cat CHARGE.txt)

# CREST's own "-T" flag parallelizes across concurrent MTDs; each
# individual MTD's GFN2-xTB/tblite calculation must stay single-threaded
# or OpenBLAS tries to multithread *within* every already-parallel MTD,
# causing massive core oversubscription (observed: constant "OpenBLAS
# Warning: Detect OpenMP Loop" spam bloating logs to 500MB-1GB+, and a
# large slowdown vs. a clean single-thread benchmark on the same
# hardware). Keep the BLAS layer single-threaded and let "-T" alone
# manage concurrency across MTDs.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

echo "Starting CREST iMTD-GC/GFN2-xTB search in $WORK_DIR"
echo "  charge=$CHARGE  threads=$SLURM_NTASKS  ewin=$EWIN kcal/mol"

crest coord --gfn2 --chrg "$CHARGE" -T "$SLURM_NTASKS" -ewin "$EWIN" \
    > crest_search.log 2>&1

echo "CREST finished (exit code $?) for $WORK_DIR"
