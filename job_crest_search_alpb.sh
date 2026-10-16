#!/bin/bash
#SBATCH -t 08:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o crest_alpb.out
#SBATCH -e crest_alpb.error
#SBATCH -J crest_alpb

# CREST iMTD-GC / GFN2-xTB with ALPB(water): builds the conformer ensemble IN SOLVENT, because the vacuum
# ensembles were shown to be biased for the localised-charge anions (Fig. 12, protonation_and_conformer_investigation).
# Same settings as job_crest_search.sh (single-threaded BLAS, -T for concurrency) except --alpb water.
#
# Usage: sbatch job_crest_search_alpb.sh <dir with crest_input.xyz (or coord) and CHARGE.txt> [ewin_kcal]
set -e
ulimit -s unlimited
module purge
module load miniconda/24.4.0-libmamba
source /opt/gridunesp/dist/miniconda/24.4.0/libmamba/etc/profile.d/conda.sh
conda activate paper2_funnel
WORK_DIR=${1:?"Usage: sbatch job_crest_search_alpb.sh <dir> [ewin]"}
EWIN=${2:-8}
cd "$WORK_DIR"
CHARGE=$(cat CHARGE.txt)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
echo "CREST/GFN2/ALPB(water) in $WORK_DIR charge=$CHARGE threads=$SLURM_NTASKS ewin=$EWIN"
crest coord --gfn2 --alpb water --chrg "$CHARGE" -T "$SLURM_NTASKS" -ewin "$EWIN" > crest_search.log 2>&1
echo "CREST finished (exit $?) for $WORK_DIR"
