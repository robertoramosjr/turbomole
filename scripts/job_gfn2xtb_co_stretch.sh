#!/bin/bash
#SBATCH -t 05:00:00
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH -n 8
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH -o slurm_gfn2xtb_costretch.out
#SBATCH -e slurm_gfn2xtb_costretch.error
#SBATCH -J gfn2xtb_costretch

# Runs batch_gfn2xtb_co_stretch.py for all 4 protonation states'
# labeled subsets -- the cheap-proxy half of B1's delta-learning data
# (nu(C=O) at GFN2-xTB, to be matched against the HSE06/def2-TZVP
# value from job_hse06_labeling_retry2.sh once that finishes). Moved
# to the SLURM queue instead of running on the login node (as it had
# been, briefly, by mistake earlier in this project -- same lesson as
# Sec. 13 of protonation_and_conformer_investigation.md).

set -u
ulimit -s unlimited

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

export OMP_NUM_THREADS=$SLURM_NTASKS
export OMP_STACKSIZE=2G

cd ~/work_turbomole/artepillin_C_studies

declare -A charges=( [neutral]=0 [monoanion_a]=-1 [monoanion_b]=-1 [dianion]=-2 )
for state in neutral monoanion_a monoanion_b dianion; do
    d="artepillin_c_${state}_d3_hse06"
    chg=${charges[$state]}
    echo "=== $state (charge $chg) ==="
    python3 ~/work_turbomole/scripts/batch_gfn2xtb_co_stretch.py \
        --subset "$d/layer_b/initial_labeled_subset.csv" \
        --r2scan3c-jobs "$d/layer_a/conformer_jobs" \
        --charge "$chg" \
        --output "$d/layer_b/gfn2xtb_co_stretch.csv"
    echo
done
echo "ALL STATES DONE"
