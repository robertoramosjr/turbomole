#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH -n 1
#SBATCH --partition=short
#SBATCH --mem=4G
#SBATCH -o slurm_stda_confidence_neutral.out
#SBATCH -e slurm_stda_confidence_neutral.error
#SBATCH -J stda_confidence_neutral

# Runs characterize_stda_selection_confidence.py on the full 170-
# conformer neutral Layer A ensemble (the only state with unlabeled
# conformers left to score for B2's oscillator-strength selection
# confidence -- monoanion_b's 11-conformer ensemble is already 100%
# labeled, monoanion_a/dianion do not use this max-f criterion at all,
# see the script's own docstring). Single-threaded subprocess work
# (tm2molden + std2 per conformer, ~4s each, ~12 min total measured on
# a 5-conformer sample) -- submitted via SLURM rather than run on the
# login node per this project's standing rule.

set -u
ulimit -s unlimited

CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

cd ~/work_turbomole/artepillin_C_studies

d="artepillin_c_neutral_d3_hse06"
ls "$d/layer_a/conformer_jobs" > /tmp/joblist_neutral_full.txt

python3 ~/work_turbomole/scripts/characterize_stda_selection_confidence.py \
    --job-list /tmp/joblist_neutral_full.txt \
    --r2scan3c-jobs "$d/layer_a/conformer_jobs" \
    --calibration-csv "$d/layer_b/stda_state_selection_check_win10.csv" \
    --ax 0.25 --window-states 10 \
    --output "$d/layer_b/stda_confidence_full_ensemble.csv"

echo "DONE"
