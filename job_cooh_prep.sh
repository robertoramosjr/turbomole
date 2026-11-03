#!/bin/bash
#SBATCH -t 00:30:00
#SBATCH -N 1 -n 1
#SBATCH --partition=short
#SBATCH --mem=4G
#SBATCH -J cooh_prep
#SBATCH -o /home/rramos/work_turbomole/artepillin_C_studies/todo_v12/refs_cooh/slurm_prep_%j.out
# Light preparation step of the COOH-reference chain: prepare_ref_water.py <reopt|w4> restricted to the new references.
# define/cosmoprep need the TURBOMOLE environment (same as job_ref_reopt.sh).  Usage: sbatch job_cooh_prep.sh <reopt|w4>
set -eu
STEP=${1:?"Usage: sbatch job_cooh_prep.sh <reopt|w4>"}
source /opt/gridunesp/dist/miniconda/24.4.0/libmamba/etc/profile.d/conda.sh; conda activate paper2_funnel
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE; source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
python /home/rramos/work_turbomole/scripts/prepare_ref_water.py "$STEP" --root /home/rramos/work_turbomole/artepillin_C_studies --only cinnamic_acid 4_coumaric_cooh
