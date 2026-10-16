#!/bin/bash
#SBATCH -J w2_sp
#SBATCH -t 04:00:00
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH -o slurm_w2_sp_%A_%a.out
# Tier-2 single points; array index = chunk (1..NCHUNKS). Usage: sbatch --array=1-40%12 job_water_tier2_run.sh <state> <NCHUNKS>
source /opt/gridunesp/dist/miniconda/24.4.0/libmamba/etc/profile.d/conda.sh; conda activate paper2_funnel
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE; source $TURBODIR/Config_turbo_env
export PARA_ARCH=SMP PARNODES=8 OMP_NUM_THREADS=8
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
ST=${1:?state}; NCH=${2:?nchunks}
R=/home/rramos/work_turbomole/artepillin_C_studies
cd $R && python3 /home/rramos/work_turbomole/scripts/water_tier2_sp.py run --xyz water_ensemble/$ST/crest_conformers.xyz \
   --template tier2/$ST/template --outdir tier2/$ST --chunk $SLURM_ARRAY_TASK_ID --nchunks $NCH
