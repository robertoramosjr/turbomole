#!/bin/bash
#SBATCH -J w2_template
#SBATCH -t 00:40:00
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --mem=4G
#SBATCH -o slurm_w2_template_%j.out
# Builds the four COSMO-ready r2SCAN-3c templates (neutral, monoanion_a from layer_a conformers via cosmoprep;
# monoanion_b, dianion reuse the cosmo-prepped solv_relax_pilot dirs).
source /opt/gridunesp/dist/miniconda/24.4.0/libmamba/etc/profile.d/conda.sh; conda activate paper2_funnel
export TURBODIR=$HOME/software/TURBOMOLE; source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
R=/home/rramos/work_turbomole/artepillin_C_studies; S=/home/rramos/work_turbomole/scripts
cd $R && mkdir -p tier2
python3 $S/water_tier2_sp.py template --source artepillin_c_neutral_d3_hse06/layer_a/conformer_jobs/job_0044 --out tier2/neutral/template
python3 $S/water_tier2_sp.py template --source artepillin_c_monoanion_a_d3_hse06/layer_a/conformer_jobs/job_0093 --out tier2/monoanion_a/template
mkdir -p tier2/monoanion_b tier2/dianion
rm -rf tier2/monoanion_b/template tier2/dianion/template
cp -r solv_relax_pilot/B_old_job_0001 tier2/monoanion_b/template
cp -r solv_relax_pilot/D_old_job_0001 tier2/dianion/template
