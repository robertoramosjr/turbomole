#!/bin/bash
#SBATCH -J ferul_prep
#SBATCH -t 00:40:00
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --mem=4G
#SBATCH -o slurm_ferul_prep_%j.out
source /opt/gridunesp/dist/miniconda/24.4.0/libmamba/etc/profile.d/conda.sh; conda activate paper2_funnel
export TURBODIR=$HOME/software/TURBOMOLE; source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
export OMP_NUM_THREADS=1
S=/home/rramos/work_turbomole/scripts
R=/home/rramos/work_turbomole/artepillin_C_studies/pka_reference
mkdir -p $R/ferulate && cd $R/ferulate
python3 $S/build_ferulate_reference.py --coum-mono ../4_coumarate/coum_mono.xyz --coum-di ../4_coumarate/coum_di.xyz \
   --xtb $TURBODIR/bin/em64t-unknown-linux-gnu/xtb --workdir xtb_build --outdir .
python3 $S/setup_hse06_opt_template.py --xyz ferul_mono.xyz --charge -1 --output-dir acid
python3 $S/setup_hse06_opt_template.py --xyz ferul_di.xyz   --charge -2 --output-dir base
printf "%s\n%s\n" $R/ferulate/acid $R/ferulate/base > joblist.txt
ls acid base | head -20
