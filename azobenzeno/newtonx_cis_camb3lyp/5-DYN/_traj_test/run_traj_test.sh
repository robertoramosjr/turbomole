#!/bin/bash
#SBATCH --job-name=nx_cis_test
#SBATCH --time=02:00:00
#SBATCH --output=%x.out
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE
export PATH=$TURBODIR/scripts:$PATH
export PARA_ARCH=SMP
export PARNODES=$SLURM_CPUS_PER_TASK
source $TURBODIR/Config_turbo_env
# NX le a versao do nome do arquivo TURBOMOLE_781TMG -> usar shadow TURBODIR (lido 7.8)
export TURBODIR=$HOME/software/TURBOMOLE_nxshadow
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
G=/opt/gridunesp/dist/gcc/13.2.0
export LD_LIBRARY_PATH=$G/lib64:$G/lib:$HOME/software/NX-NS_patched_install/lib64:$LD_LIBRARY_PATH
export NX=$HOME/software/NX-NS_patched_install/bin
export CIOVERLAP=$HOME/software/newton-x-package-26a/CIOVERLAP/bin
cd $SLURM_SUBMIT_DIR/TRAJ1
s=$(date +%s)
$NX/nx_moldyn --clean > md.log 2>&1
rc=$?
echo "nx_moldyn rc=$rc $(( $(date +%s)-s )) s" > ../timing.txt
exit $rc
