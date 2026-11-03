#!/bin/bash
#SBATCH --job-name=nx_test_fhop
#SBATCH --time=08:00:00
#SBATCH --output=TRAJECTORIES/TRAJ%a/slurm.out
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
# NX le a versao TM pelo nome do arquivo -> shadow TURBODIR (7.8) + NX patchado
export TURBODIR=$HOME/software/TURBOMOLE_nxshadow
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
G=/opt/gridunesp/dist/gcc/13.2.0
export LD_LIBRARY_PATH=$G/lib64:$G/lib:$HOME/software/NX-NS_patched_install_fhop/lib64:$LD_LIBRARY_PATH
export NX=$HOME/software/NX-NS_patched_install_fhop/bin
export CIOVERLAP=$HOME/software/newton-x-package-26a/CIOVERLAP/bin
cd $SLURM_SUBMIT_DIR/TRAJECTORIES/TRAJ${SLURM_ARRAY_TASK_ID}
# --clean only on a fresh start; if a previous attempt exists (timeout), use nx_restart manually
$NX/nx_moldyn --clean > md.log 2>&1
rc=$?
grep -q 'Normal termination' md.log || rc=1
echo "rc=$rc" > rc.txt
exit $rc
