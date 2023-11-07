#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --tasks-per-node=1
#SBATCH --ntasks=1          # Total MPI tasks
#SBATCH -p batch            # queue (partition)
#SBATCH -J GFreitasCR       # Job name
#SBATCH --time=04:00:00     # Max time - mandatory
#SBATCH --exclusive

executable="../bin/gf.x"

source env_hpctoolkit_gprof.sh

ulimit -s unlimited

export OMP_NUM_THREADS=1
cd  $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NUM_NODES = $SLURM_JOB_NUM_NODES"

echo "Running ${executable} with gprof"

${executable}
gprof --graph ./${executable} > gprof.out  #  other options:  --exec-times --graph --brief --flat-profile

exit




