#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --tasks-per-node=1
#SBATCH --ntasks=1          #Numero total de tarefas MPI
#SBATCH -p batch            #Fila (partition) a ser utilizada
#SBATCH -J GFreitasCR       #Nome job
#SBATCH --time=00:10:00     #Obrigatório
#SBATCH --exclusive

executable="../bin/gf.x"
hpct_measure="hpctoolkit-gf.x-measurements"
hpct_struct="gf.x.hpcstruct"

source env_hpctoolkit.sh

ulimit -s unlimited

export OMP_NUM_THREADS=1
cd  $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NUM_NODES = $SLURM_JOB_NUM_NODES"

echo "./${executable}"

hpcrun -t ${executable}
hpcstruct ${executable}
mv ${hpct_measure}-${SLURM_JOBID} ${hpct_measure}
hpcprof -I . -S ${hpct_struct} ${hpct_measure}

exit




