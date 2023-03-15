#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --tasks-per-node=1
#SBATCH --ntasks=1               #Numero total de tarefas MPI
#SBATCH -p proc             #Fila (partition) a ser utilizada
#SBATCH -J GFreitasCR       #Nome job
#SBATCH --time=00:30:00          #Obrigatório
######sbatch --mem-per-cpu=64000M
#SBATCH --exclusive

executable=gf.x

module list

resultdir=/mnt/beegfs/carlos.souza/GF_standalone/results/partition-${SLURM_JOB_PARTITION}/NUMNODES-$SLURM_JOB_NUM_NODES/MPI-${SLURM_NTASKS}/JOBID-${SLURM_JOBID}

mkdir -p ${resultdir}

ulimit -s unlimited

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export I_MPI_DEBUG=5
export MKL_DEBUG_CPU_TYPE=5
export I_MPI_FABRICS=shm:ofi

cd  $SLURM_SUBMIT_DIR

echo "SLURM_JOB_NUM_NODES = $SLURM_JOB_NUM_NODES"

date
echo "./${executable}"
time ./${executable}
date

cp slurm-${SLURM_JOBID}.out ${resultdir}/
#cp log.atmosphere.*.out ${resultdir}/
#cp namelist.atmosphere ${resultdir}/
cp submit*.sh ${resultdir}/
#mv x1.*.init.nc-*.lock ${resultdir}/
#mv diag* ${resultdir}/
#mv histor* ${resultdir}/ 

exit




