#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --job-name=Force-calculation-N1000
#SBATCH --nodes=24
#SBATCH --ntasks=532                  # Number of MPI tasks (i.e. processes)
#SBATCH --cpus-per-task=1            # Number of cores per MPI task 
#SBATCH --ntasks-per-node=28
#SBATCH --mem-per-cpu 1GB
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --array=0-7
#SBATCH --output=./logs/MD_N50000_trial2_%a.out
#SBATCH --error=./logs/MD_N50000_trial2_%a.out
#SBATCH -p large-28core

module load shared
module load slurm/17.11.12
module load gcc/7.1.0
module load ucx/1.8.1
module load mpich/gcc/3.4.2

N=50000
ProcArray=("125" "216" "343" "512")
echo "Number of CPU:" ${ProcArray[$SLURM_ARRAY_TASK_ID]}

for i in {1..5}
do
    mpirun -n ${ProcArray[$SLURM_ARRAY_TASK_ID]} ./a.out ${N} 72
done
