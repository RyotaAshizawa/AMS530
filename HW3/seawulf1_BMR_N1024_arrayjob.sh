#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --job-name=BMR_1024
#SBATCH --nodes=4
#SBATCH --ntasks=112                  # Number of MPI tasks (i.e. processes)
#SBATCH --cpus-per-task=1            # Number of cores per MPI task 
#SBATCH --ntasks-per-node=28
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --array=0-2
#SBATCH --output=./logs/BMR_1024_%a.out
#SBATCH --error=./logs/BMR_1024_%a.out
#SBATCH -p long-28core

module load shared
module load slurm/17.11.12
module load gcc/7.1.0
module load ucx/1.8.1
module load mpich/gcc/3.4.2

ProcArray=("4" "16" "64")
echo "Number of CPU:" ${ProcArray[$SLURM_ARRAY_TASK_ID]}

for i in {1..5}
do
    echo "Loop:" $i
    mpirun -n ${ProcArray[$SLURM_ARRAY_TASK_ID]} ./BMR.o 1024 0
done
