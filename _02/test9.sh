#!/bin/bash
#
#SBATCH --partition=compute2011
#SBATCH --exclusive
#SBATCH --output=test9_%j.out
echo we begin our execution...

# load the modules
module load mpi/openmpi/

date
mpicc exercise_04_2.c
echo runing with 8 processors 
mpirun -np 8 ./a.out < inputs.txt
date
