#!/bin/bash
#
#SBATCH --partition=compute2011
#SBATCH --exclusive
#SBATCH --output=test10_%j.out
echo we begin our execution...

# load the modules
module load mpi/openmpi/

date
mpicc exercise_04_2_b.c

for j in {1..10}
do
	for i in {2..16}
	do
		echo running with $i processors
		echo 0.0 1.0 10*$i | mpirun -np $i ./a.out >> results
	done
done
date
