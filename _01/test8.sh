#!/bin/bash
#
#SBATCH --partition=compute2011
#SBATCH --exclusive
#SBATCH --output=output.out
echo we begin our execution...

# load the modules
module load mpi/openmpi/

date
mpicc exercise_04_1.c

for j in {1..50}
do
	for i in {2..16}
	do
		echo running with $i processors
		echo 0.0 100*$i | mpirun -np $i ./a.out >> results
	done
done
date
