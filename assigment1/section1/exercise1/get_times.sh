#! /bin/sh
module load openmpi-4.1.1+gnu-9.3.0
mpicc get_ring_time.c -o ./build/get_ring_time

PROCESSORS=48
ITERATIONS=50
CWD="~/High-Performance-Computing-Assignments/assignment1/section1/exercise1/"

for ((i=1;i<=PROCESSORS;i++));
do for((j=1;j<=ITERATIONS;j++));
	do
		mpirun -np $i $CWD/build/get_ring_time $j
	done
done
