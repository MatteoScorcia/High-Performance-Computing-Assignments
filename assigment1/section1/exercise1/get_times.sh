#! /bin/sh
module load openmpi-4.1.1+gnu-9.3.0
mpicc get_ring_time.c -o ./build/get_ring_time -lm

END=48

for ((i=1;i<=END;i++));
do
	mpirun -np $i ./build/get_ring_time
done