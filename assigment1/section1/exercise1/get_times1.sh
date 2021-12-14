#! /bin/sh
mpicc ring_time.c -o ./build/ring_time

END=24

for ((i=1;i<=END;i++));
do
	mpirun -np $i ./build/ring_time
done