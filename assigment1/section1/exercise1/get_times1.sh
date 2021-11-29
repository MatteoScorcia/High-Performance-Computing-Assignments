#! /bin/sh
mpicc ring_stream.c -o ./build/ring_stream_Isend.x
mpicc ring_time.c -o ./build/ring_time.x

END=4

for ((i=1;i<=END;i++));
do
	mpirun -np $i ./build/ring_time.x
done