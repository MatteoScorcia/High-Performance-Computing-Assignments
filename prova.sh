#! /bin/sh
module load openmpi-4.1.1+gnu-9.3.0

PROCESSORS=48
ITERATIONS=50
CWD=/u/dssc/matteo/High-Performance-Computing-Assignments/assigment1/section1/exercise1
mpicc $CWD/get_ring_time.c -o $CWD/build/get_ring_time

touch /u/dssc/matteo/High-Performance-Computing-Assignments/output.txt

for ((i=1;i<=PROCESSORS;i++));
do for((j=1;j<=ITERATIONS;j++));
	do
		echo "mpirun -np $i $CWD/build/get_ring_time $j" >> "/u/dssc/matteo/High-Performance-Computing-Assignments/output.txt"
	done
done
