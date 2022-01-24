#! /bin/sh
module load openmpi-4.1.1+gnu-9.3.0

CWD=/u/dssc/matteo/High-Performance-Computing-Assignments/assigment1/section1/exercise2

mpicc $CWD/get_matrix_time.c -o $CWD/build/get_matrix_time

topology=(4 3 2)
distribution=(800 300 100)

for processors in {1..24}
do
	mpirun -np $processors $CWD/build/get_matrix_time ${distribution[@]} ${topology[@]}
done
