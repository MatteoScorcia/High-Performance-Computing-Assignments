#! /bin/sh
module load openmpi

CWD=/u/dssc/matteo/High-Performance-Computing-Assignments/assigment1/section3/

mpif77 -ffixed-line-length-none $CWD/Jacobi_MPI_vectormode.F -o $CWD/jacoby3D.x

PROCESSES=(12 24 48)

for process in "${PROCESSES[@]}"
do
	mpirun -np $process --map-by node $CWD/jacoby3D.x 2>/dev/null
done