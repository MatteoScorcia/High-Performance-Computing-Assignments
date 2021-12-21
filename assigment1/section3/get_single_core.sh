#! /bin/sh
module load openmpi

CWD=/u/dssc/matteo/High-Performance-Computing-Assignments/assigment1/section3/

mpif77 -ffixed-line-length-none $CWD/Jacobi_MPI_vectormode.F -o $CWD/jacoby3D.x

mpirun -np 1 $CWD/jacoby3D.x