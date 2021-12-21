#! /bin/sh
module load openmpi-4.1.1+gnu-9.3.0

CWD=/u/dssc/matteo/High-Performance-Computing-Assignments/assigment1/section3/

mpif77 -ffixed-line-length-none $CWD/Jacobi_MPI_vectormode.F -o $CWD/jacoby3D.x

mpirun --mca btl ^openib -np 1 $CWD/jacoby3D.x <$CWD/input.1200 >$CWD/single_core_thin.txt