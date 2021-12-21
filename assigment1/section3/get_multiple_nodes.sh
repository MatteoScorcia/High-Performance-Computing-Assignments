#! /bin/sh
module load openmpi-4.1.1+gnu-9.3.0

CWD=/u/dssc/matteo/High-Performance-Computing-Assignments/assigment1/section3/

mpif77 -ffixed-line-length-none $CWD/Jacobi_MPI_vectormode.F -o $CWD/jacoby3D.x

PROCESSES=(12 24 48)

for process in "${PROCESSES[@]}"
do
	mpirun --mca btl ^openib -np $process --map-by node $CWD/jacoby3D.x <$CWD/input.1200 >$CWD/$process-map-node.txt
done