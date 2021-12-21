#! /bin/sh
module load openmpi-4.1.1+gnu-9.3.0

CWD=/u/dssc/matteo/High-Performance-Computing-Assignments/assigment1/section3/

mpif77 -ffixed-line-length-none $CWD/Jacobi_MPI_vectormode.F -o $CWD/jacoby3D.x

PROCESSES=(4 8 12)

for process in "${PROCESSES[@]}"
do
	mpirun --mca btl ^openib -np $process --map-by core $CWD/jacoby3D.x <$CWD/input.1200 2>/dev/null
	mpirun --mca btl ^openib -np $process --map-by socket $CWD/jacoby3D.x <$CWD/input.1200 2>/dev/null
done