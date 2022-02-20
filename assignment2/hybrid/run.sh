#! /bin/sh

# qsub -l nodes=2:ppn=16 -I -q dssc -l walltime=1:00:00

CWD=/u/dssc/matteo/High-Performance-Computing-Assignments/assignment2/hybrid/

POINTS=(10000000)
NUMTHREADS=8
MPI_PROCS=(1 2 4)

export OMP_NUM_THREADS=$NUMTHREADS

export OMP_PLACES=sockets

module load openmpi-4.1.1+gnu-9.3.0

mpicc $CWD/build_kdtree.v0.c -fopenmp -lm -w -o $CWD/build_kdtree.v0

for procs in "${MPI_PROCS[@]}"
do
  mpirun --map-by node:PE=$NUMTHREADS --mca btl ^openib --mca pml ucx -np $procs ./build_kdtree.v0 ${POINTS[@]} >> out.txt
done
