#! /bin/sh

qsub -l nodes=2:16 -I -q dssc -l walltime=1:00:00

module load openmpi-4.1.1+gbu-9.3.0

CWD=/u/dssc/matteo/High-Performance-Computing-Assignments/assigment2/hybrid/

POINTS=1000

mpicc build_kdtree.v0.c -fopenmp -lm -w -o build_kdtree.v0

mpirun --map-by node --mca btl ^openib --mca pml ucx -np 2 ./build_kdtree.v0 $POINTS


#for version 0
$NUMTHREADS=8

export OMP_NUM_THREADS=$NUMTHREADS

export OMP_PLACES=sockets
