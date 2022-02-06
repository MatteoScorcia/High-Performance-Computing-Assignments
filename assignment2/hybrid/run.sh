#! /bin/sh
module load openmpi-4.1.1+gbu-9.3.0

CWD=/u/dssc/matteo/High-Performance-Computing-Assignments/assigment2/hybrid/

POINTS=1000

mpicc build_kdtree.v0.c -fopenmp -lm -w -o build_kdtree.v0

mpirun --map-by node --mca btl ^openib --mca pml ucx -np 2 ./build_kdtree.v0 $POINTS