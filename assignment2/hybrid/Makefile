CC := mpicc

BUILD_DIR := ./build
SRC_DIR := ./src

all: build_kdtree

build_kdtree:
	$(shell module load openmpi-4.1.1+gnu-9.3.0)
	$(CC) $(SRC_DIR)/$@.c -fopenmp -lm -o $(BUILD_DIR)/$@.o
	
.PHONY: clean
clean:
	rm -f $(BUILD_DIR)/*
