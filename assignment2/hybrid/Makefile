CC := mpicc

BUILD_DIR := ./build
SRC_DIR := ./src

all: build_kdtree

build_kdtree:
	$(CC) $(SRC_DIR)/$@.c -o $(BUILD_DIR)/$@.o
	
.PHONY: clean
clean:
	rm -f $(BUILD_DIR)/*
