CC := mpicc

BUILD_DIR := ./build
SRC_DIR := ./src

all: ring_stream_Isend matrix_addition

ring_stream_Isend:
	$(CC) $(SRC_DIR)/$@.c -o $(BUILD_DIR)/$@.o
	
matrix_addition:
	$(CC) $(SRC_DIR)/$@.c -o $(BUILD_DIR)/$@.o

.PHONY: clean
clean:
	rm -f $(BUILD_DIR)/*