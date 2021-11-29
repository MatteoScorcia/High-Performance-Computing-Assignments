#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
	int myid, numprocs, proc;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	// Size of the default communicator
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	// Ask MPI to decompose our processes in a 2D cartesian grid for us
	int dims[2] = {0, 0};
	MPI_Dims_create(size, 2, dims);
	// Make both dimensions periodic
	int periods[2] = {0, 1};
	// Let MPI assign arbitrary ranks if it deems it necessary
	int reorder = 1;

	// Create a communicator given the 2D torus topology.
	MPI_Comm cartesian_communicator;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartesian_communicator);
	// My rank in the new communicator
	int my_rank;
	MPI_Comm_rank(cartesian_communicator, &my_rank);
	// Get my coordinates in the new communicator
	int my_coords[2];
	MPI_Cart_coords(cartesian_communicator, my_rank, 2, my_coords);
	// Print my location in the 2D torus.
	printf("[MPI process %d] I am located at (%d, %d).\n", my_rank, my_coords[0], my_coords[1]);

	int direction = 0, disposition = 1;
	int rank_source, rank_dest;

	MPI_Cart_shift(cartesian_communicator, direction, disposition, &rank_source, &rank_dest);
	printf("Shifting -> rank_source: %d, rank_dest: %d \n", rank_source, rank_dest);

	MPI_Finalize();
	return EXIT_SUCCESS;
}