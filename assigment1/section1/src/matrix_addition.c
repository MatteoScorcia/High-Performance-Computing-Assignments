#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define SEED 35791250

void fillRandomToMatrix(int x_size, int y_size, int z_size, double *matrix);

int main(int argc, char *argv[])
{
	int numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	int my_rank, root = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if (my_rank == root)
	{
		if ((argc > 7) || (argc < 4) || (argc == 5))
		{
			printf("Usage: mpirun -np 'number_of_processors' %s 'x_dimension' 'y_dimension' 'z_dimension' 'x_topology_dimension(optional)' 'y_topology_dimension(optional)) 'z_topology_dimension(optional)'", argv[0]);
			return 1;
		}
	}

	int matrix_size, matrix_chunk_size;
	int x_size, y_size, z_size;
	int x_topo_dim, y_topo_dim, z_topo_dim;

	x_size = atoi(argv[1]);
	y_size = atoi(argv[2]);
	z_size = atoi(argv[3]);

	int ndims, *dims, *periods, reorder;

	if (argc == 4)
	{
		dims = malloc(sizeof(int));
		periods = malloc(sizeof(int));

		ndims = 1;
		dims[0] = 0;
		periods[0] = 0;
		reorder = 1;
	}
	else if (argc == 6)
	{
		x_topo_dim = atoi(argv[4]);
		y_topo_dim = atoi(argv[5]);

		dims = malloc(2 * sizeof(int));
		periods = malloc(2 * sizeof(int));

		ndims = 2;
		dims[0] = x_topo_dim;
		dims[1] = y_topo_dim;
		periods[0] = 0;
		periods[1] = 0;
		reorder = 1;
	}
	else if (argc == 7)
	{
		x_topo_dim = atoi(argv[4]);
		y_topo_dim = atoi(argv[5]);
		z_topo_dim = atoi(argv[6]);

		dims = malloc(3 * sizeof(int));
		periods = malloc(3 * sizeof(int));

		ndims = 3;
		dims[0] = x_topo_dim;
		dims[1] = y_topo_dim;
		dims[2] = z_topo_dim;
		periods[0] = 0;
		periods[1] = 0;
		periods[2] = 0;
		reorder = 1;
	}

	matrix_size = x_size * y_size * z_size;
	matrix_chunk_size = matrix_size / numprocs;

	double *matrixA = malloc(matrix_size * sizeof(double));
	double *matrixB = malloc(matrix_size * sizeof(double));
	double *matrixC = malloc(matrix_size * sizeof(double));

	double *chunk_matrixA = malloc(matrix_chunk_size * sizeof(double));
	double *chunk_matrixB = malloc(matrix_chunk_size * sizeof(double));
	double *chunk_matrixC = malloc(matrix_chunk_size * sizeof(double));

	MPI_Comm matrix_communicator;
	MPI_Dims_create(numprocs, ndims, dims);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &matrix_communicator);

	if (my_rank == root)
	{
		srand48(SEED);
		fillRandomToMatrix(x_size, y_size, z_size, matrixA);
		fillRandomToMatrix(x_size, y_size, z_size, matrixB);
	}
	MPI_Scatter(matrixA, matrix_chunk_size, MPI_DOUBLE, chunk_matrixA, matrix_chunk_size, MPI_DOUBLE, root, matrix_communicator);
	MPI_Scatter(matrixB, matrix_chunk_size, MPI_DOUBLE, chunk_matrixB, matrix_chunk_size, MPI_DOUBLE, root, matrix_communicator);

	for (int i = 0; i < matrix_chunk_size; i++)
	{
		chunk_matrixC[i] = chunk_matrixA[i] + chunk_matrixB[i];
	}

	MPI_Gather(chunk_matrixC, matrix_chunk_size, MPI_DOUBLE, matrixC, matrix_chunk_size, MPI_DOUBLE, root, matrix_communicator);

	free(matrixA);
	free(matrixB);
	free(matrixC);

	free(chunk_matrixA);
	free(chunk_matrixB);
	free(chunk_matrixC);

	MPI_Finalize();
	return 0;
}

void fillRandomToMatrix(int x_size, int y_size, int z_size, double *matrix)
{
	int i, j, k;

	for (i = 0; i < x_size; i++)
		for (j = 0; j < y_size; j++)
			for (k = 0; k < z_size; k++)
				matrix[(y_size * z_size * i) + (z_size * j) + k] = drand48();
}
