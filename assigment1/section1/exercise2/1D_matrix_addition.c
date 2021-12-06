#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define SEED 35791250

void fillRandomToMatrix(int x_size, int y_size, int z_size, double *matrix);

int main(int argc, char *argv[])
{
	if (argc != 4)
	{
		printf("Usage: mpirun -np 'number_of_processors' %s 'x_dimension' 'y_dimension' 'z_dimension'", argv[0]);
		return 1;
	}

	int numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	int my_rank, root = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	int matrix_size, matrix_chunk_size;
	int x_size, y_size, z_size;

	x_size = atoi(argv[1]);
	y_size = atoi(argv[2]);
	z_size = atoi(argv[3]);

	matrix_size = x_size * y_size * z_size;
	matrix_chunk_size = matrix_size / numprocs;

	double *matrixA = malloc(matrix_size * sizeof(double));
	double *matrixB = malloc(matrix_size * sizeof(double));
	double *matrixC = malloc(matrix_size * sizeof(double));

	double *chunk_matrixA = malloc(matrix_chunk_size * sizeof(double));
	double *chunk_matrixB = malloc(matrix_chunk_size * sizeof(double));
	double *chunk_matrixC = malloc(matrix_chunk_size * sizeof(double));

	int ndims = 1;
	int dims[1] = {0};
	int periods[1] = {0};
	int reorder = 1;

	MPI_Comm matrix_communicator;
	MPI_Dims_create(numprocs, ndims, dims);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &matrix_communicator);

	if (my_rank == root)
	{
		srand48(SEED);
		fillRandomToMatrix(x_size, y_size, z_size, matrixA);
		fillRandomToMatrix(x_size, y_size, z_size, matrixB);

		// for (int i = 0; i < matrix_chunk_size; i++)
		// {
		// 	printf("%f ", matrixA[i]);
		// }
		// printf("\n");

		// for (int i = 0; i < matrix_chunk_size; i++)
		// {
		// 	printf("%f ", matrixB[i]);
		// }
		// printf("\n");
	}
	MPI_Scatter(matrixA, matrix_chunk_size, MPI_DOUBLE, chunk_matrixA, matrix_chunk_size, MPI_DOUBLE, root, matrix_communicator);
	MPI_Scatter(matrixB, matrix_chunk_size, MPI_DOUBLE, chunk_matrixB, matrix_chunk_size, MPI_DOUBLE, root, matrix_communicator);

	for (int i = 0; i < matrix_chunk_size; i++)
	{
		chunk_matrixC[i] = chunk_matrixA[i] + chunk_matrixB[i];
	}

	MPI_Gather(chunk_matrixC, matrix_chunk_size, MPI_DOUBLE, matrixC, matrix_chunk_size, MPI_DOUBLE, root, matrix_communicator);

	// if (my_rank == root)
	// {
	// 	for (int i = 0; i < matrix_chunk_size; i++)
	// 	{
	// 		printf("%f ", matrixC[i]);
	// 	}
	// 	printf("\n");
	// }

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
