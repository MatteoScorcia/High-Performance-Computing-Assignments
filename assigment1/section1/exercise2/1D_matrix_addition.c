#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define SEED 35791250

void fillRandomToMatrix(int x_size, int y_size, int z_size, double result_matrix[x_size][y_size][z_size]);

int main(int argc, char const *argv[])
{
	if (argc != 4)
	{
		printf("Usage: mpirun -np 'number_of_processors' %s 'x_dimension' 'y_dimension' 'z_dimension' ", argv[0]);
		return 1;
	}

	int x_size = atoi(argv[1]);
	int y_size = atoi(argv[2]);
	int z_size = atoi(argv[3]);

	double matrixA[x_size][y_size][z_size];
	double matrixB[x_size][y_size][z_size];

	fillRandomToMatrix(x_size, y_size, z_size, matrixA);
	fillRandomToMatrix(x_size, y_size, z_size, matrixB);

	int numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	return 0;
}

void fillRandomToMatrix(int x_size, int y_size, int z_size, double result_matrix[x_size][y_size][z_size])
{
	srand48(SEED); // always the same for each mpi process..
	int i, j, k;

	for (i = 0; i < x_size; i++)
		for (j = 0; j < y_size; j++)
			for (k = 0; k < z_size; k++)
				result_matrix[i][j][k] = drand48();
}