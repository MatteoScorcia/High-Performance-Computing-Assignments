#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#define SEED 35791246

int main(int argc, char *argv[])
{
	double start_time, end_time;
	int myid, numprocs, proc;
	MPI_Status status;
	MPI_Request request;

	int master = 0;
	int tag = 123;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	if (argc <= 1)
	{
		fprintf(stderr, " Usage : mpi -np n %s number_of_iterations \n", argv[0]);
		MPI_Finalize();
		exit(-1);
	}

	long long int N = atoll(argv[1]);
	long long int number_local_operations = N / numprocs;
	long long int local_array[number_local_operations];
	long long int local_sum;

	long long int i;

	srand48(SEED * (myid + 1));

	for (i = 0; i < number_local_operations; i++)
	{
		local_array[i] = rand();
	}

	for (i = 0; i < number_local_operations; i++)
	{
		local_sum += local_array[i];
	}

	if (myid == 0)
	{
		long long int global_sum = local_sum;

		for (proc = 1; proc < numprocs; proc++)
		{
			MPI_Recv(&local_sum, 1, MPI_LONG_LONG, proc, tag, MPI_COMM_WORLD, &status);

			global_sum += local_sum;
		}

		printf("\nfinal sum = %llu \n", global_sum);
	}
	else
	{
		MPI_Send(&local_array, 1, MPI_LONG_LONG, master, tag, MPI_COMM_WORLD);
	}

	MPI_Finalize();
}