#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <math.h>

#if !defined(DOUBLE_PRECISION)
#define float_t float
#else
#define float_t double
#endif
#define NDIM 2

typedef struct kpoint_s {
  float_t coords[NDIM];
} kpoint;

kpoint *generate_dataset(int len);

int main(int argc, char *argv[])
{
	int numprocs;
  int *provided = malloc(sizeof(int *));
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, provided);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
 
  printf("mpi process %d init!\n", my_rank);
  
  kpoint dataset[9] = {{2, 3}, {5, 4}, {9, 6}, {6, 22}, {4, 7},
                       {8, 1}, {7, 2}, {8, 9}, {1, 1}};
  int len = sizeof(dataset) / sizeof(dataset[0]);

  kpoint extremes[2] = {{1,9}, {1,22}};

  kpoint computed_median = {5, 11};

  #pragma omp parallel 
  {
    #pragma omp single 
    {
      printf("i am thread %d, there are %d threads\n", omp_get_thread_num(), omp_get_num_threads());
    }
  }

  if (my_rank == 0) {
    MPI_Send(&dataset[2], 4, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
  }

  if (my_rank == 1) {
    kpoint received[2];
    MPI_Status status;
    MPI_Recv(&received, 4, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    printf("received: %f,%f and %f,%f\n", received[0].coords[0], received[0].coords[1], received[1].coords[0], received[1].coords[1]);
  }

  if (my_rank == 0) {
    printf("test parallel for...\n");
    kpoint current_median = dataset[0];
    int current_median_idx = 0;

    #pragma omp parallel for shared(dataset, computed_median, current_median) schedule(static) proc_bind(close)
    for (int i = 1; i < len; i++) {
      if (fabs(dataset[i].coords[1] - computed_median.coords[1]) < fabs(dataset[i].coords[1] - current_median.coords[1])) {
        printf("|%f[%d] - %f| = %f \n", dataset[i].coords[1], i,  computed_median.coords[1], fabs(dataset[i].coords[1] - computed_median.coords[1]));
        current_median = dataset[i];
        current_median_idx = i;
      }
    }
    printf("median for y axis of whole dataset is %f, index is %d\n", current_median.coords[1], current_median_idx);
  }

	MPI_Finalize();
	return 0;
}

kpoint *generate_dataset(int len) {
  srand((unsigned int)time(NULL));

  kpoint *dataset = malloc(len * sizeof(kpoint));
  for (int i = 0; i < len; i++) {
    dataset[i].coords[0] = (float_t)drand48();
    dataset[i].coords[1] = (float_t)drand48();
  }
  return dataset;
}
