#if defined(__STDC__)
#if (__STDC_VERSION__ >= 199901L)
#define _XOPEN_SOURCE 700
#endif
#endif
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#if defined(_OPENMP)
#define CPU_TIME                                                               \
  (clock_gettime(CLOCK_REALTIME, &ts),                                         \
   (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9)

#define CPU_TIME_th                                                            \
  (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &myts),                              \
   (double)myts.tv_sec + (double)myts.tv_nsec * 1e-9)

#else

#define CPU_TIME                                                               \
  (clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts),                               \
   (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9)
#endif

int main(int argc, char *argv[]) {
  int sum = 0;
  int N = 800000;
  int nthreads;

  double *array = malloc(N * sizeof(double));

  for (int ii = 0; ii < N; ii++)
    array[ii] = (double)ii;

#pragma omp parallel
#pragma omp master
  nthreads = omp_get_num_threads();

  printf("starting with %d threads\n", nthreads);

  struct timespec ts;
  double start_time = CPU_TIME;

  double res1 = 0, res2 = 0, res3 = 0;

#pragma omp task shared(array, N, res1) untied
{
  int local_result = 0;
  for (int i = 0; i < N; i++) {
    local_result += i;
  }
  #pragma omp atomic
  res1 += local_result;
}

#pragma omp task shared(array, N, res2) untied
{
  int local_result = 0;
  for (int i = 0; i < N; i++) {
    local_result -= i;
  }
  #pragma omp atomic
  res2 += local_result;
}

#pragma omp taskwait
printf("res1: %f\n", res1);
printf("res2: %f\n", res2);
printf("tasks have finished the work!!\n");

  double finish_time = CPU_TIME;
  printf("elapsed time: %f\n", finish_time - start_time);
  return 0;
}
