#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

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

#if !defined(DOUBLE_PRECISION)
#define float_t float
#else
#define float_t double
#endif
#define NDIM 2

#define x_axis 0
#define y_axis 1

typedef struct {
  float_t coords[NDIM];
} kpoint;

struct kdnode {
  int axis;
  kpoint split;
  struct kdnode *left, *right;
};

int cmpfunc_x_axis(const void *a, const void *b);
int cmpfunc_y_axis(const void *a, const void *b);

void insertion_sort(kpoint **data, int start, int end,
                    int (*comparator)(const void *, const void *));

int partitioning(kpoint **data, int start, int end,
                 int (*comparator)(const void *, const void *));

void pqsort(kpoint **data, int start, int end,
            int (*comparator)(const void *, const void *));

kpoint *generate_dataset(int len);
void get_dataset_ptrs(kpoint *dataset, kpoint **dataset_ptrs, int len);
void print_array(kpoint **arr, int len);
int cmp_double(const void *a, const void *b);

int main(int argc, char *argv[]) {

  struct timespec ts;

  int numprocs;
  int *provided = malloc(sizeof(int *));
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, provided);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


  int len = 800000;
  kpoint *dataset = generate_dataset(len);

  kpoint **dataset_ptrs = malloc(len * sizeof(kpoint *));
  get_dataset_ptrs(dataset, dataset_ptrs, len);

  double tstart = CPU_TIME;
#pragma omp parallel
  {
#pragma omp single
    pqsort(dataset_ptrs, 0, len, cmpfunc_x_axis);
  }
  double telapsed = CPU_TIME - tstart;

  printf("elapsed time parallel qsort: %f\n", telapsed);

  kpoint *dataset2 = generate_dataset(len);
  kpoint **dataset_ptrs2 = malloc(len * sizeof(kpoint *));

  get_dataset_ptrs(dataset2, dataset_ptrs2, len);
  tstart = CPU_TIME;
  qsort(dataset_ptrs2, len, sizeof(kpoint *), cmpfunc_x_axis);
  telapsed = CPU_TIME - tstart;

  printf("elapsed time serial qsort: %f\n", telapsed);

  return 0;
}

kpoint *generate_dataset(int len) {
  srand48(time(NULL));

  kpoint *dataset = malloc(len * sizeof(kpoint));
  for (int i = 0; i < len; i++) {
    dataset[i].coords[0] = drand48();
    dataset[i].coords[1] = drand48();
  }
  return dataset;
}

int cmp_double(const void *a, const void *b) {
  double *A = (double *)a;
  double *B = (double *)b;

  return (*A >= *B);
}

#define SWAP(A, B, SIZE)                                                       \
  do {                                                                         \
    int sz = (SIZE);                                                           \
    char *a = (A);                                                             \
    char *b = (B);                                                             \
    do {                                                                       \
      char _temp = *a;                                                         \
      *a++ = *b;                                                               \
      *b++ = _temp;                                                            \
    } while (--sz);                                                            \
  } while (0)

void insertion_sort(kpoint **data, int start, int end,
                    int (*comparator)(const void *, const void *)) {
  {
    int min_idx = start;
    for (int i = start + 1; i < end; i++) {
      if (comparator(&data[min_idx], &data[i])) {
        min_idx = i;
      }
    }
    SWAP(&data[start], &data[min_idx], sizeof(&data[start]));
  }

  for (int head = start + 1, run = start + 1; (run = ++head) < end;) {
    while ((run > 0) && comparator(&data[run - 1], &data[run])) {
      SWAP(&data[run - 1], &data[run], sizeof(&data[run]));
      --run;
    }
  }
}

int cmpfunc_x_axis(const void *a, const void *b) {
  kpoint **ptr_a = (kpoint **)a;
  kpoint **ptr_b = (kpoint **)b;

  return ((*ptr_a)->coords[x_axis] >= (*ptr_b)->coords[x_axis]);
}

int cmpfunc_y_axis(const void *a, const void *b) {
  kpoint **ptr_a = (kpoint **)a;
  kpoint **ptr_b = (kpoint **)b;

  return ((*ptr_a)->coords[y_axis] >= (*ptr_b)->coords[y_axis]);
}

void get_dataset_ptrs(kpoint *dataset, kpoint **dataset_ptrs, int len) {
  for (int i = 0; i < len; i++) {
    dataset_ptrs[i] = &dataset[i];
  }
}

inline int compare_ge_x_axis(const void *A, const void *B) {
  kpoint **a = (kpoint **)A;
  kpoint **b = (kpoint **)B;

  return ((*a)->coords[x_axis] >= (*b)->coords[x_axis]);
}

inline int compare_ge_y_axis(const void *A, const void *B) {
  kpoint **a = (kpoint **)A;
  kpoint **b = (kpoint **)B;

  return ((*a)->coords[x_axis] >= (*b)->coords[x_axis]);
}

void print_array(kpoint **arr, int len) {
  for (int i = 0; i < len; i++) {
    printf("element[%d]: (%f,%f)\n", i, arr[i]->coords[0], arr[i]->coords[1]);
  }
}

kpoint **median_of_three(kpoint **a, kpoint **b, kpoint **c,
                         int (*comparator)(const void *, const void *)) {
  if (comparator(b, a) && comparator(c, b))
    return b; // a b c
  if (comparator(c, a) && comparator(b, c))
    return c; // a c b
  if (comparator(a, b) && comparator(c, a))
    return a; // b a c
  if (comparator(c, b) && comparator(a, c))
    return c; // b c a
  if (comparator(a, c) && comparator(b, a))
    return a; // c a b
  return b;   // c b a
}

int partitioning(kpoint **data, int start, int end,
                 int (*comparator)(const void *, const void *)) {

  // pick up the meadian of [0], [mid] and [end] as pivot
  //
  /* to be done */

  // pick up the last element as pivot
  //
  --end;
  void *pivot = &data[end];
  // int mid = ceil((end - start) / 2.0);
  // --end;
  // void *pivot = median_of_three(&data[0], &data[mid], &data[end],
  // comparator);

  int pointbreak = end - 1;
  for (int i = start; i <= pointbreak; i++)
    if (comparator(&data[i], pivot)) {
      while ((pointbreak > i) && comparator(&data[pointbreak], pivot))
        pointbreak--;
      if (pointbreak > i)
        SWAP(&data[i], &data[pointbreak--], sizeof(&data[i]));
    }
  pointbreak += !comparator(&data[pointbreak], pivot);
  SWAP(&data[pointbreak], pivot, sizeof(pivot));

  return pointbreak;
}

#define CHECKSWAP(a, b, comparator)                                            \
  {                                                                            \
    if (comparator(&data[start + (a)], &data[start + (b)]))                    \
      SWAP(&data[start + (a)], &data[start + (b)],                             \
           sizeof(&data[start + (a)]));                                        \
  }

#define task_cutoff 4
#define insertion_cutoff task_cutoff / 2

void pqsort(kpoint **data, int start, int end,
            int (*comparator)(const void *, const void *)) {
  int size = end - start;

  printf("qsort by thread %d\n", omp_get_thread_num());

  if (size == 1) {
    return;
  } else if (size == 2) {
    if (comparator(&data[start], &data[end - 1]))
      SWAP(&data[start], &data[end - 1], sizeof(&data[start]));
  } else if (size == 3) {
    CHECKSWAP(1, 2, comparator);
    CHECKSWAP(0, 2, comparator);
    CHECKSWAP(0, 1, comparator);
  }

  if (size < insertion_cutoff) {
    insertion_sort(data, start, end, comparator);
  } else {
    int mid = partitioning(data, start, end, comparator);

    int mid_start = mid - start;
    if (mid_start > 0)
#pragma omp task default(none) final(mid_start < task_cutoff) mergeable        \
shared(data) firstprivate(start, mid, comparator) untied
      pqsort(data, start, mid, comparator);

    int end_mid = end - (mid + 1);
    if (end_mid)
#pragma omp task default(none) final(end_mid < task_cutoff) mergeable shared(  \
                                                                               \
    data) firstprivate(mid, end, comparator) untied
      pqsort(data, mid + 1, end, comparator);
  }
}
