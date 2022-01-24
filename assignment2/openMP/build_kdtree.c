#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#if defined(_OPENMP)
#define CPU_TIME (clock_gettime( CLOCK_REALTIME, &ts ), (double)ts.tv_sec + \
		  (double)ts.tv_nsec * 1e-9)

#define CPU_TIME_th (clock_gettime( CLOCK_THREAD_CPUTIME_ID, &myts ), (double)myts.tv_sec +	\
		     (double)myts.tv_nsec * 1e-9)

#else

#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
		  (double)ts.tv_nsec * 1e-9)
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

// utility functions
void pqsort(kpoint **data, int start, int end,
            int (*comparator)(const void *, const void *));
int compare_ge_x_axis(const void *a, const void *b);
int compare_ge_y_axis(const void *a, const void *b);
kpoint *generate_dataset(int len);

void get_dataset_ptrs(kpoint *dataset, kpoint **dataset_ptrs, int len);
void copy_dataset_ptrs(kpoint **dataset_ptrs, kpoint **new_arr, int len);
void print_dataset(kpoint *dataset, int len);
void print_dataset_ptr(kpoint **dataset_ptr, int len);

// kd-tree build functions
struct kdnode *build_kdtree(kpoint **dataset_ptr, int len, int ndim, int axis, int nthreads);
int choose_splitting_dimension(kpoint **dataset_ptrs, int len);
kpoint *choose_splitting_point(kpoint **dataset_ptrs, int len, int chosen_axis);
float_t get_dataset_extent(kpoint **arr, int len, int axis);
float_t get_max_value_dataset(kpoint **arr, int len, int axis);
float_t get_min_value_dataset(kpoint **arr, int len, int axis);

int main(int argc, char *argv[]) {

  struct timespec ts;

  int len = 1000000;
  kpoint *dataset = generate_dataset(len);

  kpoint **dataset_ptrs = malloc(len * sizeof(kpoint *));
  get_dataset_ptrs(dataset, dataset_ptrs, len);
  
  int nthreads;
  #pragma omp parallel shared(nthreads)
  {
    #pragma omp single
    {
      int nthreads = omp_get_num_threads();
      printf("num threads: %d\n", nthreads);
    }
  }
  
  get_dataset_ptrs(dataset, dataset_ptrs, len);
  
  double tstart = CPU_TIME;
  build_kdtree(dataset_ptrs, len, 2, 0, nthreads);
  double telapsed = CPU_TIME - tstart;

  printf("elapsed time: %f\n", telapsed);
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

#define task_cutoff 16 

struct kdnode *build_kdtree(kpoint **dataset_ptrs, int len, int ndim,
                            int axis, int nthreads) {
  printf("input len: %d\n", len);
  if (len == 1) {
    struct kdnode *leaf = malloc(sizeof(struct kdnode));
    leaf->axis = axis;
    leaf->split = *dataset_ptrs[0];
    leaf->left = NULL;
    leaf->right = NULL;

    return leaf;
  } else if (len == 0) {
    return NULL;
  }

  struct kdnode *node = malloc(sizeof(struct kdnode));

  int chosen_axis = choose_splitting_dimension(dataset_ptrs, len);

  if (chosen_axis == x_axis) {
    pqsort(dataset_ptrs, 0, len, compare_ge_x_axis);
  } else {
    pqsort(dataset_ptrs, 0, len, compare_ge_y_axis);
  }

  kpoint *split_point = choose_splitting_point(dataset_ptrs, len, chosen_axis);

  kpoint **left_points, **right_points;

  int median = ceil(len / 2.0);
  int len_left = median - 1;    // length of the left points
  int len_right = len - median; // length of the right points

  left_points = &dataset_ptrs[0];       // starting pointer of left_points
  right_points = &dataset_ptrs[median]; // starting pointer of right_points

  printf("len_left = %d, len_right = %d\n", len_left, len_right);

  node->axis = chosen_axis;
  node->split = *split_point;
  #pragma omp task shared(nthreads, ndim, chosen_axis) firstprivate(left_points, len_left) mergeable untied 
    node->left = build_kdtree(left_points, len_left, ndim, chosen_axis, nthreads);
  #pragma omp task shared(nthreads, ndim, chosen_axis) firstprivate(right_points, len_right) mergeable untied
    node->right = build_kdtree(right_points, len_right, ndim, chosen_axis, nthreads);
  return node;
}

kpoint *choose_splitting_point(kpoint **ordered_dataset, int len,
                               int chosen_axis) {
  int median_idx = ceil(len / 2.0) - 1;
  return ordered_dataset[median_idx];
}

int choose_splitting_dimension(kpoint **dataset_ptrs, int len) {
  float_t x_extent = get_dataset_extent(dataset_ptrs, len, x_axis);
  float_t y_extent = get_dataset_extent(dataset_ptrs, len, y_axis);

  if (x_extent > y_extent) {
    return x_axis;
  }
  return y_axis;
}

float_t get_dataset_extent(kpoint **arr, int len, int axis) {
  float_t max = get_max_value_dataset(arr, len, axis);
  float_t min = get_min_value_dataset(arr, len, axis);

  return (max - min);
}

float_t get_max_value_dataset(kpoint **arr, int len, int axis) {
  float_t max_value = (*arr[0]).coords[axis];
  #pragma omp parallel for reduction(max:max_value) schedule(static) proc_bind(close)
  for (int i = 1; i < len; i++) {
    max_value = max_value > arr[i]->coords[axis] ? max_value : arr[i]->coords[axis];
  }

  return max_value;
}

float_t get_min_value_dataset(kpoint **arr, int len, int axis) {
  float_t min_value = (*arr[0]).coords[axis];
  #pragma omp parallel for reduction(min:min_value) schedule(static) proc_bind(close)
  for (int i = 1; i < len; i++) {
    min_value = min_value < arr[i]->coords[axis] ? min_value : arr[i]->coords[axis];
  }

  return min_value;
}

void copy_dataset_ptrs(kpoint **dataset_ptrs, kpoint **new_arr, int len) {
  for (int i = 0; i < len; i++) {
    new_arr[i] = dataset_ptrs[i];
  }
}

void get_dataset_ptrs(kpoint *dataset, kpoint **dataset_ptrs, int len) {
  for (int i = 0; i < len; i++) {
    dataset_ptrs[i] = &dataset[i];
  }
}

void print_dataset_ptr(kpoint **dataset_ptr, int len) {
  for (int i = 0; i < len; i++) {
    printf("dataset[%d] -> (%f,%f)\n", i, dataset_ptr[i]->coords[0],
           dataset_ptr[i]->coords[1]);
  }
}

void print_dataset(kpoint *dataset, int len) {
  for (int i = 0; i < len; i++) {
    printf("dataset[%d] -> (%f,%f)\n", i, dataset[i].coords[0],
           dataset[i].coords[1]);
  }
  printf("\n");
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

#define CHECKSWAP(a, b, comparator)                                            \
  {                                                                            \
    if (comparator(&data[start + (a)], &data[start + (b)]))                    \
      SWAP(&data[start + (a)], &data[start + (b)],                             \
           sizeof(&data[start + (a)]));                                        \
  }

inline int compare_ge_x_axis(const void *A, const void *B) {
  kpoint **a = (kpoint **)A;
  kpoint **b = (kpoint **)B;

  return ((*a)->coords[x_axis] >= (*b)->coords[x_axis]);
}

inline int compare_ge_y_axis(const void *A, const void *B) {
  kpoint **a = (kpoint **)A;
  kpoint **b = (kpoint **)B;

  return ((*a)->coords[y_axis] >= (*b)->coords[y_axis]);
}

kpoint **median_of_three(kpoint **a, kpoint **b, kpoint **c, int (*comparator)(const void *, const void *)) {
  if (comparator(b,a) && comparator(c,b)) return b;  // a b c
  if (comparator(c,a) && comparator(b,c)) return c;  // a c b
  if (comparator(a,b) && comparator(c,a)) return a;  // b a c
  if (comparator(c,b) && comparator(a,c)) return c;  // b c a
  if (comparator(a,c) && comparator(b,a)) return a;  // c a b
  return b;                                          // c b a
}

int partitioning(kpoint **data, int start, int end,
                 int (*comparator)(const void *, const void *)) {
  // pick up the meadian of [0], [mid] and [end] as pivot
  // int mid = ceil((end - start) / 2.0);
  // --end;
  // void *pivot = median_of_three(&data[0], &data[mid], &data[end], comparator); 
  --end;
  void *pivot = &data[end];

  int pointbreak = end - 1;
  for (int i = start; i <= pointbreak; i++)
    if (comparator(&data[i], pivot)) {
      while ((pointbreak > i) && comparator(&data[pointbreak], pivot))
        pointbreak--;
      if (pointbreak > i)
        SWAP(&data[i], &data[pointbreak--], sizeof(&data[i]));
    }
  pointbreak += !comparator(&data[pointbreak], pivot);
  SWAP(&data[pointbreak], pivot, sizeof(&data[pointbreak]));

  return pointbreak;
}

void insertion_sort(kpoint **data, int start, int end,
                    int (*comparator)(const void *, const void *)) {
  {
    int min_idx = start;
    for (int i = start + 1; i < end; i++)
      if (comparator(&data[min_idx], &data[i]))
        min_idx = i;

    SWAP(&data[start], &data[min_idx], sizeof(&data[start]));
  }

  for (int head = start + 1, run = start + 1; (run = ++head) < end;) {
    while ((run > 0) && comparator(&data[run - 1], &data[run])) {
      SWAP(&data[run - 1], &data[run], sizeof(&data[run - 1]));
      --run;
    }
  }
}
#define task_cutoff 4
#define insertion_cutoff task_cutoff / 2

void pqsort(kpoint **data, int start, int end,
            int (*comparator)(const void *, const void *)) {
  int size = end - start;

  switch (size) {
  case 1:
    break;
  case 2: {
    if (comparator(&data[start], &data[end - 1]))
      SWAP(&data[start], &data[end - 1], sizeof(&data[start]));
  } break;
  case 3: {
    CHECKSWAP(1, 2, comparator);
    CHECKSWAP(0, 2, comparator);
    CHECKSWAP(0, 1, comparator);
  } break;
  default: {
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
    data) firstprivate(mid, end, comparator) untied
        pqsort(data, mid + 1, end, comparator);
    }
  } break;
  }
}
