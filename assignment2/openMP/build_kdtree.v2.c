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
void copy_dataset_ptrs(kpoint **dataset_ptrs, kpoint **new_dataset, int len);

// kd-tree build functions
struct kdnode *build_kdtree(kpoint **dataset_ptrs, int len, int axis, int level);
int choose_splitting_dimension(kpoint **dataset_ptrs, int len);
kpoint *choose_splitting_point(kpoint **dataset_ptrs, int len, int chosen_axis);
float_t get_dataset_extent(kpoint **dataset, int len, int axis);
float_t get_max_value_dataset(kpoint **dataset, int len, int axis);
float_t get_min_value_dataset(kpoint **dataset, int len, int axis);


int main(int argc, char *argv[]) {

  struct timespec ts;

  int len = 800000;
  kpoint *dataset = generate_dataset(len);
  
  // kpoint dataset[9] = {{2, 3}, {5, 4}, {9, 6}, {6, 22}, {4, 7},
  //                      {8, 1}, {7, 2}, {8, 9}, {1, 1}};
  // int len = sizeof(dataset) / sizeof(dataset[0]);

  kpoint **dataset_ptrs = malloc(len * sizeof(kpoint *));
  get_dataset_ptrs(dataset, dataset_ptrs, len);

  int nthreads;
  struct kdnode *root;

  #pragma omp parallel shared(nthreads)
  {
    #pragma omp single
    {
      nthreads = omp_get_num_threads();
      printf("num threads: %d\n", nthreads);
    }
  }
  
  get_dataset_ptrs(dataset, dataset_ptrs, len);
  
  double tstart = CPU_TIME;
  #pragma omp parallel shared(dataset_ptrs, root) firstprivate(len, nthreads) 
  {
    #pragma omp single nowait
    {
      root = build_kdtree(dataset_ptrs, len, 0, 0);
    }
  }
  double telapsed = CPU_TIME - tstart;

  printf("elapsed time: %f\n", telapsed);
  printf("root node: split->(%f,%f)\n", root->split.coords[0], root->split.coords[1] );
  printf("left node: split->(%f,%f)\n", root->left->split.coords[0], root->left->split.coords[1] );
  printf("right node: split->(%f,%f)\n", root->right->split.coords[0], root->right->split.coords[1] );

  // free(dataset);
  free(dataset_ptrs);
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

#define build_cutoff 16 

struct kdnode *build_kdtree(kpoint **dataset_ptrs, int len, int axis, int level) {
  // printf("level: %d, thread: %d\n", level, omp_get_thread_num());
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

  #pragma omp taskgroup
  {
    if (chosen_axis == x_axis) {
      pqsort(dataset_ptrs, 0, len, compare_ge_x_axis);
    } else {
      pqsort(dataset_ptrs, 0, len, compare_ge_y_axis);
    }
  }

  // if (chosen_axis == x_axis) {
  //   qsort(dataset_ptrs, len, sizeof(kpoint *), compare_ge_x_axis);
  // } else {
  //   qsort(dataset_ptrs, len, sizeof(kpoint *), compare_ge_y_axis);
  // }

  kpoint *split_point = choose_splitting_point(dataset_ptrs, len, chosen_axis);

  kpoint **left_points, **right_points;

  int median = ceil(len / 2.0);
  int len_left = median - 1;    // length of the left points
  int len_right = len - median; // length of the right points
  
  #pragma omp task shared(dataset_ptrs, left_points, right_points, median) depend(in:dataset_ptrs) depend(out:left_points) depend(out:right_points)
  {
    left_points = &dataset_ptrs[0];       // starting pointer of left_points
    right_points = &dataset_ptrs[median]; // starting pointer of right_points
  }
  
  node->axis = chosen_axis;
  node->split = *split_point;
  #pragma omp task shared(left_points) depend(in:left_points) firstprivate(len_left, chosen_axis, level) if(level <= build_cutoff) mergeable untied
    node->left = build_kdtree(left_points, len_left, chosen_axis, level+1);
  #pragma omp task shared(right_points) depend(in:right_points) firstprivate(len_right, chosen_axis, level) if(level <= build_cutoff) mergeable untied
    node->right = build_kdtree(right_points, len_right, chosen_axis, level+1);

  #pragma omp taskwait
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

float_t get_dataset_extent(kpoint **dataset, int len, int axis) {
  float_t max = get_max_value_dataset(dataset, len, axis);
  float_t min = get_min_value_dataset(dataset, len, axis);

  return (max - min);
}

float_t get_max_value_dataset(kpoint **dataset, int len, int axis) {
  float_t max_value = (*dataset[0]).coords[axis];
  #pragma omp parallel for reduction(max:max_value) schedule(static) proc_bind(close)
  for (int i = 1; i < len; i++) {
    max_value = max_value > dataset[i]->coords[axis] ? max_value : dataset[i]->coords[axis];
  }

  return max_value;
}

float_t get_min_value_dataset(kpoint **dataset, int len, int axis) {
  float_t min_value = (*dataset[0]).coords[axis];
  #pragma omp parallel for reduction(min:min_value) schedule(static) proc_bind(close)
  for (int i = 1; i < len; i++) {
    min_value = min_value < dataset[i]->coords[axis] ? min_value : dataset[i]->coords[axis];
  }

  return min_value;
}

void copy_dataset_ptrs(kpoint **dataset_ptrs, kpoint **new_dataset, int len) {
  for (int i = 0; i < len; i++) {
    new_dataset[i] = dataset_ptrs[i];
  }
}

void get_dataset_ptrs(kpoint *dataset, kpoint **dataset_ptrs, int len) {
  for (int i = 0; i < len; i++) {
    dataset_ptrs[i] = &dataset[i];
  }
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

int partitioning(kpoint **data, int start, int end,
                 int (*comparator)(const void *, const void *)) {
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
#define task_cutoff 64
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
