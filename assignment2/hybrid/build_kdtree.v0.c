#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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
#define FLOAT_T MPI_FLOAT
#else
#define float_t double
#define FLOAT_T MPI_DOUBLE
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
            int (*comparator)(const void *, const void *),
            int (*comparator_insort)(const void *, const void *));
int compare_ge_x_axis(const void *a, const void *b);
int compare_ge_y_axis(const void *a, const void *b);
int compare_g_x_axis(const void *a, const void *b);
int compare_g_y_axis(const void *a, const void *b);

kpoint *generate_dataset(int len);
void get_dataset_ptrs(kpoint *dataset, kpoint **dataset_ptrs, int len);
void copy_dataset_ptrs(kpoint **dataset_ptrs, kpoint **new_dataset, int len);
void copy_dataset_from_ptrs(kpoint *new_dataset, kpoint **dataset_ptrs,
                            int len);

// kd-tree build functions
struct kdnode *build_kdtree(kpoint **dataset_ptrs, kpoint extremes[NDIM],
                            int len, int axis, int level);
struct kdnode *build_kdtree_until_level(kpoint **dataset_ptrs,
                                        kpoint extremes[NDIM], int len,
                                        int axis, int level, int final_level,
                                        int is_root_proc, int *is_proc_free);
int choose_splitting_dimension(kpoint extremes[NDIM]);
kpoint *choose_splitting_point(kpoint **dataset_ptrs, int len, int chosen_axis);
void get_dataset_extremes(kpoint **dataset_ptrs, kpoint extremes[NDIM], int len,
                          int axis);
void copy_extremes(kpoint old_extremes[NDIM], kpoint new_extremes[NDIM]);
void send_dataset_to_free_process(int dataset_len, kpoint **dataset_ptrs,
                                  kpoint extremes[NDIM], int previous_axis,
                                  int final_level, int *is_proc_free);
void recv_dataset_from_root_process(int *recv_len, int *recv_axis,
                                    kpoint *recv_dataset,
                                    kpoint *recv_extremes);

int main(int argc, char *argv[]) {

  struct timespec ts;
  double tstart;
  int len;
  int nthreads;

  int numprocs;
  int *provided = malloc(sizeof(int *));
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, provided);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  tstart = CPU_TIME;

#pragma omp parallel shared(nthreads)
  {
#pragma omp single
    {
      nthreads = omp_get_num_threads();
      printf("I am mpi process %d and I have %d threads\n\n", my_rank,
             nthreads);
    }
  }

  if (my_rank == 0) {
    // kpoint dataset[16] = {{2, 3}, {5, 4}, {9, 6}, {6, 22}, {4, 7},
    //                      {8, 1}, {7, 2}, {8, 9}, {1, 1}, {0.55, 6},
    //                      {1.11111, 8.432}, {9, 2.7}, {4.3253, 42}, {0.1224,
    //                      0.4635}, {0.124, 0.773}, {0.4336, 0.7456}};
    // len = sizeof(dataset) / sizeof(dataset[0]);

    if (argc == 2) {
      len = atoi(argv[1]);
    } else if (argc == 1) {
      len = 100000000;
    } else {
      return 1;
    }

    kpoint *dataset = generate_dataset(len);

    printf("len: %d\n", len);

    kpoint **dataset_ptrs = malloc(len * sizeof(kpoint *));
    get_dataset_ptrs(dataset, dataset_ptrs, len);

    struct kdnode *root;

    tstart = CPU_TIME;

    // extremes[chosen_dimension] = {min_value, max_value} in each dimension
    // NDIM
    kpoint extremes[NDIM];
    get_dataset_extremes(dataset_ptrs, extremes, len, x_axis);
    get_dataset_extremes(dataset_ptrs, extremes, len, y_axis);

    int chosen_axis = choose_splitting_dimension(extremes);

    double sort_start = CPU_TIME;
    printf("starting pre-sorting..\n\n");
#pragma omp parallel shared(dataset_ptrs) firstprivate(chosen_axis, len)
#pragma omp single nowait
    {
      if (chosen_axis == x_axis) {
        pqsort(dataset_ptrs, 0, len, compare_ge_x_axis, compare_g_x_axis);
      } else {
        pqsort(dataset_ptrs, 0, len, compare_ge_y_axis, compare_ge_y_axis);
      }
    }
    printf("pre-sorting done in %f [s]\n", CPU_TIME - sort_start);

    // during the building of the tree, we want to perform an atomic scatter of
    // the dataset to the mpi processes when we reach the "final_level"
    int final_level = log2(numprocs);

    // a variable that is needed to keep track of which process is still free
    // during the atomic scatter
    int is_proc_free[numprocs];
    for (int i = 0; i < numprocs; i++) {
      is_proc_free[i] = 1;
    }

#pragma omp parallel shared(dataset_ptrs, root, is_proc_free)                  \
    firstprivate(extremes, chosen_axis, len, final_level)
#pragma omp single nowait
    {
      int current_level = 0, is_root_proc = 1;
      root = build_kdtree_until_level(dataset_ptrs, extremes, len, chosen_axis,
                                      current_level, final_level, is_root_proc,
                                      is_proc_free);
      printf("finished build kd_tree until level %d\n\n", final_level);
    }

    printf("mpi process %d has root node is %f,%f\n\n", my_rank,
           root->split.coords[0], root->split.coords[1]);

    free(dataset);
    free(dataset_ptrs);
  } else {
    // receiving the dataset chunk that the root process has sent
    int recv_len;
    int recv_axis;
    kpoint *recv_dataset = malloc(recv_len * sizeof(kpoint));
    kpoint *recv_extremes = malloc(NDIM * sizeof(kpoint));

    recv_dataset_from_root_process(&recv_len, &recv_axis, recv_dataset,
                                   recv_extremes);

    kpoint **recv_dataset_ptrs = malloc(recv_len * sizeof(kpoint *));
    get_dataset_ptrs(recv_dataset, recv_dataset_ptrs, recv_len);

    struct kdnode *chunk_root;

    printf("i am mpi process %d, start building my kd-tree..\n\n", my_rank);

#pragma omp parallel shared(recv_dataset_ptrs, chunk_root)                     \
    firstprivate(recv_extremes, recv_axis, recv_len)
#pragma omp single
    {
      int level = 0;
      // chunk_root = build_kdtree(recv_dataset_ptrs, recv_extremes, recv_len,
      //                           recv_axis, level);
    }

    // printf("i am mpi process %d, my chunk root node is %f,%f\n\n", my_rank,
    //        chunk_root->split.coords[0], chunk_root->split.coords[1]);
    free(recv_dataset);
    free(recv_dataset_ptrs);
  }

  double telapsed = CPU_TIME - tstart;
  printf("elapsed time of mpi process %d: %f\n\n", my_rank, telapsed);

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

#define build_cutoff 128

struct kdnode *build_kdtree(kpoint **dataset_ptrs, kpoint extremes[NDIM],
                            int len, int previous_axis, int level) {
  // elementary condition to finish the recursion
  if (len == 1) {
    struct kdnode *leaf = malloc(sizeof(struct kdnode));
    leaf->axis = previous_axis;
    leaf->split = *dataset_ptrs[0];
    leaf->left = NULL;
    leaf->right = NULL;

    return leaf;
  }

  // choose the splitting dimension using the previous extremes of the dataset
  int chosen_axis = choose_splitting_dimension(extremes);

  // sort the array of pointers if there is a change of chosen axis respect to
  // previous chosen axis
#pragma omp taskgroup
  {
    if (chosen_axis != previous_axis) {
      if (chosen_axis == x_axis) {
        pqsort(dataset_ptrs, 0, len, compare_ge_x_axis, compare_g_x_axis);
      } else {
        pqsort(dataset_ptrs, 0, len, compare_ge_y_axis, compare_g_y_axis);
      }
    }
  }

  // take the median of ordered dataset
  kpoint *split_point = choose_splitting_point(dataset_ptrs, len, chosen_axis);

  struct kdnode *node = malloc(sizeof(struct kdnode));
  node->axis = chosen_axis;
  node->split = *split_point;

  kpoint **left_points, **right_points;

  int median_idx = ceil(len / 2.0);
  int len_left = median_idx - 1;    // length of the left points
  int len_right = len - median_idx; // length of the right points

  if (len_left != 0) {
    left_points = &dataset_ptrs[0]; // starting pointer of left_points

    // min value of chosen axis for left points
    extremes[chosen_axis].coords[0] = dataset_ptrs[0]->coords[chosen_axis];
    // max value of chosen axis for left points
    extremes[chosen_axis].coords[1] =
        dataset_ptrs[len_left - 1]->coords[chosen_axis];

#pragma omp task shared(left_points)                                           \
    firstprivate(extremes, len_left, chosen_axis,                              \
                 level) if (len_left >= build_cutoff) mergeable untied
    node->left =
        build_kdtree(left_points, extremes, len_left, chosen_axis, level + 1);
  }

  right_points = &dataset_ptrs[median_idx]; // starting pointer of right_points

  // min value of chosen axis for right points
  extremes[chosen_axis].coords[0] =
      dataset_ptrs[median_idx]->coords[chosen_axis];
  // max value of chosen axis for right points
  extremes[chosen_axis].coords[1] = dataset_ptrs[len - 1]->coords[chosen_axis];

#pragma omp task shared(right_points)                                          \
    firstprivate(extremes, len_right, chosen_axis,                             \
                 level) if (len_right >= build_cutoff) mergeable untied
  node->right =
      build_kdtree(right_points, extremes, len_right, chosen_axis, level + 1);

#pragma omp taskwait
  return node;
}

struct kdnode *build_kdtree_until_level(kpoint **dataset_ptrs,
                                        kpoint extremes[NDIM], int len,
                                        int previous_axis, int current_level,
                                        int final_level, int is_root_proc,
                                        int *is_proc_free) {
  // condition to break the recursion and scatter to other processes
  if ((current_level == final_level) && (is_root_proc != 1)) {
    send_dataset_to_free_process(len, dataset_ptrs, extremes, previous_axis,
                                 final_level, is_proc_free);

    return NULL;
  }

  if (len == 1) {
    struct kdnode *leaf = malloc(sizeof(struct kdnode));
    leaf->axis = previous_axis;
    leaf->split = *dataset_ptrs[0];
    leaf->left = NULL;
    leaf->right = NULL;

    return leaf;
  }

  int chosen_axis = choose_splitting_dimension(extremes);

#pragma omp parallel shared(dataset_ptrs)                                      \
    firstprivate(chosen_axis, previous_axis, len)
#pragma omp single nowait
  {
    if (chosen_axis != previous_axis) {
      if (chosen_axis == x_axis) {
        pqsort(dataset_ptrs, 0, len, compare_ge_x_axis, compare_g_x_axis);
      } else {
        pqsort(dataset_ptrs, 0, len, compare_ge_y_axis, compare_g_y_axis);
      }
    }
  }

  kpoint *split_point = choose_splitting_point(dataset_ptrs, len, chosen_axis);

  struct kdnode *node = malloc(sizeof(struct kdnode));
  node->axis = chosen_axis;
  node->split = *split_point;

  kpoint **left_points, **right_points;

  int median_idx = ceil(len / 2.0) - 1;
  int len_left = median_idx;
  int len_right = len - (median_idx + 1);

  if (len_left != 0) {
    left_points = &dataset_ptrs[0];

    extremes[chosen_axis].coords[0] = dataset_ptrs[0]->coords[chosen_axis];
    extremes[chosen_axis].coords[1] =
        dataset_ptrs[len_left - 1]->coords[chosen_axis];

#pragma omp task shared(left_points, is_root_proc, is_proc_free)               \
    firstprivate(extremes, len_left, chosen_axis, current_level,               \
                 final_level) if (len_left >= build_cutoff) mergeable untied
    node->left = build_kdtree_until_level(
        left_points, extremes, len_left, chosen_axis, current_level + 1,
        final_level, is_root_proc + 0, is_proc_free);
  }

  right_points = &dataset_ptrs[median_idx];

  extremes[chosen_axis].coords[0] =
      dataset_ptrs[median_idx]->coords[chosen_axis];
  extremes[chosen_axis].coords[1] = dataset_ptrs[len - 1]->coords[chosen_axis];

#pragma omp task shared(right_points, is_root_proc, is_proc_free)              \
    firstprivate(extremes, len_right, chosen_axis, current_level,              \
                 final_level) if (len_right >= build_cutoff) mergeable untied
  node->right = build_kdtree_until_level(
      right_points, extremes, len_right, chosen_axis, current_level + 1,
      final_level, is_root_proc + 1, is_proc_free);

#pragma omp taskwait
  return node;
}

kpoint *choose_splitting_point(kpoint **ordered_dataset, int len,
                               int chosen_axis) {
  int median_idx = ceil(len / 2.0) - 1;
  return ordered_dataset[median_idx];
}

int choose_splitting_dimension(kpoint extremes[NDIM]) {
  float_t x_extent =
      extremes[x_axis].coords[1] -
      extremes[x_axis].coords[0]; // max value - min value of x axis
  float_t y_extent =
      extremes[y_axis].coords[1] -
      extremes[y_axis].coords[0]; // max value - min value of y axis

  if (x_extent > y_extent) {
    return x_axis;
  }
  return y_axis;
}

void get_dataset_extremes(kpoint **dataset_ptrs, kpoint extremes[NDIM], int len,
                          int chosen_axis) {
  float_t max_value = dataset_ptrs[0]->coords[chosen_axis];
  float_t min_value = max_value;
#pragma omp parallel for reduction(max                                         \
                                   : max_value) reduction(min                  \
                                                          : min_value)         \
    schedule(static) proc_bind(close)
  for (int i = 1; i < len; i++) {
    max_value = max_value > dataset_ptrs[i]->coords[chosen_axis]
                    ? max_value
                    : dataset_ptrs[i]->coords[chosen_axis];
    min_value = min_value < dataset_ptrs[i]->coords[chosen_axis]
                    ? min_value
                    : dataset_ptrs[i]->coords[chosen_axis];
  }
  extremes[chosen_axis].coords[0] = min_value;
  extremes[chosen_axis].coords[1] = max_value;
}

void send_dataset_to_free_process(int dataset_len, kpoint **dataset_ptrs,
                                  kpoint extremes[NDIM], int previous_axis,
                                  int final_level, int *is_proc_free) {
  // just an hack to compute a power with integer numbers to retrieve numprocs
  int numprocs = (int)(pow(2, final_level) + 1e-9);

  // we have to perform an atomic send inside a thread, searching for the first
  // free process
#pragma omp critical
  for (int mpi_process = 1; mpi_process < numprocs; mpi_process++) {
    if (is_proc_free[mpi_process] == 1) {
      is_proc_free[mpi_process] = 0;
      MPI_Send(&dataset_len, 1, MPI_INT, mpi_process, 0, MPI_COMM_WORLD);
      printf("sent chunk_length from mpi process 0 to mpi process %d, len "
             "%d\n",
             mpi_process, dataset_len);

      kpoint *chunk = malloc(dataset_len * sizeof(kpoint));
      copy_dataset_from_ptrs(chunk, dataset_ptrs, dataset_len);

      MPI_Send(chunk, dataset_len * sizeof(kpoint), MPI_BYTE, mpi_process, 0,
               MPI_COMM_WORLD);
      MPI_Send(extremes, NDIM * sizeof(kpoint), MPI_BYTE, mpi_process, 0,
               MPI_COMM_WORLD);
      MPI_Send(&previous_axis, 1, MPI_INT, mpi_process, 0, MPI_COMM_WORLD);

      printf("sent chunk from mpi process 0, thread %d, to mpi process %d\n",
             omp_get_thread_num(), mpi_process);

      free(chunk);
      break;
    }
  }
}

void recv_dataset_from_root_process(int *recv_len, int *recv_axis,
                                    kpoint *recv_dataset,
                                    kpoint *recv_extremes) {
  int mpi_root_process = 0;
  MPI_Status status;
  MPI_Recv(recv_len, 1, MPI_INT, mpi_root_process, 0, MPI_COMM_WORLD, &status);

  MPI_Recv(recv_dataset, (*recv_len) * sizeof(kpoint), MPI_BYTE,
           mpi_root_process, 0, MPI_COMM_WORLD, &status);

  MPI_Recv(recv_extremes, NDIM * sizeof(kpoint), MPI_BYTE, mpi_root_process, 0,
           MPI_COMM_WORLD, &status);

  MPI_Recv(recv_axis, 1, MPI_INT, mpi_root_process, 0, MPI_COMM_WORLD, &status);
}

void copy_extremes(kpoint old_extremes[NDIM], kpoint new_extremes[NDIM]) {
  for (int dim = 0; dim < NDIM; dim++) {
    new_extremes[dim].coords[0] = old_extremes[dim].coords[0];
    new_extremes[dim].coords[1] = old_extremes[dim].coords[1];
  }
}

void copy_dataset_from_ptrs(kpoint *new_dataset, kpoint **dataset_ptrs,
                            int len) {
  for (int i = 0; i < len; i++) {
    new_dataset[i].coords[0] = dataset_ptrs[i]->coords[0];
    new_dataset[i].coords[1] = dataset_ptrs[i]->coords[1];
  }
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
      SWAP((void *)&data[start + (a)], (void *)&data[start + (b)],             \
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

inline int compare_g_x_axis(const void *A, const void *B) {
  kpoint **a = (kpoint **)A;
  kpoint **b = (kpoint **)B;

  return ((*a)->coords[x_axis] > (*b)->coords[x_axis]);
}

inline int compare_g_y_axis(const void *A, const void *B) {
  kpoint **a = (kpoint **)A;
  kpoint **b = (kpoint **)B;

  return ((*a)->coords[y_axis] > (*b)->coords[y_axis]);
}

int partitioning(kpoint **data, int start, int end,
                 int (*comparator)(const void *, const void *)) {
  --end;
  void *pivot = &data[end];

  int pointbreak = end - 1;
  for (int i = start; i <= pointbreak; i++)
    if (comparator((void *)&data[i], pivot)) {
      while ((pointbreak > i) && comparator(&data[pointbreak], pivot))
        pointbreak--;
      if (pointbreak > i)
        SWAP((void *)&data[i], (void *)&data[pointbreak--], sizeof(&data[i]));
    }
  pointbreak += !comparator(&data[pointbreak], pivot);
  SWAP((void *)&data[pointbreak], pivot, sizeof(&data[pointbreak]));

  return pointbreak;
}

void insertion_sort(kpoint **data, int start, int end,
                    int (*comparator)(const void *, const void *)) {
  {
    int min_idx = start;
    for (int i = start + 1; i < end; i++)
      if (comparator(&data[min_idx], &data[i]))
        min_idx = i;

    SWAP((void *)&data[start], (void *)&data[min_idx], sizeof(&data[start]));
  }

  for (int head = start + 1, run = start + 1; (run = ++head) < end;) {
    while ((run > 0) && comparator(&data[run - 1], &data[run])) {
      SWAP((void *)&data[run - 1], (void *)&data[run], sizeof(&data[run - 1]));
      --run;
    }
  }
}
#define task_cutoff 64
#define insertion_cutoff task_cutoff / 2

void pqsort(kpoint **data, int start, int end,
            int (*comparator)(const void *, const void *),
            int (*comparator_insort)(const void *, const void *)) {
  int size = end - start;

  switch (size) {
  case 1:
    break;
  case 2: {
    if (comparator(&data[start], &data[end - 1]))
      SWAP((void *)&data[start], (void *)&data[end - 1], sizeof(&data[start]));
  } break;
  case 3: {
    CHECKSWAP(1, 2, comparator);
    CHECKSWAP(0, 2, comparator);
    CHECKSWAP(0, 1, comparator);
  } break;
  default: {
    if (size < insertion_cutoff) {
      insertion_sort(data, start, end, comparator_insort);
    } else {

      int mid = partitioning(data, start, end, comparator);

      int mid_start = mid - start;
      if (mid_start > 0)
#pragma omp task default(none) final(mid_start < task_cutoff) mergeable        \
shared(data) firstprivate(start, mid, comparator, comparator_insort) untied
        pqsort(data, start, mid, comparator, comparator_insort);

      int end_mid = end - (mid + 1);
      if (end_mid)
#pragma omp task default(none) final(end_mid < task_cutoff) mergeable shared(  \
    data) firstprivate(mid, end, comparator, comparator_insort) untied
        pqsort(data, mid + 1, end, comparator, comparator_insort);
    }
  } break;
  }
}
