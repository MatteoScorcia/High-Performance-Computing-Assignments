#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <mpi.h>

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
#define MPI_FLOAT_T MPI_FLOAT
#else
#define float_t double
#define MPI_FLOAT_T MPI_DOUBLE
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
            int (*comparator)(const void *, const void *), int(*comparator_insort)(const void *, const void *));
int compare_ge_x_axis(const void *a, const void *b);
int compare_ge_y_axis(const void *a, const void *b);
int compare_g_x_axis(const void *a, const void *b);
int compare_g_y_axis(const void *a, const void *b);

kpoint *generate_dataset(int len);
void get_dataset_ptrs(kpoint *dataset, kpoint **dataset_ptrs, int len);
void copy_dataset_ptrs(kpoint **dataset_ptrs, kpoint **new_dataset, int len);
void copy_dataset_from_ptrs(kpoint *new_dataset, kpoint **dataset_ptrs, int len);

// kd-tree build functions
struct kdnode *build_kdtree(kpoint **dataset_ptrs, float_t extremes[NDIM][2], int len, int axis, int level);
struct kdnode *build_kdtree_until_level_then_scatter(kpoint **dataset_ptrs, float_t extremes[NDIM][2], int len, int axis, int level, int final_level, int counter);
int choose_splitting_dimension(float_t extremes[NDIM][2]);
kpoint *choose_splitting_point(kpoint **dataset_ptrs, int len, int chosen_axis);
void get_dataset_extremes(kpoint **dataset, float_t extremes[NDIM][2], int len, int axis);
void copy_extremes(float_t old_extremes[NDIM][2], float_t new_extremes[NDIM][2]);

int main(int argc, char *argv[]) {

  struct timespec ts;
  int len;

	int numprocs;
  int *provided = malloc(sizeof(int *));
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, provided);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0) {
    kpoint dataset[11] = {{2, 3}, {5, 4}, {9, 6}, {6, 22}, {4, 7},
                         {8, 1}, {7, 2}, {8, 9}, {1, 1}, {0.55, 6}, {1.11111, 8.432}};
    len = sizeof(dataset) / sizeof(dataset[0]);

    // if (argc == 2 ) {
    //   len = atoi(argv[1]); 
    // } else if (argc == 1){
    //   len = 100000000;
    // } else {
    //   return 1;
    // }
    //
    // kpoint *dataset = generate_dataset(len);

    printf("len: %d\n", len);

    kpoint **dataset_ptrs = malloc(len * sizeof(kpoint *));
    get_dataset_ptrs(dataset, dataset_ptrs, len);

    int nthreads;
    struct kdnode *root;

    double tstart = CPU_TIME;
    #pragma omp parallel shared(nthreads)
    {
      #pragma omp single
      {
        nthreads = omp_get_num_threads();
        printf("num threads: %d\n", nthreads);
      }
    }

    printf("choosing splitting dimension..\n");
    
    float_t extremes[NDIM][2]; //min_value (index 0) and max value (index 1) in each dimension NDIM

    get_dataset_extremes(dataset_ptrs, extremes, len, x_axis);
    get_dataset_extremes(dataset_ptrs, extremes, len, y_axis);
   
    int chosen_axis = choose_splitting_dimension(extremes);

    printf("starting pre-sorting..\n");
    #pragma omp parallel shared(dataset_ptrs) firstprivate(chosen_axis, len)
    {
      #pragma omp single nowait
      {
        if(chosen_axis == x_axis) {
          pqsort(dataset_ptrs, 0, len, compare_ge_x_axis, compare_g_x_axis);
        } else {
          pqsort(dataset_ptrs, 0, len, compare_ge_y_axis, compare_ge_y_axis);
        }
      }
    }
    printf("pre-sorting done!\n");

    #pragma omp parallel shared(dataset_ptrs) firstprivate(extremes, chosen_axis, len) 
    {
      #pragma omp single nowait
      {
        int current_level = 0, counter = 0;
        int final_level = log2(numprocs);
        root = build_kdtree_until_level_then_scatter(dataset_ptrs, extremes, len, chosen_axis, current_level, final_level, counter);
        printf("finished build kd_tree until level %d\n", final_level);
      }
    }

    double telapsed = CPU_TIME - tstart;

    printf("elapsed time: %f\n", telapsed);

    // free(dataset);
    free(dataset_ptrs);
  } 

  int recv_len;
  double len2;
  MPI_Status status;
  MPI_Recv(&recv_len, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
  printf("i am processor %d, received len %d\n", my_rank, recv_len);
  MPI_Recv(&len2, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
  printf("i am processor %d, received len2 %f\n", my_rank, len2);

  // double prova;
  // MPI_Recv(&prova, 1, MPI_FLOAT_T, 0, 0, MPI_COMM_WORLD, &status);
  // printf("received prova %f\n", prova);
  // kpoint *recv_dataset = malloc(len * sizeof(kpoint));
  // MPI_Recv(&recv_dataset, len*2, MPI_FLOAT_T, 0, 0, MPI_COMM_WORLD, &status);
  
  // printf("i am processor %d, first point received: (%f,%f)\n", my_rank, recv_dataset[0].coords[0], recv_dataset[0].coords[1]);
   
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

struct kdnode *build_kdtree(kpoint **dataset_ptrs, float_t extremes[NDIM][2], int len, int previous_axis, int level) {
  if (len == 1) {
    struct kdnode *leaf = malloc(sizeof(struct kdnode));
    leaf->axis = previous_axis;
    leaf->split = *dataset_ptrs[0];
    leaf->left = NULL;
    leaf->right = NULL;

    return leaf;
  }

  struct kdnode *node = malloc(sizeof(struct kdnode));

  int chosen_axis = choose_splitting_dimension(extremes);

  if((level >= 0) && (level <= 3)) {
    printf("reached level %d, len is %d\n", level, len);
  }

  // printf("y axis max: %f, min: %f\n", extremes[y_axis][1], extremes[y_axis][0]);
  // printf("x axis max: %f, min: %f\n", extremes[x_axis][1], extremes[x_axis][0]);

  #pragma omp taskgroup
  {
    if (chosen_axis != previous_axis) {
      if(chosen_axis == x_axis) {
        pqsort(dataset_ptrs, 0, len, compare_ge_x_axis, compare_g_x_axis);
      } else {
        pqsort(dataset_ptrs, 0, len, compare_ge_y_axis, compare_g_y_axis);
      }
    }   
  }

  kpoint *split_point = choose_splitting_point(dataset_ptrs, len, chosen_axis);
  node->axis = chosen_axis;
  node->split = *split_point;

  kpoint **left_points, **right_points;

  int median_idx = ceil(len / 2.0);
  int len_left = median_idx - 1;    // length of the left points
  int len_right = len - median_idx; // length of the right points
  
  if (len_left != 0) {
    left_points = &dataset_ptrs[0];       // starting pointer of left_points

    extremes[chosen_axis][0] = dataset_ptrs[0]->coords[chosen_axis]; //min value of chosen axis for left points
    extremes[chosen_axis][1] = dataset_ptrs[len_left - 1]->coords[chosen_axis]; //max value of chosen axis for left points

    #pragma omp task shared(left_points) firstprivate(extremes, len_left, chosen_axis, level) if(len_left >= build_cutoff) mergeable untied
      node->left = build_kdtree(left_points, extremes, len_left, chosen_axis, level+1);
  }

  right_points = &dataset_ptrs[median_idx]; // starting pointer of right_points
  
  extremes[chosen_axis][0] = dataset_ptrs[median_idx]->coords[chosen_axis]; //min value of chosen axis for right points
  extremes[chosen_axis][1] = dataset_ptrs[len - 1]->coords[chosen_axis]; //max value of chosen axis for right points

  #pragma omp task shared(right_points) firstprivate(extremes, len_right, chosen_axis, level) if(len_right >= build_cutoff) mergeable untied
    node->right = build_kdtree(right_points, extremes, len_right, chosen_axis, level+1);
  
  #pragma omp taskwait
  return node;
}

struct kdnode *build_kdtree_until_level_then_scatter(kpoint **dataset_ptrs, float_t extremes[NDIM][2], int len, int previous_axis, int current_level, int final_level, int counter) {
  if (current_level == final_level) {
    #pragma omp critical 
    {
      printf("sending to mpi process %d, dataset chunk\n", counter);
      MPI_Send(&len, 1, MPI_INT, counter, 0, MPI_COMM_WORLD);
      kpoint *chunk = malloc(len * sizeof(kpoint));
      copy_dataset_from_ptrs(chunk, dataset_ptrs, len);
      double len2 = chunk[0].coords[0];
      MPI_Send(&len2, 1, MPI_DOUBLE, counter, 0, MPI_COMM_WORLD);
      // kpoint *chunk = malloc(len * sizeof(kpoint));
      // copy_dataset_from_ptrs(chunk, dataset_ptrs, len);
      // printf("first kpoint sent is (%f,%f)\n", chunk[0].coords[0], chunk[0].coords[1]);
      // MPI_Send(&chunk, 2, MPI_FLOAT_T, counter, 0, MPI_COMM_WORLD);
      // free(chunk);
      // double prova = dataset_ptrs[0]->coords[0];
      // MPI_Send(&prova, 1, MPI_FLOAT_T, counter, 0, MPI_COMM_WORLD);
      counter++;
    }
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

  struct kdnode *node = malloc(sizeof(struct kdnode));

  int chosen_axis = choose_splitting_dimension(extremes);

  if((current_level >= 0) && (current_level <= 3)) {
    printf("reached level %d, len is %d\n", current_level, len);
  }

  #pragma omp taskgroup
  {
    if (chosen_axis != previous_axis) {
      if(chosen_axis == x_axis) {
        pqsort(dataset_ptrs, 0, len, compare_ge_x_axis, compare_g_x_axis);
      } else {
        pqsort(dataset_ptrs, 0, len, compare_ge_y_axis, compare_g_y_axis);
      }
    }   
  }

  kpoint *split_point = choose_splitting_point(dataset_ptrs, len, chosen_axis);
  node->axis = chosen_axis;
  node->split = *split_point;

  kpoint **left_points, **right_points;

  int median_idx = ceil(len / 2.0);
  int len_left = median_idx - 1;    // length of the left points
  int len_right = len - median_idx; // length of the right points
  
  if (len_left != 0) {
    left_points = &dataset_ptrs[0];       // starting pointer of left_points

    extremes[chosen_axis][0] = dataset_ptrs[0]->coords[chosen_axis]; //min value of chosen axis for left points
    extremes[chosen_axis][1] = dataset_ptrs[len_left - 1]->coords[chosen_axis]; //max value of chosen axis for left points

    #pragma omp task shared(left_points) firstprivate(extremes, len_left, chosen_axis, current_level, final_level, counter) if(len_left >= build_cutoff) mergeable untied
      node->left = build_kdtree_until_level_then_scatter(left_points, extremes, len_left, chosen_axis, current_level+1, final_level, counter);
  }

  right_points = &dataset_ptrs[median_idx]; // starting pointer of right_points
  
  extremes[chosen_axis][0] = dataset_ptrs[median_idx]->coords[chosen_axis]; //min value of chosen axis for right points
  extremes[chosen_axis][1] = dataset_ptrs[len - 1]->coords[chosen_axis]; //max value of chosen axis for right points

  #pragma omp task shared(right_points) firstprivate(extremes, len_right, chosen_axis, current_level, final_level, counter) if(len_right >= build_cutoff) mergeable untied
    node->right = build_kdtree_until_level_then_scatter(right_points, extremes, len_right, chosen_axis, current_level+1, final_level, counter);
  
  #pragma omp taskwait
  return node;
}

kpoint *choose_splitting_point(kpoint **ordered_dataset, int len,
                               int chosen_axis) {
  int median_idx = ceil(len / 2.0) - 1;
  return ordered_dataset[median_idx];
}

int choose_splitting_dimension(float_t extremes[NDIM][2]) {
  float_t x_extent = extremes[x_axis][1] - extremes[x_axis][0]; //max value - min value
  float_t y_extent = extremes[y_axis][1] - extremes[y_axis][0]; //max value - min value

  if (x_extent > y_extent) {
    return x_axis;
  }
  return y_axis;
}

void get_dataset_extremes(kpoint **dataset, float_t extremes[NDIM][2], int len, int axis) {
  float_t max_value = dataset[0]->coords[axis];
  float_t min_value = max_value;
  #pragma omp parallel for reduction(max:max_value) reduction(min:min_value) schedule(static) proc_bind(close)
  for (int i = 1; i < len; i++) {
    max_value = max_value > dataset[i]->coords[axis] ? max_value : dataset[i]->coords[axis];
    min_value = min_value < dataset[i]->coords[axis] ? min_value : dataset[i]->coords[axis];
  }
  extremes[axis][0] = min_value;
  extremes[axis][1] = max_value;
}

void copy_extremes(float_t old_extremes[NDIM][2], float_t new_extremes[NDIM][2]) {
  for (int dim = 0; dim < NDIM; dim++) {
    new_extremes[dim][0] = old_extremes[dim][0];
    new_extremes[dim][1] = old_extremes[dim][1];  
  }
}

void copy_dataset_from_ptrs(kpoint *new_dataset, kpoint **dataset_ptrs, int len) {
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
      SWAP((void *)&data[start + (a)], (void *)&data[start + (b)],                             \
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
            int (*comparator)(const void *, const void *), int(*comparator_insort)(const void *, const void *)) {
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
