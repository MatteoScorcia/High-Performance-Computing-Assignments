#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

kpoint *generate_dataset(int len);
struct kdnode *build_kdtree(kpoint *dataset, int len, int current_axis, int level);

int main(int argc, char *argv[]) {
  
	int numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0) {
  int len;

  if (argc == 2 ) {
    len = atoi(argv[1]);  
  } else if (argc == 1){
    len = 80000000;
  } else {
    return 1;
  }

    printf("len: %d\n", len);
    kpoint *dataset = generate_dataset(len);

    //qsort(dataset);
    //build kd_tree(dataset);
  }
  
  free(dataset);
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

// #define build_cutoff 128
//
// struct kdnode *build_kdtree(kpoint *dataset, int len, int current_axis, int level) {
//   if (len == 1) {
//     struct kdnode *leaf = malloc(sizeof(struct kdnode));
//     leaf->axis = current_axis;
//     leaf->split = dataset[0];
//     leaf->left = NULL;
//     leaf->right = NULL;
//
//     return leaf;
//   } else if (len == 0) {
//     return NULL;
//   }
//
//   struct kdnode *node = malloc(sizeof(struct kdnode));
//
//   int chosen_axis = choose_splitting_dimension(dataset_ptrs, len);
//
//   if((level >= 0) && (level <= 3)) {
//     printf("reached level %d\n", level);
//   }
//
//   #pragma omp taskgroup
//   {
//     if (chosen_axis != current_axis) {
//       if(chosen_axis == x_axis) {
//         pqsort(dataset_ptrs, 0, len, compare_ge_x_axis, compare_g_x_axis);
//       } else {
//         pqsort(dataset_ptrs, 0, len, compare_ge_y_axis, compare_g_y_axis);
//       }
//     }   
//   }
//
//   kpoint *split_point = choose_splitting_point(dataset_ptrs, len, chosen_axis);
//
//   kpoint **left_points, **right_points;
//
//   int median_idx = ceil(len / 2.0);
//   int len_left = median_idx - 1;    // length of the left points
//   int len_right = len - median_idx; // length of the right points
//   
//   left_points = &dataset_ptrs[0];       // starting pointer of left_points
//   right_points = &dataset_ptrs[median_idx]; // starting pointer of right_points
//   
//   node->axis = chosen_axis;
//   node->split = *split_point;
//   #pragma omp task shared(left_points) firstprivate(len_left, chosen_axis, level) if(len_left >= build_cutoff) mergeable untied
//     node->left = build_kdtree(left_points, len_left, chosen_axis, level+1);
//   #pragma omp task shared(right_points) firstprivate(len_right, chosen_axis, level) if(len_right >= build_cutoff) mergeable untied
//     node->right = build_kdtree(right_points, len_right, chosen_axis, level+1);
//
//   #pragma omp taskwait
//   return node;
// }
//
// kpoint *choose_splitting_point(kpoint **ordered_dataset, int len,
//                                int chosen_axis) {
//   int median_idx = ceil(len / 2.0) - 1;
//   return ordered_dataset[median_idx];
// }
//
// int choose_splitting_dimension(kpoint **dataset_ptrs, int len) {
//   float_t x_extent = get_dataset_extent(dataset_ptrs, len, x_axis);
//   float_t y_extent = get_dataset_extent(dataset_ptrs, len, y_axis);
//
//   if (x_extent > y_extent) {
//     return x_axis;
//   }
//   return y_axis;
// }
//
// float_t get_dataset_extent(kpoint **dataset, int len, int axis) {
//   float_t max = get_max_value_dataset(dataset, len, axis);
//   float_t min = get_min_value_dataset(dataset, len, axis);
//
//   return (max - min);
// }
//
// float_t get_max_value_dataset(kpoint **dataset, int len, int axis) {
//   float_t max_value = dataset[0]->coords[axis];
//   #pragma omp parallel for reduction(max:max_value) schedule(static) proc_bind(close)
//   for (int i = 1; i < len; i++) {
//     max_value = max_value > dataset[i]->coords[axis] ? max_value : dataset[i]->coords[axis];
//   }
//
//   return max_value;
// }
//
// float_t get_min_value_dataset(kpoint **dataset, int len, int axis) {
//   float_t min_value = dataset[0]->coords[axis];
//   #pragma omp parallel for reduction(min:min_value) schedule(static) proc_bind(close)
//   for (int i = 1; i < len; i++) {
//     min_value = min_value < dataset[i]->coords[axis] ? min_value : dataset[i]->coords[axis];
//   }
//
//   return min_value;
// }
//
// void copy_dataset_ptrs(kpoint **dataset_ptrs, kpoint **new_dataset, int len) {
//   for (int i = 0; i < len; i++) {
//     new_dataset[i] = dataset_ptrs[i];
//   }
// }
//
// void get_dataset_ptrs(kpoint *dataset, kpoint **dataset_ptrs, int len) {
//   for (int i = 0; i < len; i++) {
//     dataset_ptrs[i] = &dataset[i];
//   }
// }
//
// #define SWAP(A, B, SIZE)                                                       \
//   do {                                                                         \
//     int sz = (SIZE);                                                           \
//     char *a = (A);                                                             \
//     char *b = (B);                                                             \
//     do {                                                                       \
//       char _temp = *a;                                                         \
//       *a++ = *b;                                                               \
//       *b++ = _temp;                                                            \
//     } while (--sz);                                                            \
//   } while (0)
//
// #define CHECKSWAP(a, b, comparator)                                            \
//   {                                                                            \
//     if (comparator(&data[start + (a)], &data[start + (b)]))                    \
//       SWAP((void *)&data[start + (a)], (void *)&data[start + (b)],                             \
//            sizeof(&data[start + (a)]));                                        \
//   }
//
// inline int compare_ge_x_axis(const void *A, const void *B) {
//   kpoint **a = (kpoint **)A;
//   kpoint **b = (kpoint **)B;
//
//   return ((*a)->coords[x_axis] >= (*b)->coords[x_axis]);
// }
//
// inline int compare_ge_y_axis(const void *A, const void *B) {
//   kpoint **a = (kpoint **)A;
//   kpoint **b = (kpoint **)B;
//
//   return ((*a)->coords[y_axis] >= (*b)->coords[y_axis]);
// }
//
// inline int compare_g_x_axis(const void *A, const void *B) {
//   kpoint **a = (kpoint **)A;
//   kpoint **b = (kpoint **)B;
//
//   return ((*a)->coords[x_axis] > (*b)->coords[x_axis]);
// }
//
// inline int compare_g_y_axis(const void *A, const void *B) {
//   kpoint **a = (kpoint **)A;
//   kpoint **b = (kpoint **)B;
//
//   return ((*a)->coords[y_axis] > (*b)->coords[y_axis]);
// }
//
// int partitioning(kpoint **data, int start, int end,
//                  int (*comparator)(const void *, const void *)) {
//   --end;
//   void *pivot = &data[end];
//
//   int pointbreak = end - 1;
//   for (int i = start; i <= pointbreak; i++)
//     if (comparator((void *)&data[i], pivot)) {
//       while ((pointbreak > i) && comparator(&data[pointbreak], pivot))
//         pointbreak--;
//       if (pointbreak > i)
//         SWAP((void *)&data[i], (void *)&data[pointbreak--], sizeof(&data[i]));
//     }
//   pointbreak += !comparator(&data[pointbreak], pivot);
//   SWAP((void *)&data[pointbreak], pivot, sizeof(&data[pointbreak]));
//
//   return pointbreak;
// }
//
// void insertion_sort(kpoint **data, int start, int end,
//                     int (*comparator)(const void *, const void *)) {
//   {
//     int min_idx = start;
//     for (int i = start + 1; i < end; i++)
//       if (comparator(&data[min_idx], &data[i]))
//         min_idx = i;
//
//     SWAP((void *)&data[start], (void *)&data[min_idx], sizeof(&data[start]));
//   }
//
//   for (int head = start + 1, run = start + 1; (run = ++head) < end;) {
//     while ((run > 0) && comparator(&data[run - 1], &data[run])) {
//       SWAP((void *)&data[run - 1], (void *)&data[run], sizeof(&data[run - 1]));
//       --run;
//     }
//   }
// }
// #define task_cutoff 64
// #define insertion_cutoff task_cutoff / 2
//
// void pqsort(kpoint **data, int start, int end,
//             int (*comparator)(const void *, const void *), int(*comparator_insort)(const void *, const void *)) {
//   int size = end - start;
//
//   switch (size) {
//   case 1:
//     break;
//   case 2: {
//     if (comparator(&data[start], &data[end - 1]))
//       SWAP((void *)&data[start], (void *)&data[end - 1], sizeof(&data[start]));
//   } break;
//   case 3: {
//     CHECKSWAP(1, 2, comparator);
//     CHECKSWAP(0, 2, comparator);
//     CHECKSWAP(0, 1, comparator);
//   } break;
//   default: {
//     if (size < insertion_cutoff) {
//       insertion_sort(data, start, end, comparator_insort);
//     } else {
//
//       int mid = partitioning(data, start, end, comparator);
//
//       int mid_start = mid - start;
//       if (mid_start > 0)
// #pragma omp task default(none) final(mid_start < task_cutoff) mergeable        \
// shared(data) firstprivate(start, mid, comparator, comparator_insort) untied
//         pqsort(data, start, mid, comparator, comparator_insort);
//
//       int end_mid = end - (mid + 1);
//       if (end_mid)
// #pragma omp task default(none) final(end_mid < task_cutoff) mergeable shared(  \
//     data) firstprivate(mid, end, comparator, comparator_insort) untied
//         pqsort(data, mid + 1, end, comparator, comparator_insort);
//     }
//   } break;
//   }
// }
