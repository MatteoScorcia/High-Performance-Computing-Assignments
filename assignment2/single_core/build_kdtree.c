#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

// kd-tree build functions
struct kdnode *build_kdtree(kpoint **dataset_ptr, int len, int ndim, int axis);
int choose_splitting_dimension(kpoint **dataset_ptrs, int len);
kpoint *choose_splitting_point(kpoint **dataset_ptrs, int len, int chosen_axis);
float_t get_dataset_extent(kpoint **arr, int len, int axis);
float_t get_max_value_dataset(kpoint **arr, int len, int axis);
float_t get_min_value_dataset(kpoint **arr, int len, int axis);
int cmpfunc_x_axis(const void *a, const void *b);
int cmpfunc_y_axis(const void *a, const void *b);

// utility functions
kpoint *generate_dataset(int len);
void get_dataset_ptrs(kpoint *dataset, kpoint **dataset_ptrs, int len);
void copy_dataset_ptrs(kpoint **dataset_ptrs, kpoint **new_arr, int len);
void print_dataset(kpoint *dataset, int len);
void print_dataset_ptr(kpoint **dataset_ptr, int len);

int main(int argc, char *argv[]) {

  struct  timespec ts;
  int len = 1000000;
  kpoint *dataset = generate_dataset(len);

  kpoint **dataset_ptrs = malloc(len * sizeof(kpoint *));
  get_dataset_ptrs(dataset, dataset_ptrs, len);

  printf("extent of x components of dataset: %f\n",
         get_dataset_extent(dataset_ptrs, len, x_axis));

  printf("extent of y components of dataset: %f\n",
         get_dataset_extent(dataset_ptrs, len, y_axis));

  int chosen_axis = choose_splitting_dimension(dataset_ptrs, len);

  printf("choose splitting dimension (0=x, 1=y) -> %d\n", chosen_axis);

  if (chosen_axis == x_axis) {
    qsort(dataset_ptrs, len, sizeof(kpoint *), cmpfunc_x_axis);
  } else if (chosen_axis == y_axis) {
    qsort(dataset_ptrs, len, sizeof(kpoint *), cmpfunc_y_axis);
  }

  kpoint *split_point = choose_splitting_point(dataset_ptrs, len, chosen_axis);
  printf("chosen splitting point: (%f,%f)\n", (*split_point).coords[x_axis],
         (*split_point).coords[y_axis]);

  printf("\n--- end of testing ---\n\n");

  get_dataset_ptrs(dataset, dataset_ptrs, len);

  double tstart = CPU_TIME;
  build_kdtree(dataset_ptrs, len, 2, 0);
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

struct kdnode *build_kdtree(kpoint **dataset_ptrs, int len, int ndim,
                            int axis) {
  printf("len: %d\n", len);
  if (len == 1) {
    struct kdnode *leaf = malloc(sizeof(struct kdnode));
    leaf->axis = axis;
    leaf->split = *dataset_ptrs[0];
    leaf->left = NULL;
    leaf->right = NULL;

    printf("chosen axis: %d\n", axis);
    printf("chosen splitting point: (%f,%f)\n\n", dataset_ptrs[0]->coords[0],
           dataset_ptrs[0]->coords[1]);
    return leaf;
  } else if (len == 0) {
    printf("\n");
    return NULL;
  }

  printf("input dataset: \n");
  print_dataset_ptr(dataset_ptrs, len);
  printf("\n");

  struct kdnode *node = malloc(sizeof(struct kdnode));

  int chosen_axis = choose_splitting_dimension(dataset_ptrs, len);

  if (chosen_axis == x_axis) {
    qsort(dataset_ptrs, len, sizeof(kpoint *), cmpfunc_x_axis);
  } else {
    qsort(dataset_ptrs, len, sizeof(kpoint *), cmpfunc_y_axis);
  }

  printf("sorted dataset\n");
  print_dataset_ptr(dataset_ptrs, len);
  printf("\n");

  kpoint *split_point = choose_splitting_point(dataset_ptrs, len, chosen_axis);

  printf("chosen axis: %d\n", chosen_axis);
  printf("chosen splitting point: (%f,%f)\n", (*split_point).coords[0],
         (*split_point).coords[1]);

  kpoint **left_points, **right_points;

  int median = ceil(len / 2.0);
  int len_left = median - 1;    // length of the left points
  int len_right = len - median; // length of the right points

  printf("median: %d, len_left: %d, len_right: %d\n\n", median, len_left,
         len_right);

  left_points = &dataset_ptrs[0];       // starting pointer of left_points
  right_points = &dataset_ptrs[median]; // starting pointer of right_points

  node->axis = chosen_axis;
  node->split = *split_point;
  node->left = build_kdtree(left_points, len_left, ndim, chosen_axis);
  node->right = build_kdtree(right_points, len_right, ndim, chosen_axis);

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

int cmpfunc_x_axis(const void *a, const void *b) {
  kpoint **ptr_a = (kpoint **)a;
  kpoint **ptr_b = (kpoint **)b;

  return (((**ptr_a).coords[x_axis] > (**ptr_b).coords[x_axis]) -
          ((**ptr_a).coords[x_axis] < (**ptr_b).coords[x_axis]));
}

int cmpfunc_y_axis(const void *a, const void *b) {
  kpoint **ptr_a = (kpoint **)a;
  kpoint **ptr_b = (kpoint **)b;

  return (((**ptr_a).coords[y_axis] > (**ptr_b).coords[y_axis]) -
          ((**ptr_a).coords[y_axis] < (**ptr_b).coords[y_axis]));
}

float_t get_max_value_dataset(kpoint **arr, int len, int axis) {
  float_t temp = (*arr[0]).coords[axis];
  for (int i = 1; i < len; i++) {
    if ((*arr[i]).coords[axis] > temp) {
      temp = (*arr[i]).coords[axis];
    }
  }

  return temp;
}

float_t get_min_value_dataset(kpoint **arr, int len, int axis) {
  float_t temp = (*arr[0]).coords[axis];
  for (int i = 1; i < len; i++) {
    if ((*arr[i]).coords[axis] < temp) {
      temp = (*arr[i]).coords[axis];
    }
  }

  return temp;
}

// utility func
void copy_dataset_ptrs(kpoint **dataset_ptrs, kpoint **new_arr, int len) {
  for (int i = 0; i < len; i++) {
    new_arr[i] = dataset_ptrs[i];
  }
}

// utility func
void get_dataset_ptrs(kpoint *dataset, kpoint **dataset_ptrs, int len) {
  for (int i = 0; i < len; i++) {
    dataset_ptrs[i] = &dataset[i];
  }
}

// utility func
void print_dataset_ptr(kpoint **dataset_ptr, int len) {
  for (int i = 0; i < len; i++) {
    printf("dataset[%d] -> (%f,%f)\n", i, dataset_ptr[i]->coords[0],
           dataset_ptr[i]->coords[1]);
  }
}

// utility func
void print_dataset(kpoint *dataset, int len) {
  for (int i = 0; i < len; i++) {
    printf("dataset[%d] -> (%f,%f)\n", i, dataset[i].coords[0],
           dataset[i].coords[1]);
  }
  printf("\n");
}

// useless function, here just to understand how to
// write a cmpfunc that uses pointers
int cmpfunc_ptr(const void *a, const void *b) {
  float_t **ptr_a = (float_t **)a;
  float_t **ptr_b = (float_t **)b;

  return ((**ptr_a > **ptr_b) - (**ptr_a < **ptr_b));
}
