#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

// useless functions
float_t nlogn_median(float_t *arr, int len);
void get_axis_coord(kpoint *dataset, float_t *arr, int len, int axis);
int cmpfunc(const void *a, const void *b);

// useless functions
void get_axis_coord_ptr(kpoint *dataset, float_t **arr, int len, int axis);
float_t get_array_variance(float_t **arr, int len);
int choose_splitting_dimension_variance(float_t **x_dataset_ptr,
                                        float_t **y_dataset_ptr, int len);
int cmpfunc_ptr(const void *a, const void *b);

// actual functions used to work with dataset pointers
void get_dataset_ptrs(kpoint *dataset, kpoint **dataset_ptrs, int len);
float_t get_dataset_extent(kpoint **arr, int len, int axis);
int choose_splitting_dimension(kpoint **x_ordered_dataset,
                               kpoint **y_ordered_dataset, int len);
kpoint *choose_splitting_point(kpoint **x_ordered_dataset,
                               kpoint **y_ordered_dataset, int len,
                               int chosen_axis);
struct kdnode *build_kdtree(kpoint *dataset, int len, int ndim, int axis);
int cmpfunc_x_axis(const void *a, const void *b);
int cmpfunc_y_axis(const void *a, const void *b);

void print_dataset(kpoint **dataset, int len);

int main(int argc, char *argv[]) {

  kpoint dataset[9] = {{2, 3}, {5, 4}, {9, 6}, {6, 22}, {4, 7},
                       {8, 1}, {7, 2}, {8, 9}, {1, 1}};
  int len = sizeof(dataset) / sizeof(dataset[0]);

  kpoint *x_ordered_dataset[len];
  kpoint *y_ordered_dataset[len];

  get_dataset_ptrs(dataset, x_ordered_dataset, len);
  get_dataset_ptrs(dataset, y_ordered_dataset, len);

  qsort(x_ordered_dataset, len, sizeof(kpoint *), cmpfunc_x_axis);
  qsort(y_ordered_dataset, len, sizeof(kpoint *), cmpfunc_y_axis);

  printf("last x of the ordered dataset by x: %f\n",
         (*x_ordered_dataset[len - 1]).coords[x_axis]);

  printf("last y of the ordered dataset by x: %f\n",
         (*y_ordered_dataset[len - 1]).coords[y_axis]);

  printf("extent of x components of dataset: %f\n",
         get_dataset_extent(x_ordered_dataset, len, x_axis));

  printf("extent of y components of dataset: %f\n",
         get_dataset_extent(y_ordered_dataset, len, y_axis));

  int chosen_axis =
      choose_splitting_dimension(x_ordered_dataset, y_ordered_dataset, len);

  printf("choose splitting dimension (0=x, 1=y) -> %d\n", chosen_axis);

  kpoint *split_point = choose_splitting_point(
      x_ordered_dataset, y_ordered_dataset, len, chosen_axis);

  printf("chosen splitting point: (%f,%f)\n", (*split_point).coords[x_axis],
         (*split_point).coords[y_axis]);

  printf("\n--- end of testing ---\n\n");

  build_kdtree(dataset, len, 2, 0);

  return 0;
}

struct kdnode *build_kdtree(kpoint *dataset, int len, int ndim, int axis) {
  printf("len: %d\n", len);
  if (len == 1) {
    struct kdnode *leaf = malloc(sizeof(struct kdnode));
    leaf->axis = axis;
    leaf->split = dataset[0];
    leaf->left = NULL;
    leaf->right = NULL;

    printf("chosen axis: %d\n", axis);
    printf("chosen splitting point: (%f,%f)\n\n", dataset[0].coords[0],
           dataset[0].coords[1]);
    return leaf;
  }else if (len == 0) {
    printf("\n");
    return NULL;
  }

  struct kdnode *node = malloc(sizeof(struct kdnode));

  kpoint *x_ordered_dataset[len];
  kpoint *y_ordered_dataset[len];

  get_dataset_ptrs(dataset, x_ordered_dataset, len);
  get_dataset_ptrs(dataset, y_ordered_dataset, len);

  qsort(x_ordered_dataset, len, sizeof(kpoint *), cmpfunc_x_axis);
  qsort(y_ordered_dataset, len, sizeof(kpoint *), cmpfunc_y_axis);

  printf("sorting done\n");

  int chosen_axis =
      choose_splitting_dimension(x_ordered_dataset, y_ordered_dataset, len);
  kpoint *split_point = choose_splitting_point(
      x_ordered_dataset, y_ordered_dataset, len, chosen_axis);

  printf("chosen axis: %d\n", chosen_axis);
  printf("chosen splitting point: (%f,%f)\n", (*split_point).coords[0],
         (*split_point).coords[1]);

  kpoint *left_points, *right_points;

  int median = ceil(len / 2.0); 
  int len_left = median - 1; // length of the left points 
  int len_right = len - median; // length of the right points 

  printf("median: %d, len_left: %d, len_right: %d\n\n", median, len_left,
         len_right);

  if (chosen_axis == x_axis) {
    left_points = x_ordered_dataset[0]; // sorted points from 0 to median - 1
    right_points =
        x_ordered_dataset[median + 1]; // sorted points from median + 1 to len
        print_dataset(x_ordered_dataset, len);
  } else {
    left_points = y_ordered_dataset[0]; // sorted points from 0 to len/2 - 1
    right_points =
        y_ordered_dataset[median + 1]; // sorted points from len/2 + 1 to len
        print_dataset(y_ordered_dataset, len);
  }

  node->axis = chosen_axis;
  node->split = *split_point;
  node->left = build_kdtree(left_points, len_left, ndim, chosen_axis);
  node->right = build_kdtree(right_points, len_right, ndim, chosen_axis);

  return node;
}

//utility func
void print_dataset(kpoint **dataset, int len) {
  for (int i = 0; i < len; i++) {
    printf("dataset[%d] -> (%f,%f)\n", i, dataset[i]->coords[0], dataset[i]->coords[1]);
  }
  printf("\n");
}

kpoint *choose_splitting_point(kpoint **x_ordered_dataset,
                               kpoint **y_ordered_dataset, int len,
                               int chosen_axis) {
  int median_idx = ceil(len / 2.0) - 1;
  if (chosen_axis == 0) {
    return x_ordered_dataset[median_idx];
  }
  return y_ordered_dataset[median_idx];
}

int choose_splitting_dimension(kpoint **x_ordered_dataset,
                               kpoint **y_ordered_dataset, int len) {
  float_t x_extent = get_dataset_extent(x_ordered_dataset, len, x_axis);
  float_t y_extent = get_dataset_extent(y_ordered_dataset, len, y_axis);

  if (x_extent > y_extent) {
    return x_axis;
  }
  return y_axis;
}

float_t
get_dataset_extent(kpoint **arr, int len,
                   int axis) { // suppose datased ordered in axis dimension
  return ((*arr[len - 1]).coords[axis] -
          (*arr[0]).coords[0]); // extent = max - min
}

void get_dataset_ptrs(kpoint *dataset, kpoint **dataset_ptrs, int len) {
  for (int i = 0; i < len; i++) {
    dataset_ptrs[i] = &dataset[i];
  }
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

int choose_splitting_dimension_variance(float_t **x_dataset_ptr,
                                        float_t **y_dataset_ptr, int len) {
  float_t var_x = get_array_variance(x_dataset_ptr, len);
  float_t var_y = get_array_variance(y_dataset_ptr, len);

  if (var_x >= var_y)
    return x_axis;
  return y_axis;
}

float_t get_max_value_dataset(kpoint **arr, int len,
                              int axis) { // not necessary :(
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

float_t get_array_variance(float_t **arr, int len) {
  float_t sum = 0.0, sum1 = 0.0;
  for (int i = 0; i < len; i++) {
    sum += *arr[i];
  }

  float_t average = sum / (float_t)len;
  for (int i = 0; i < len; i++) {
    sum1 = sum1 + pow((*arr[i] - average), 2);
  }

  return sum1 / (float_t)len;
}

void get_axis_coord_ptr(kpoint *dataset, float_t **arr, int len, int axis) {
  for (int i = 0; i < len; i++) {
    arr[i] = &dataset[i].coords[axis];
  }
}

void get_axis_coord(kpoint *dataset, float_t *arr, int len, int axis) {
  for (int i = 0; i < len; i++) {
    arr[i] = dataset[i].coords[axis];
  }
}

int cmpfunc_ptr(const void *a, const void *b) {
  float_t **ptr_a = (float_t **)a;
  float_t **ptr_b = (float_t **)b;

  return ((**ptr_a > **ptr_b) - (**ptr_a < **ptr_b));
}

int cmpfunc(const void *a, const void *b) {
  float_t *ptr_a = (float_t *)a;
  float_t *ptr_b = (float_t *)b;

  return ((*ptr_a > *ptr_b) - (*ptr_a < *ptr_b));
}

float_t nlogn_median(float_t *arr, int len) {
  qsort(arr, len, sizeof(float_t), cmpfunc);
  return arr[len / 2];
}
