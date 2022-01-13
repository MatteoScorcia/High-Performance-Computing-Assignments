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
int cmpfunc_x_axis(const void *a, const void *b);
int cmpfunc_y_axis(const void *a, const void *b);

// actual functions used to work with dataset pointers
void get_dataset_ptrs(kpoint *dataset, kpoint **dataset_ptrs, int len);
float_t get_dataset_extent(kpoint **arr, int len, int axis);
int choose_splitting_dimension(kpoint **x_ordered_dataset,
                               kpoint **y_ordered_dataset, int len);
kpoint *choose_splitting_point(kpoint **x_ordered_dataset,
                               kpoint **y_ordered_dataset, int len,
                               int chosen_axis);
struct kdnode *build_kdtree(kpoint *dataset, int len, int ndim, int axis);
int cmpfunc_ptr(const void *a, const void *b);

int main(int argc, char *argv[]) {

  kpoint dataset[7] = {{2, 3}, {5, 4}, {9, 6}, {6, 22}, {4, 7}, {8, 1}, {7, 2}};
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

  // --------------------------------   working separately on the 2 dimensions
  // printf("\nworking separately on the 2 dimensions ----- \n");
  // float_t *x_coords_dataset[len];
  // float_t *y_coords_dataset[len];
  //
  // get_axis_coord_ptr(dataset, x_coords_dataset, len, x_axis);
  // qsort(x_coords_dataset, len, sizeof(float_t *), cmpfunc_ptr);
  // printf("last x of the ordered dataset by x: %f\n",
  //        *x_coords_dataset[len - 1]);
  //
  // get_axis_coord_ptr(dataset, y_coords_dataset, len, y_axis);
  // qsort(y_coords_dataset, len, sizeof(float_t *), cmpfunc_ptr);
  // printf("last y of the ordered dataset by y: %f\n",
  //        *y_coords_dataset[len - 1]);
  //
  // printf("variance of x components of dataset: %f\n",
  //        get_array_variance(x_coords_dataset, len));
  // printf("variance of y components of dataset: %f\n",
  //        get_array_variance(y_coords_dataset, len));
  //
  // printf("choose splitting dimension (0=x, 1=y) -> %d\n",
  //        choose_splitting_dimension_variance(x_coords_dataset,
  //        y_coords_dataset,
  //                                            len));

  return 0;
}

struct kdnode *build_kdtree(kpoint *dataset, int len, int ndim, int axis) {
  if (len == 1) {
    struct kdnode leaf = {axis, dataset[0], NULL, NULL};
    return &leaf;
  }
}

kpoint *choose_splitting_point(kpoint **x_ordered_dataset,
                               kpoint **y_ordered_dataset, int len,
                               int chosen_axis) {
  if (chosen_axis == 0) {
    return x_ordered_dataset[len / 2];
  }
  return y_ordered_dataset[len / 2];
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

int cmpfunc_ptr(const void *a, const void *b) {
  float_t **ptr_a = (float_t **)a;
  float_t **ptr_b = (float_t **)b;

  return ((**ptr_a > **ptr_b) - (**ptr_a < **ptr_b));
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

int cmpfunc(const void *a, const void *b) {
  float_t *ptr_a = (float_t *)a;
  float_t *ptr_b = (float_t *)b;

  return ((*ptr_a > *ptr_b) - (*ptr_a < *ptr_b));
}

float_t nlogn_median(float_t *arr, int len) {
  qsort(arr, len, sizeof(float_t), cmpfunc);
  return arr[len / 2];
}
