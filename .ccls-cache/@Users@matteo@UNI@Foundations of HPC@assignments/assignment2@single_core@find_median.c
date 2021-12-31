#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

float_t nlogn_median(float_t *arr, int len);
void get_axis_coord(kpoint *dataset, float_t *arr, int len, int axis);
void get_axis_coord_ptr(kpoint *dataset, float_t **arr, int len, int axis);
int cmpfunc(const void *a, const void *b);
int cmpfunc_ptr(const void *a, const void *b);
float_t get_array_variance(float_t **arr,int len);

int choose_splitting_dimension(float_t **x_dataset_ptr, float_t **y_dataset_ptr, int len);


int main(int argc, char *argv[]) {

  kpoint dataset[7] = {{2, 3}, {5, 4}, {9, 6}, {6, 22}, {4, 7}, {8, 1}, {7, 2}};
  int len = sizeof(dataset) / sizeof(dataset[0]);
  float_t x_axis_dataset[len];
  float_t y_axis_dataset[len];

  // get_axis_coord(dataset, x_axis_dataset, len, x_axis);
  // get_axis_coord(dataset, y_axis_dataset, len, y_axis);

  // printf("median x_axis: %f\n", nlogn_median(x_axis_dataset, len));
  // printf("median y_axis: %f\n", nlogn_median(y_axis_dataset, len));
  
  //try sorting only pointers
  float_t *x_ordered_dataset_ptr[len];
  float_t *y_ordered_dataset_ptr[len];

  get_axis_coord_ptr(dataset, x_ordered_dataset_ptr, len, x_axis);
  qsort(x_ordered_dataset_ptr, len, sizeof(float_t *), cmpfunc_ptr);
  printf("last x of the ordered dataset by x: %f\n", *x_ordered_dataset_ptr[len-1]);

  get_axis_coord_ptr(dataset, y_ordered_dataset_ptr, len, y_axis);
  qsort(y_ordered_dataset_ptr, len, sizeof(float_t *), cmpfunc_ptr);
  printf("last y of the ordered dataset by y: %f\n", *y_ordered_dataset_ptr[len-1]);

  printf("variance of x components of dataset: %f\n", get_array_variance(x_ordered_dataset_ptr, len));
  printf("variance of y components of dataset: %f\n", get_array_variance(y_ordered_dataset_ptr, len));

  printf("choose splitting dimension (0=x, 1=y) -> %d\n", choose_splitting_dimension(x_ordered_dataset_ptr, y_ordered_dataset_ptr, len));

  return 0;
}

kpoint* choose_splitting_point(kpoint *dataset, float_t **x_dataset_ptr, float_t **y_dataset_ptr, int len, int chosen_axis) {
  if (chosen_axis == 0) {
    return x_dataset_ptr[len/2];
  }
  return y_dataset_ptr[len/2];
}

int choose_splitting_dimension(float_t **x_dataset_ptr, float_t **y_dataset_ptr, int len) {
  float_t var_x = get_array_variance(x_dataset_ptr, len);
  float_t var_y = get_array_variance(y_dataset_ptr, len);

  if (var_x >= var_y) 
    return x_axis;  
  return y_axis;
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
