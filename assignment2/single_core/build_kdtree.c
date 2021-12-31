#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if !defined(DOUBLE_PRECISION)
#define float_t float
#else
#define float_t double
#endif
#define NDIM 2

typedef struct {
  float_t coords[NDIM];
} kpoint;

struct kdnode {
  int axis;
  kpoint split;
  struct kdnode *left, *right;
};

struct kdnode *build_kdtree(kpoint *points, int N, int ndim, int axis);
void quicksort(kpoint *data, int start, int end);
void insertion_sort(kpoint *data, int start, int end);
int compare_ge(const void *A, const void *B);
int compare_ptr(const void *A, const void *B);

int main(int argc, char *argv[]) {

  kpoint root = {{0, 0}};

  kpoint left = {{1, 0}};
  kpoint right = {{0, 2}};

  struct kdnode node1 = {0, left, NULL, NULL};
  struct kdnode node2 = {1, right, NULL, NULL};
  struct kdnode root_node = {0, root, &node1, &node2};

  printf("root node -> %f,%f\n", root_node.split.coords[0],
         root_node.split.coords[1]);
  printf("root node left pointer -> %f,%f\n", root_node.left->split.coords[0],
         root_node.left->split.coords[1]);
  printf("root node right pointer -> %f,%f\n", root_node.right->split.coords[0],
         root_node.right->split.coords[1]);

  kpoint dataset[6] = {{2, 3}, {5, 4}, {9, 6}, {4, 7}, {8, 1}, {7, 2}};

  int dataset_len = (sizeof(dataset) / sizeof(kpoint));

  float_t *ptrs_ordered_by_x[dataset_len];
  float_t *ptrs_ordered_by_y[dataset_len];

  float_t arr[3] = {3.0, 3.5, 1.8};
  float_t *ptr_arr[3] = {&arr[0], &arr[1], &arr[2]};

  printf("%p %p %p \n", ptr_arr[0], ptr_arr[1], ptr_arr[2]);

  qsort(ptr_arr, 3, sizeof(float_t), compare_ptr);

  printf("%p %p %p \n", ptr_arr[0], ptr_arr[1], ptr_arr[2]);
  printf("%f %f %f \n", arr[0], arr[1], arr[2]);
  printf("%d \n", compare_ptr(ptr_arr[2], ptr_arr[0]));

  return 0;
}


void sort_by_axis(int axis, kpoint *dataset) {}

int compare_ge(const void *A, const void *B) {
  float_t *a = (float_t *)A;
  float_t *b = (float_t *)B;

  return (a >= b);
}

int compare_ptr(const void *a, const void *b) {
  float_t *ptr_a = (float_t *)a;
  float_t *ptr_b = (float_t *)b;

  return ((*ptr_a > *ptr_b) - (*ptr_a < *ptr_b));
}

float_t nlogn_median(float_t *arr, int len) {
  qsort(arr, len, sizeof(float_t), compare_ge);
  return arr[len / 2];
}

