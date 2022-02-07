#include <mpi.h>
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

kpoint *generate_dataset(int len);
void copy_dataset(kpoint *dataset, kpoint *new_dataset, int len);

// kd-tree build functions
struct kdnode *build_kdtree(kpoint *dataset, float_t extremes[NDIM][2], int len, int axis, int level);
struct kdnode *build_kdtree_until_level_then_scatter(kpoint *dataset, float_t extremes[NDIM][2], int len, int axis, int level, int final_level, int counter);
int choose_splitting_dimension(float_t extremes[NDIM][2]);
int choose_splitting_point(kpoint *dataset, float_t extremes[NDIM][2], int len, int chosen_axis);
void get_dataset_extremes(kpoint *dataset, float_t extremes[NDIM][2], int len, int axis);
void copy_extremes(float_t old_extremes[NDIM][2], float_t new_extremes[NDIM][2]);

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

  printf("print provided %d\n", *provided);

  tstart = CPU_TIME;

  #pragma omp parallel shared(nthreads)
  {
    #pragma omp single
    {
      nthreads = omp_get_num_threads();
      printf("I am mpi process %d an I have %d threads\n", my_rank, nthreads);
    }
  }

  if (my_rank == 0) {
    // kpoint dataset[16] = {{2, 3}, {5, 4}, {9, 6}, {6, 22}, {4, 7},
    //                      {8, 1}, {7, 2}, {8, 9}, {1, 1}, {0.55, 6},
    //                      {1.11111, 8.432}, {9, 2.7}, {4.3253, 42}, {0.1224, 0.4635}, {0.124, 0.773}, {0.4336, 0.7456}};
    // len = sizeof(dataset) / sizeof(dataset[0]);

    if (argc == 2 ) {
      len = atoi(argv[1]);
    } else if (argc == 1){
      len = 100000000;
    } else {
      return 1;
    }

    kpoint *dataset = generate_dataset(len);

    printf("len: %d\n", len);

    struct kdnode *root;

    tstart = CPU_TIME;

    printf("choosing splitting dimension..\n");
    
    float_t extremes[NDIM][2]; //min_value (index 0) and max value (index 1) in each dimension NDIM

    get_dataset_extremes(dataset, extremes, len, x_axis);
    get_dataset_extremes(dataset, extremes, len, y_axis);
   
    int chosen_axis = choose_splitting_dimension(extremes);
   
    int final_level = log2(numprocs);
    
    printf("start building kdtree..\n");
    #pragma omp parallel shared(dataset, root) firstprivate(extremes, chosen_axis, len, final_level) 
    {
      #pragma omp master
      {
        int current_level = 0, counter = 0;
        root = build_kdtree_until_level_then_scatter(dataset, extremes, len, chosen_axis, current_level, final_level, counter);
        printf("finished build kd_tree on process %d \n", my_rank);
      }
    }

    printf("mpi process %d has root node is %f,%f\n", my_rank, root->split.coords[0], root->split.coords[1]);

    free(dataset);
  } else {
    int recv_len;
    MPI_Status status;

    MPI_Recv(&recv_len, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    printf("mpi process %d received len %d\n", my_rank,  recv_len);

    kpoint *recv_dataset = malloc(recv_len * sizeof(kpoint));
    MPI_Recv(recv_dataset, recv_len * sizeof(kpoint), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);

    float_t recv_extremes[NDIM][2] = {};
    MPI_Recv(recv_extremes, NDIM * 2 * sizeof(float_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);

    int recv_axis;
    MPI_Recv(&recv_axis, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

    int level = 0;
    struct kdnode *chunk_root;
    
    printf("i am mpi process %d, start building my kd-tree..\n", my_rank);

    #pragma omp parallel shared(recv_dataset, chunk_root) firstprivate(recv_extremes, recv_axis, recv_len) 
    {
      #pragma omp master
      {
        int current_level = 0;
        chunk_root = build_kdtree(recv_dataset, recv_extremes, recv_len, recv_axis, current_level);
      }
    }

    printf("i am mpi process %d, my chunk root node is %f,%f\n", my_rank, chunk_root->split.coords[0], chunk_root->split.coords[1]);
    free(recv_dataset);
  } 

  double telapsed = CPU_TIME - tstart;
  printf("elapsed time of mpi process %d: %f\n", my_rank, telapsed);

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

struct kdnode *build_kdtree(kpoint *dataset, float_t extremes[NDIM][2], int len, int previous_axis, int level) {
  if (len == 1) {
    struct kdnode *leaf = malloc(sizeof(struct kdnode));
    leaf->axis = previous_axis;
    leaf->split = dataset[0];
    leaf->left = NULL;
    leaf->right = NULL;

    return leaf;
  }

  struct kdnode *node = malloc(sizeof(struct kdnode));

  int chosen_axis = choose_splitting_dimension(extremes);

  int median_idx = choose_splitting_point(dataset, extremes, len, chosen_axis);
  node->axis = chosen_axis;
  node->split = dataset[median_idx];

  kpoint *left_points, *right_points;

  int len_left = median_idx - 1;    // length of the left points
  int len_right = len - median_idx; // length of the right points
  
  if (len_left != 0) {
    left_points = &dataset[0];       // starting pointer of left_points

    extremes[chosen_axis][0] = dataset[0].coords[chosen_axis]; //min value of chosen axis for left points
    extremes[chosen_axis][1] = dataset[len_left - 1].coords[chosen_axis]; //max value of chosen axis for left points

    #pragma omp task shared(left_points) firstprivate(extremes, len_left, chosen_axis, level) if(len_left >= build_cutoff) mergeable untied
      node->left = build_kdtree(left_points, extremes, len_left, chosen_axis, level+1);
  }

  right_points = &dataset[median_idx]; // starting pointer of right_points
  
  extremes[chosen_axis][0] = dataset[median_idx].coords[chosen_axis]; //min value of chosen axis for right points
  extremes[chosen_axis][1] = dataset[len - 1].coords[chosen_axis]; //max value of chosen axis for right points

  #pragma omp task shared(right_points) firstprivate(extremes, len_right, chosen_axis, level) if(len_right >= build_cutoff) mergeable untied
    node->right = build_kdtree(right_points, extremes, len_right, chosen_axis, level+1);
  
  #pragma omp taskwait
  return node;
}

struct kdnode *build_kdtree_until_level_then_scatter(kpoint *dataset, float_t extremes[NDIM][2], int len, int previous_axis, int current_level, int final_level, int counter) {
  if ((current_level == final_level) && (counter != 0)) {

    //TODO: works only for numprocs = 2 for now 
    MPI_Send(&len, 1, MPI_INT, counter, 0, MPI_COMM_WORLD);
    printf("sent chunk_length from mpi process 0 to mpi process %d, len %d\n", counter, len);

    kpoint *chunk = malloc(len * sizeof(kpoint));
    copy_dataset(chunk, dataset, len);

    MPI_Send(chunk, len * sizeof(kpoint), MPI_BYTE, counter, 0, MPI_COMM_WORLD);
    MPI_Send(extremes, NDIM * 2 * sizeof(float_t), MPI_BYTE, counter, 0, MPI_COMM_WORLD);
    MPI_Send(&previous_axis, 1, MPI_INT, counter, 0, MPI_COMM_WORLD);

    printf("sent chunk from mpi process 0, thread %d, to mpi process %d\n", omp_get_thread_num(), counter);

    free(chunk);
    return NULL;
  }

  if (len == 1) {
    struct kdnode *leaf = malloc(sizeof(struct kdnode));
    leaf->axis = previous_axis;
    leaf->split = dataset[0];
    leaf->left = NULL;
    leaf->right = NULL;

    return leaf;
  }

  struct kdnode *node = malloc(sizeof(struct kdnode));

  int chosen_axis = choose_splitting_dimension(extremes);

  int median_idx = choose_splitting_point(dataset, extremes, len, chosen_axis);
  
  node->axis = chosen_axis;
  node->split = dataset[median_idx];

  kpoint *left_points, *right_points;

  int len_left = median_idx;    // length of the left points
  int len_right = len - (median_idx + 1); // length of the right points
  
  if (len_left != 0) {
    left_points = &dataset[0];       // starting pointer of left_points

    extremes[chosen_axis][0] = dataset[0].coords[chosen_axis]; //min value of chosen axis for left points
    extremes[chosen_axis][1] = dataset[len_left - 1].coords[chosen_axis]; //max value of chosen axis for left points

    #pragma omp task shared(left_points, counter) firstprivate(extremes, len_left, chosen_axis, current_level, final_level) if(len_left >= build_cutoff) mergeable untied
      node->left = build_kdtree_until_level_then_scatter(left_points, extremes, len_left, chosen_axis, current_level+1, final_level, counter+1);
  }
  
  if(len_right != 0) {
    right_points = &dataset[median_idx]; // starting pointer of right_points

    extremes[chosen_axis][0] = dataset[median_idx].coords[chosen_axis]; //min value of chosen axis for right points
    extremes[chosen_axis][1] = dataset[len - 1].coords[chosen_axis]; //max value of chosen axis for right points

    #pragma omp task shared(right_points, counter) firstprivate(extremes, len_right, chosen_axis, current_level, final_level) if(len_right >= build_cutoff) mergeable untied
      node->right = build_kdtree_until_level_then_scatter(right_points, extremes, len_right, chosen_axis, current_level+1, final_level, counter+0);
  }
  
  #pragma omp taskwait
  return node;
}

int choose_splitting_point(kpoint *dataset, float_t extremes[NDIM][2], int len, int chosen_axis) {

  float_t *distances = malloc(len * sizeof(float_t));
  float_t computed_median = (extremes[chosen_axis][1] + extremes[chosen_axis][0]) / 2.0;

  #pragma omp parallel for shared(dataset, distances) firstprivate(computed_median, chosen_axis) schedule(static) proc_bind(close)
    for (int i = 0; i < len; i++) {
      distances[i] = fabs(dataset[i].coords[chosen_axis] - computed_median);
    }

    struct Compare { float_t val; int index; };    
    struct Compare min_distance = {10, 0};
    #pragma omp declare reduction(minimum : struct Compare : omp_out = omp_in.val < omp_out.val ? omp_in : omp_out)\
                          initializer(omp_priv = {100, 0})
   
    #pragma omp parallel for shared(distances) reduction(minimum:min_distance) schedule(static) proc_bind(close)
      for (int i = 0; i < len; i++) {
        if (distances[i] < min_distance.val) {
          min_distance.val = distances[i];
          min_distance.index = i;
        }
      }
  free(distances);
  printf("min distance index: %d\n", min_distance.index);
  return min_distance.index;
}

int choose_splitting_dimension(float_t extremes[NDIM][2]) {
  float_t x_extent = extremes[x_axis][1] - extremes[x_axis][0]; //max value - min value
  float_t y_extent = extremes[y_axis][1] - extremes[y_axis][0]; //max value - min value

  if (x_extent > y_extent) {
    return x_axis;
  }
  return y_axis;
}

void get_dataset_extremes(kpoint *dataset, float_t extremes[NDIM][2], int len, int axis) {
  float_t max_value = dataset[0].coords[axis];
  float_t min_value = max_value;
  #pragma omp parallel for reduction(max:max_value) reduction(min:min_value) schedule(static) proc_bind(close)
  for (int i = 1; i < len; i++) {
    max_value = max_value > dataset[i].coords[axis] ? max_value : dataset[i].coords[axis];
    min_value = min_value < dataset[i].coords[axis] ? min_value : dataset[i].coords[axis];
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

void copy_dataset(kpoint *new_dataset, kpoint *dataset, int len) {
  for (int i = 0; i < len; i++) {
    new_dataset[i].coords[0] = dataset[i].coords[0];
    new_dataset[i].coords[1] = dataset[i].coords[1];
  }
}
