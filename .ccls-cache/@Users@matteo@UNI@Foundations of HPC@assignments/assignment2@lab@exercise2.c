#if defined(__STDC__)
#if (__STDC_VERSION__ >= 199901L)
#define _XOPEN_SOURCE 700
#endif
#endif
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define THRESHOLD 2000

int more_data_arriving(int);
int getting_data(int **);
double heavy_work(int);

int main(int argc, char **argv) {

  srand48(time(NULL));

  // you are free to decide which
  // variable are shared and (thread)private
  //
  int nthreads;
  int *data;
  int ndata;
  int bunch = 10;
  int next_bunch = 0;
  int iterations = 0;
  int data_are_arriving = more_data_arriving(iterations);

  nthreads = omp_get_num_threads();
  data = calloc(nthreads, sizeof(int));

#pragma omp parallel
  {
    int me = omp_get_thread_num();

    while (data_are_arriving) {
#pragma omp single
      {
        ndata = getting_data(&data);
        printf("iteration %d: thread %d got %d data\n", iterations, me, ndata);
      }
      // code of version 2 here
      int mystart;
      do {
        int mystop;
#pragma omp atomic capture
        {
          mystart = next_bunch;
          next_bunch += bunch;
        }
        mystop = mystart + bunch;
        mystop  = (mystop > ndata ? ndata : mystop);
      } while (mystart < ndata);

#pragma omp single // remember that there is an implicit barrier at the end of
                   // single region
      {
        data_are_arriving = more_data_arriving(iterations + 1);
        if (data_are_arriving) {
          iterations++;
        } else {
          printf("\t>>> iteration %d : thread %d got the news that "
                 "no more data will arrive\n",
                 iterations, me);
        }
      }
    }
  }

  return 0;
}

int more_data_arriving(int i) {
  // it is increasingly probable that
  // no more data arrive when i approaches
  // THRESHOLD
  //
  double p = (double)(THRESHOLD - i) / THRESHOLD;
  return (drand48() < p);
}

int getting_data(int **data) {
#define MIN 10
#define MAX 25
#define MAX_DATA 40

  // produces no more than n-1
  // data
  int howmany = lrand48() % MAX_DATA;
  howmany = (howmany == 0 ? 1 : howmany);

  // be sure that the data
  // array has enough room
  // to host up to n-1 data
  *data = (int *)calloc(howmany, sizeof(int));

  for (int j = 0; j < howmany; j++)
    (*data)[j] = 1024 + lrand48() % (MAX - MIN); // values will range
                                                 // from MIN up to MAX

  return howmany;
}

double heavy_work(int N) {
  double guess = 3.141572 / 3 * N;

  for (int i = 0; i < N; i++) {
    guess = exp(guess);
    guess = sin(guess);
  }
  return guess;
}
