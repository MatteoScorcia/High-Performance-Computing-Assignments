
/* ────────────────────────────────────────────────────────────────────────── *
 │                                                                            │
 │ This file is part of the exercises for the Lectures on                     │
 │   "Foundations of High Performance Computing"                              │
 │ given at                                                                   │
 │   Master in HPC and                                                        │
 │   Master in Data Science and Scientific Computing                          │
 │ @ SISSA, ICTP and University of Trieste                                    │
 │                                                                            │
 │ contact: luca.tornatore@inaf.it                                            │
 │                                                                            │
 │     This is free software; you can redistribute it and/or modify           │
 │     it under the terms of the GNU General Public License as published by   │
 │     the Free Software Foundation; either version 3 of the License, or      │
 │     (at your option) any later version.                                    │
 │     This code is distributed in the hope that it will be useful,           │
 │     but WITHOUT ANY WARRANTY; without even the implied warranty of         │
 │     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          │
 │     GNU General Public License for more details.                           │
 │                                                                            │
 │     You should have received a copy of the GNU General Public License      │
 │     along with this program.  If not, see <http://www.gnu.org/licenses/>   │
 │                                                                            │
 * ────────────────────────────────────────────────────────────────────────── */

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
int getting_data(int, int *);
double heavy_work(int);

int main(int argc, char **argv) {

  srand48(time(NULL));

  // you are free to decide which
  // variable are shared and (thread)private
  //
  int nthreads;
  int *data;
  int ndata;
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
        ndata = getting_data(nthreads, data);
        printf("iteration %d: thread %d got %d data\n", iterations, me, ndata);
      }

      if (me < ndata) {
        heavy_work(data[me]);
      }

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

int getting_data(int n, int *data) {
#define MIN 1000
#define MAX 10000

  // produces no more than n-1
  // data
  int howmany = lrand48() % n;
  howmany = (howmany == 0 ? 1 : howmany);

  // be sure that the data
  // array has enough room
  // to host up to n-1 data
  for (int j = 0; j < howmany; j++)
    data[j] = 1024 + lrand48() % (MAX - MIN); // values will range
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
