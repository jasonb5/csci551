#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_RAND 1.0e6

#ifdef DEBUG
#ifdef _OPENMP
#define dprintf(X, ...) fprintf(stdout, "[DEBUG][%i] " X, omp_get_thread_num(), ##__VA_ARGS__)
#else
#define dprintf(X, ...) fprintf(stdout, "[DEBUG] " X, ##__VA_ARGS__)
#endif
#else
#define dprintf(X, ...)
#endif

void swap_pivot(double **matrix, int col, int n);
void generate_aug_matrix(double **matrix, int n);
void print_aug_matrix(double **matrix, int n);
void print_usage(char *prog);

int main(int argc, char **argv) {
  if (argc != 3) {
    print_usage(argv[0]); 

    exit(1);
  }

  int n;
  int i, j, k;
  int threads;
  double factor;
  double interm;
  double l2norm;
  double start, end;
  double *result;
  double *result_res;
  double **matrix;
  double **matrix_a;

  srand(time(NULL));
  srand48(time(NULL));

  n = atoi(argv[1]);

  threads = atoi(argv[2]);

  result = malloc(n * sizeof(double));

  result_res = malloc(n * sizeof(double));

  memset(result, 0, n * sizeof(double));

  matrix = malloc(n * sizeof(double));

  matrix_a = malloc(n * sizeof(double));

  for (i = 0; i < n; ++i) {
    matrix[i] = malloc((n + 1) * sizeof(double));

    matrix_a[i] = malloc((n + 1) * sizeof(double));
  }

  generate_aug_matrix(matrix, n);

  for (i = 0; i < n; ++i) {
    memcpy(matrix_a[i], matrix[i], (n + 1) * sizeof(double));
  }
    
  for (i = 0; i < n - 1; ++i) {
    swap_pivot(matrix, i, n);

  start = omp_get_wtime();

#pragma omp parallel for private(j, k, factor, interm) \
  shared(i, n, matrix) default(none) num_threads(threads)
    for (j = i + 1; j < n; ++j) {
      factor = matrix[j][i] / matrix[i][i];

      for (k = 0; k < n + 1; ++k) {
        interm = factor * matrix[i][k];

        matrix[j][k] -= interm;
      }
    }
  }

  for (i = n - 1; i >= 0; --i) {
    interm = matrix[i][n];

#pragma omp parallel for private(j) shared(n, i, matrix, result) \
  default(none) num_threads(threads) reduction(-:interm)
    for (j = n - 1; j > i; --j) {
      interm -= matrix[i][j] * result[j]; 
    }

    result[i] = interm /  matrix[i][i];
  }

  end = omp_get_wtime();

#pragma omp parallel default(none) num_threads(threads) \
  private(j, interm) shared(i, n, matrix_a, result, result_res) 
{
  int id = omp_get_thread_num();

  if (id == 0) {
    printf("number of cores - %i\n", omp_get_num_procs());
    printf("number of threads - %i\n", omp_get_num_threads());
  }

#pragma omp for   
  for (i = 0; i < n; ++i) {
    interm = 0;

    for (j = 0; j < n; ++j) {
      interm += matrix_a[i][j] * result[j];
    }

    result_res[i] = interm - matrix_a[i][n];
  }
}
  l2norm = 0;

  for (i = 0; i < n; ++i) {
    l2norm += pow(fabs(result_res[i]), 2); 
  }

  l2norm = sqrt(l2norm);

  printf("l2-norm %.16f\n", l2norm);
  printf("elapsed time %f\n", end-start);

  free(result);

  free(result_res);

  for (i = 0; i < n; ++i) {
    free(matrix[i]);

    free(matrix_a[i]);
  }

  free(matrix);

  free(matrix_a);

  return 0;
}

void swap_pivot(double **matrix, int col, int n) {
  int i;
  int pivot;
  double max;
  double value;
  double *temp;

  pivot = 0;

  for (i = col; i < n; ++i) {
    value = fabs(matrix[i][col]);

    if (value > max) {
      pivot = i;

      max = value;
    } 
  }

  if (pivot == col) return;

  dprintf("Swapping row %i with %i\n", pivot, col);

  temp = matrix[pivot];

  matrix[pivot] = matrix[col];

  matrix[col] = temp;
}

void generate_aug_matrix(double **matrix, int n) {
  int i, j;

  for (i = 0; i < n; ++i) {
    for (j = 0; j < n + 1; ++j) {
      if (rand() % 2 == 0) {
        matrix[i][j] = drand48() * MAX_RAND;
      } else {
        matrix[i][j] = -1.0 * drand48() * MAX_RAND;
      }
    }
  }
}

#ifdef DEBUG
void print_aug_matrix(double **matrix, int n) {
  int i, j, t;
  static char buf[1024];

  dprintf("\n");

  for (i = 0; i < n; ++i) {
    for (j = 0, t = 0; j < n; ++j) {
      t += sprintf(&buf[t], "%.2f\t", matrix[i][j]);
    }
   
    dprintf("%s\tx%i\t%.2f\n", buf, i+1, matrix[i][n]);
  }

  dprintf("\n");
}
#else
void print_aug_matrix(double **matrix, int n) { return; }
#endif

void print_usage(char *prog) {
  printf("%s n threads\n", prog); 
  printf("n - Size of the matrix\n");
  printf("threads - Number of threads to run on\n");
}
