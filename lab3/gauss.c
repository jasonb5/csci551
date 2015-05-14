/**
 * @file gauss.c
 * @author Jason Boutte <jboutte@mail.csuchico.edu>
 * 
 * @brief Parallelization of gaussian elimination.
 */
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

/**
 * @brief Finds pivot and performs a swap.
 *
 * @param matrix n x n matrix
 * @param col Starting column for search
 * @param n Size of matrix
 */
void swap_pivot(double **matrix, int col, int n);

/**
 * @brief Generates augmented matrix n x n 
 *
 * @param matrix n x n matrix
 * @param n Size of matrix
 */ 
void generate_aug_matrix(double **matrix, int n);

/**
 * @brief Prints augmented matrix
 *
 * @param matrix n x n matrix
 * @param n Size of matrix
 */
void print_aug_matrix(double **matrix, int n);

/**
 * @brief Prints program usage
 *
 * @param prog Program name in a C-style string
 */
void print_usage(char *prog);

int main(int argc, char **argv) {
  // Check if we've got the correct number of arguments
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

  // Seed rng for signing numbers
  srand(time(NULL));

  // Seed rng for double precision floating point
  srand48(time(NULL));

  // Get size of the matrix
  n = atoi(argv[1]);

  // Get the number of threads to parallelize with
  threads = atoi(argv[2]);

  // Initialize vectors
  result = malloc(n * sizeof(double));

  result_res = malloc(n * sizeof(double));

  memset(result, 0, n * sizeof(double));

  // Initialize matrices
  matrix = malloc(n * sizeof(double));

  matrix_a = malloc(n * sizeof(double));

  for (i = 0; i < n; ++i) {
    matrix[i] = malloc((n + 1) * sizeof(double));

    matrix_a[i] = malloc((n + 1) * sizeof(double));
  }

  // Generate augmented matrix
  generate_aug_matrix(matrix, n);

  // Create copy of original matrix
  for (i = 0; i < n; ++i) {
    memcpy(matrix_a[i], matrix[i], (n + 1) * sizeof(double));
  }

  // Start timining
  start = omp_get_wtime();
   
  // Perform forward elimination 
  for (i = 0; i < n - 1; ++i) {
    swap_pivot(matrix, i, n);

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

  // Perform back substitution
  for (i = n - 1; i >= 0; --i) {
    interm = matrix[i][n];

#pragma omp parallel for private(j) shared(n, i, matrix, result) \
  default(none) num_threads(threads) reduction(-:interm)
    for (j = n - 1; j > i; --j) {
      interm -= matrix[i][j] * result[j]; 
    }

    result[i] = interm /  matrix[i][i];
  }

  // Stop timing
  end = omp_get_wtime();

// Enter parallel region
#pragma omp parallel default(none) num_threads(threads) \
  private(j, interm) shared(i, n, matrix_a, result, result_res) 
{
  // Get current thrad number
  int id = omp_get_thread_num();

  // Print if we're the master thread
  if (id == 0) {
    printf("number of cores - %i\n", omp_get_num_procs());
    
    // omp_get_num_threads must be called in a parallel region
    // otherwise it will always return one
    printf("number of threads - %i\n", omp_get_num_threads());
  }

// Compute the residual vector as Ax-b
// A = matrix_a
// x = result stored as n x [0..n-1] 
// b = result stored as n x [n]
#pragma omp for   
  for (i = 0; i < n; ++i) {
    interm = 0;

    for (j = 0; j < n; ++j) {
      interm += matrix_a[i][j] * result[j];
    }

    result_res[i] = interm - matrix_a[i][n];
  }
}
  // Calculate the l2-norm
  // as sqrt(x0^2+x1^2+...+xn^2)
  l2norm = 0;

  for (i = 0; i < n; ++i) {
    l2norm += pow(fabs(result_res[i]), 2); 
  }

  l2norm = sqrt(l2norm);

  // Print the results
  printf("l2-norm %.16f\n", l2norm);
  printf("elapsed time %f\n", end-start);

  // Clean up the memory allocations
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

/**
 * @brief Finds the pivot and swaps rows.
 *
 * Using col as the index we start at the diaganol defined by the index
 * and continue down the column looking for the largest absolute value.
 * Once this value has been found we swap the row containing it and the 
 * row of the diaganol index.
 */ 
void swap_pivot(double **matrix, int col, int n) {
  int i;
  int pivot;
  double max;
  double value;
  double *temp;

  pivot = 0;

  // Iterate down the column look for the largest absolute value
  for (i = col; i < n; ++i) {
    value = fabs(matrix[i][col]);

    if (value > max) {
      pivot = i;

      max = value;
    } 
  }

  // No need to swap of the pivot is the row we started at
  if (pivot == col) return;

  dprintf("Swapping row %i with %i\n", pivot, col);

  // Swap pivot and index row
  temp = matrix[pivot];

  matrix[pivot] = matrix[col];

  matrix[col] = temp;
}

/**
 * @brief Generates augmented matrix.
 *
 * Generate an augmented matrix of size n x n+1.
 * To determine the sign of a value we use a  call to rand().
 */
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
/**
 * @brief Prints augmented matrix
 *
 */
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

/**
 * @brief Prints program usage.
 *
 */
void print_usage(char *prog) {
  printf("%s n threads\n", prog); 
  printf("n - Size of the matrix\n");
  printf("threads - Number of threads to run on\n");
}
