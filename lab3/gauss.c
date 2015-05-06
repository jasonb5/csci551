#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_RAND 1.0e2 //1.0e6

#define DEBUG

#ifdef DEBUG
#define dprintf(X, ...) fprintf(stdout, "[DEBUG] " X, ##__VA_ARGS__)
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
  double sum;
  double factor;
  double interm;
  double l2norm;
  double *result;
  double **matrix;
  double **matrix_a;

  srand(time(NULL));
  srand48(time(NULL));

  n = atoi(argv[1]);

  result = malloc(n * sizeof(double));

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

  dprintf("Forward elimination\n");
    
  for (i = 0; i < n - 1; ++i) {
    print_aug_matrix(matrix, n);

    swap_pivot(matrix, i, n);

    print_aug_matrix(matrix, n);

    for (j = i + 1; j < n; ++j) {
      factor = matrix[j][i] / matrix[i][i];

      dprintf("Factor %f = %f / %f\n", factor, matrix[j][i], matrix[i][i]);

      for (k = 0; k < n + 1; ++k) {
        interm = factor * matrix[i][k];

        dprintf("Interm %f = %f * %f\n", interm, factor, matrix[i][k]);
        dprintf("\t%f = %f - %f\n", matrix[j][k]-interm, matrix[j][k], interm);

        matrix[j][k] -= interm;
      }
    }

    print_aug_matrix(matrix, n);
  }

  dprintf("Backwards substitution\n"); 

  for (i = n - 1; i >= 0; --i) {
    result[i] = matrix[i][n];

    for (j = n - 1; j > i; --j) {
      dprintf("%f = %f * %f\n", matrix[i][j] * result[j], matrix[i][j], result[j]);

      result[i] -= matrix[i][j] * result[j]; 
    }

    result[i] /= matrix[i][i];
      
    dprintf("%f\n", result[i]);
  }

  for (i = 0; i < n; ++i) {
    free(matrix[i]);
  }

  free(matrix);

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
