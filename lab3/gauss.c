#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_RAND 1.0e2 //1.0e6

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
  double factor;
  double interm;
  double *result;
  double **matrix;

  srand(time(NULL));
  srand48(time(NULL));

  n = atoi(argv[1]);

  result = malloc(n * sizeof(double));

  memset(result, 0, n * sizeof(double));

  matrix = malloc(n * sizeof(double));

  for (i = 0; i < n; ++i) {
    matrix[i] = malloc((n + 1) * sizeof(double));
  }

  generate_aug_matrix(matrix, n);

  for (i = 0; i < n; ++i) {
    swap_pivot(matrix, i, n);

    print_aug_matrix(matrix, n);

    for (j = i + 1; j < n; ++j) {
      factor = matrix[j][i] / matrix[i][i];

      printf("Factor %f\n", factor);

      for (k = 0; k < n + 1; ++k) {
        interm = factor * matrix[i][k];

        matrix[j][k] -= interm;
      }
    }
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

  for (i = col; i < n; ++i) {
    value = fabs(matrix[i][col]);

    if (value > max) {
      pivot = i;

      max = value;
    } 
  }

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

void print_aug_matrix(double **matrix, int n) {
  int i, j;

  printf("\n");

  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      printf("%.2f\t", matrix[i][j]);
    }
   
    printf("\tx%i\t%.2f\n", i+1, matrix[i][n]);
  }

  printf("\n");
}

void print_usage(char *prog) {
  printf("%s n threads\n", prog); 
  printf("n - Size of the matrix\n");
  printf("threads - Number of threads to run on\n");
}
