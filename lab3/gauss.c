#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_RAND 1.0e2 //1.0e6

void swap_pivot(double **matrix, int col, int n);
void generate_matrix(double **matrix, int n);
void generate_vector(double *vector, int n);
void print_matrix(double **matrix, int n);
void print_vector(double *vector, int n);
void print_usage(char *prog);

int main(int argc, char **argv) {
  if (argc != 3) {
    print_usage(argv[0]); 

    exit(1);
  }

  int n, i, j, k;
  double *vector;
  double **matrix;
  double factor, interm;

  srand(time(NULL));
  srand48(time(NULL));

  n = atoi(argv[1]);

  matrix = malloc(n * sizeof(double));

  for(i = 0; i < n; ++i) {
    matrix[i] = malloc(n * sizeof(double));
  } 

  vector = malloc(n * sizeof(double));

  generate_matrix(matrix, n);
  generate_vector(vector, n);

  for (i = 0; i < n; ++i) { 
    print_matrix(matrix, n);
    printf("\n");
    swap_pivot(matrix, i, n);
    printf("\n");
    print_matrix(matrix, n);

    for (j = i + 1; j < n; ++j) {
      factor = matrix[j][i] / matrix[i][i];

      printf("\n%f\n", factor); 

      for (k = 0; k < n; ++k) {
        interm = matrix[i][k] * factor;

        matrix[j][k] -= interm;
      }

      printf("\n");
      print_matrix(matrix, n);
    }
  } 

  for (i = 0; i < n; ++i) {
    free(matrix[i]);
  }

  free(matrix);
  free(vector);

  return 0;
}

void swap_pivot(double **matrix, int col, int n) {
  int i;
  int pivot; 
  double max;
  double value;
  double *temp;

  max = 0;

  for (i = col; i < n; ++i) {
    value = fabs(matrix[i][col]);

    if (value > max) {
      pivot = i;

      max = value;
    } 
  }

  printf("Swapping row %i with row %i\n", pivot+1, col+1);

  temp = matrix[col];

  matrix[col] = matrix[pivot];

  matrix[pivot] = temp;
}

void generate_matrix(double **matrix, int n) {
  int i, j;
  
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      if (rand() % 2 == 0) {
        matrix[i][j] = -1.0 * drand48() * MAX_RAND; 
      } else {
        matrix[i][j] = drand48() * MAX_RAND; 
      }
    }
  }
}

void generate_vector(double *vector, int n) {
  int i;

  for (i = 0; i < n; ++i) {
    if (rand() %  2 == 0) {
      vector[i] = -1.0 * drand48() * MAX_RAND;
    } else {
      vector[i] = drand48() * MAX_RAND;
    }
  }
}

void print_matrix(double **matrix, int n) {
  int i, j;

  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      printf("%.2f\t", matrix[i][j]);
    }

    printf("\n");
  }
}

void print_vector(double *vector, int n) {
  int i;

  for (i = 0; i < n; ++i) {
    printf("%.2f\n", vector[i]);
  }
}

void print_usage(char *prog) {
  printf("%s n threads\n", prog); 
  printf("n - Size of the matrix\n");
  printf("threads - Number of threads to run on\n");
}
