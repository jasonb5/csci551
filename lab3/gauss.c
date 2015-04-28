#include <stdio.h>
#include <stdlib.h>

#define MAX_RAND 1.0e6

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

  int n, i;
  double *vector;
  double **matrix;

  srand48(time(NULL));

  n = atoi(argv[1]);

  matrix = malloc(n * sizeof(double));

  for(i = 0; i < n; ++i) {
    matrix[i] = malloc(n * sizeof(double));
  } 

  vector = malloc(n * sizeof(double));

  generate_matrix(matrix, n);
  generate_vector(vector, n);

  print_matrix(matrix, n);
  printf("\n");
  print_vector(vector, n);

  for (i = 0; i < n; ++i) {
    free(matrix[i]);
  }

  free(matrix);
  free(vector);

  return 0;
}

void generate_matrix(double **matrix, int n) {
  int i, j;
  
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      matrix[i][j] = drand48() * MAX_RAND; 
    }
  }
}

void generate_vector(double *vector, int n) {
  int i;

  for (i = 0; i < n; ++i) {
    vector[i] = drand48() * MAX_RAND;
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
