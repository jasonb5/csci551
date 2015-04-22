#include <stdio.h>
#include <stdlib.h>

#define MAX_RAND 1.0e6

void generate_matrix(double *matrix, int n, int m);
void print_matrix(double *matrix, int n, int m);
void print_usage(char *prog);

int main(int argc, char **argv) {
  if (argc != 3) {
    print_usage(argv[0]); 

    exit(1);
  }

  int n;
  double *matrix;
  double *vector;

  srand48(time(NULL));

  n = atoi(argv[1]);

  matrix = malloc(n * n * sizeof(double));
  vector = malloc(n * sizeof(double));

  generate_matrix(matrix, n, n);
  generate_matrix(vector, 1, n);

  print_matrix(matrix, n, n);
  printf("\n");
  print_matrix(vector, 1, n);

  free(matrix);
  free(vector);

  return 0;
}

void generate_matrix(double *matrix, int n, int m) {
  int i, j;

  for (i = 0; i < n; ++i) {
    for (j = 0; j < m; ++j) {
      matrix[i * n + j] = drand48() * MAX_RAND;
    }
  }
}

void print_matrix(double *matrix, int n, int m) {
  int i, j;

  for (i = 0; i < n; ++i) {
    for (j = 0; j < m; ++j) {
      printf("%.2f\t", matrix[i * n + j]);
    }

    printf("\n");
  }
}

void print_usage(char *prog) {
  printf("%s n threads\n", prog); 
  printf("n - Size of the matrix\n");
  printf("threads - Number of threads to run on\n");
}
