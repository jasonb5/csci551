#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_RAND 1.0e2 //1.0e6

void gen_augmented_matrix(double **matrix, int n);
void print_usage(char *prog);

int main(int argc, char **argv) {
  if (argc != 3) {
    print_usage(argv[0]); 

    exit(1);
  }

  int n, i, j, k;
  double **matrix;

  srand(time(NULL));
  srand48(time(NULL));

  n = atoi(argv[1]);

  matrix = malloc(n * sizeof(double));

  for (i = 0; i < n; ++i) {
    matrix[i] = malloc((n + 1) * sizeof(double));
  }

  gen_augmented_matrix(matrix, n);

  for (i = 0; i < n; ++i) {
    free(matrix[i]);
  }

  free(matrix);

  return 0;
}

void gen_augmented_matrix(double **matrix, int n) {
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

void print_usage(char *prog) {
  printf("%s n threads\n", prog); 
  printf("n - Size of the matrix\n");
  printf("threads - Number of threads to run on\n");
}
