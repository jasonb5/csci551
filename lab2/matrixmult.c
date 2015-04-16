#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#define MAX_RAND 100

#define DEBUG(r, f, ...) fprintf(stdout,"Rank: %i : " f, r, ##__VA_ARGS__);

#define TRUE 1
#define FALSE 0

void transpose(int **matrix, int n);
void get_user_input(int rank, char form[], char flag[], int *n);
void get_random_matrix(int *matrix_a, int * matrix_b, int n, int transpose);
void get_user_matrix(int *matrix_a, int * matrix_b, int n, int transpose);
void print_usage(char *prog, char *error);
void print_matrix(int rank, int *matrix, int n);

int main(int argc, char **argv) {
  int comm_sz, comm_rank;

	MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  int n, i, j, k;
  char flag[2];
  char form[4];
  int *matrix_a;
  int *matrix_b;
  int *matrix_r;
  int local_n;
  int sum;

  get_user_input(comm_rank, form, flag, &n);

  matrix_a = malloc(n * n * sizeof(int));
  matrix_b = malloc(n * n * sizeof(int));
  matrix_r = malloc(n * n * sizeof(int)); 

  memset(matrix_a, 0, n * n * sizeof(int));
  memset(matrix_b, 0, n * n * sizeof(int));
  memset(matrix_r, 0, n * n * sizeof(int));
  
  local_n = (n * n)/ comm_sz;
  
  if (comm_rank == 0) {
    int transpose = FALSE;

    if (strcmp(form, "kij") == 0) {
      transpose = TRUE;
    }

    if (strcmp(flag, "I") == 0) {
      get_user_matrix(matrix_a, matrix_b, n, transpose);
    } else if (strcmp(flag, "R") == 0) {
      get_random_matrix(matrix_a, matrix_b, n, transpose);
    } else {
      printf("Invalid value for <flag> '%s'\n", flag);

      MPI_Finalize();

      return 1;
    }

    printf("Matrix A:\n");
    print_matrix(comm_rank, matrix_a, n);
    printf("\nMatrix B:\n");
    print_matrix(comm_rank, matrix_b, n);
    printf("\n");
  }

  //  a = [[0, 1, 2], b = [[0, 1],
  //       [3, 4, 5],      [2, 3],
  //       [6, 7, 8]]      [4, 5]]
  //  a = 3x3 = am x an b = 3x2 = bm x bn
  //
  //  ijk
  //  for i[0..am-1]
  //    for j[0..bn-1]
  //      for k[0..an-1]
  //        sum += a[i,k] * b[k,j]
  //
  //  ikj
  //  for i[0..am-1]
  //    for k[0..an-1]
  //      for j[0..bn-1]
  //        sum += a[i,k] * b[k,j] 
  //
  //  kij
  //  for k[0..an-1]
  //    for i[0..am-1]
  //      for j[0..bn-1]
  //        sum += a[i,k] * b[k,j]
  //

  MPI_Barrier(MPI_COMM_WORLD);

  if (strcmp(form, "ijk") == 0) {
    MPI_Scatter(matrix_a, local_n, MPI_INT, matrix_a, local_n, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(matrix_b, n * n, MPI_INT, 0, MPI_COMM_WORLD);

    for (i = 0; i < local_n / n; ++i) {
      for (j = 0; j < n; ++j) {
        sum = 0;

        for (k = 0; k < n; ++k) {
          sum += matrix_a[i * n + k] * matrix_b[k * n + j];    
        }

        matrix_r[i * n + j] = sum;
      }
    }

    MPI_Gather(matrix_r, local_n, MPI_INT, matrix_r, local_n, MPI_INT, 0, MPI_COMM_WORLD);
  } else if (strcmp(form, "ikj") == 0) {
    MPI_Scatter(matrix_a, local_n, MPI_INT, matrix_a, local_n, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(matrix_b, n * n, MPI_INT, 0, MPI_COMM_WORLD);

    for (i = 0; i < local_n / n; ++i) {
      for (k = 0; k < n; ++k) {
        for (j = 0; j < n; ++j) {
          matrix_r[i * n + j] += matrix_a[i * n + k] * matrix_b[k * n + j];    
        }
      }
    }
    
    MPI_Gather(matrix_r, local_n, MPI_INT, matrix_r, local_n, MPI_INT, 0, MPI_COMM_WORLD);
  } else if (strcmp(form, "kij") == 0) {
    MPI_Scatter(matrix_a, local_n, MPI_INT, matrix_a, local_n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(matrix_b, local_n, MPI_INT, matrix_b, local_n, MPI_INT, 0, MPI_COMM_WORLD);

    for (k = 0; k < n; ++k) {
      for (i = 0; i < local_n / n; ++i) {
        sum = 0;

        for (j = 0; j < n; ++j) {
          sum += matrix_a[i * n + j] * matrix_b[k * n + j]; 
        }

        matrix_r[i * n + k] = sum;
      }
    }

    MPI_Gather(matrix_r, local_n, MPI_INT, matrix_r, local_n, MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    printf("Invalid value for <form> '%s'\n", form);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_rank == 0) {
    print_matrix(comm_rank, matrix_r, n);
  }

  free(matrix_a);
  free(matrix_b);
  free(matrix_r);

	MPI_Finalize();

	return 0;
}

void get_user_input(int rank, char form[], char flag[], int *n) {
  if (rank == 0) {
    scanf("%s", form);
    scanf("%s", flag);
    scanf("%i", n);
  }

  MPI_Bcast(form, 4, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(n, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void get_random_matrix(int *matrix_a, int *matrix_b, int n, int transpose) {
  int i, j;

  srand(time(NULL));

  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      if (transpose) {
        matrix_a[j * n + i] = rand() % MAX_RAND;
      } else {
        matrix_a[i * n + j] = rand() % MAX_RAND;
      }

      matrix_b[i * n + j] = rand() % MAX_RAND;
    }
  }
}

void get_user_matrix(int *matrix_a, int *matrix_b, int n, int transpose) {
  int i, j;

  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      if (transpose) {
        scanf("%i", &matrix_a[j * n + i]);
      } else {
        scanf("%i", &matrix_a[i * n + j]);
      }
    }
  }

  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      scanf("%i", &matrix_b[i * n + j]);
    }
  }
}

void print_usage(char *prog, char *error) {
  printf("%s <form> <flag> <n> <A> <B>\n", prog);
  printf("\t%s\n\n", error);
  printf("\t<form> ijk, ikj, or kij\n");
  printf("\t<flag> R (Random maxtrices), or I (Input matrices)\n");
  printf("\t<n> size of the matrix n x n\n");
  printf("\t<A> n x n matrix (only required if flag == I)\n");
  printf("\t<B> n x n matrix (only required if flag == I)\n");
}

void print_matrix(int rank, int *matrix, int n) {
  int x, y;

  for (x = 0; x < n; ++x) {
    for (y = 0; y < n; ++y) {
      DEBUG(rank, "%i\t", matrix[x * n + y]);
    }

    printf("\n");
  }
}
