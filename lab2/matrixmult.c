#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#define MAX_RAND 100

#define DEBUG(r, f, ...) fprintf(stdout,"Rank: %i : " f, r, ##__VA_ARGS__);

void get_user_input(int rank, char form[], char flag[], int *n);
void get_random_matrix(int *matrix_a, int * matrix_b, int n);
void get_user_matrix(int *matrix_a, int * matrix_b, int n);
void print_usage(char *prog, char *error);
void print_matrix(int rank, int *matrix, int n);

int main(int argc, char **argv) {
  int comm_sz, comm_rank;

	MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  int n, i, j, l;
  char flag[2];
  char form[4];
  int *matrix_a;
  int *matrix_b;
  int local_n;

  get_user_input(comm_rank, form, flag, &n);

  matrix_a = malloc(n * n * sizeof(int));
  matrix_b = malloc(n * n * sizeof(int));
  
  local_n = (n * n)/ comm_sz;

  printf("%i\n", local_n);
  
  if (comm_rank == 0) {
    if (strcmp(flag, "I") == 0) {
      get_user_matrix(matrix_a, matrix_b, n);
    } else if (strcmp(flag, "R") == 0) {
      get_random_matrix(matrix_a, matrix_b, n);
    } else {
      printf("Invalid value for <flag> '%s'\n", flag);

      MPI_Finalize();

      return 1;
    }
  }

  if (strcmp(form, "ijk") == 0) {

  } else if (strcmp(form, "ikj") == 0) {

  } else if (strcmp(form, "kij") == 0) {

  } else {
    printf("Invalid value for <form> '%s'\n", form);
  }

  free(matrix_a);
  free(matrix_b);

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

void get_random_matrix(int *matrix_a, int *matrix_b, int n) {
  int i, j;

  srand(time(NULL));

  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      matrix_a[i * n + j] = rand() % MAX_RAND;
      matrix_b[i * n + j] = rand() % MAX_RAND;
    }
  }
}

void get_user_matrix(int *matrix_a, int *matrix_b, int n) {
  int i, j;

  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      scanf("%i", &matrix_a[i * n + j]);
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
