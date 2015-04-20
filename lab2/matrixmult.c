/**
 * @file matrixmult.c
 * @author Jason Boutte <jboutte@mail.csuchico.edu>
 *
 * @brief Parallelization of ijk, ikj, and kij dot product forms.
 */
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#define MAX_RAND 100 // Maximum value for random numbers

#define TRUE 1 // defines true
#define FALSE 0 // defines false

/**
 * @brief Gets user input for dot product. 
 *
 * @param rank Rank of current process
 * @param form ijk form to be used
 * @param flag I for input or R for random
 * @param n Stride size of matrix
 */ 
void get_user_input(int rank, char form[], char flag[], int *n);

/**
 * @brief Generates two matrices in row major format.
 *
 * @param matrix_a,matrix_b Pointers to arrays of length n^2
 * @param n Stride of the matrix
 * @param transpose Transpose matrix_a
 */
void get_random_matrix(int *matrix_a, int * matrix_b, int n, int transpose);

/**
 * @brief Gets user input of two matrices in row major format.
 *
 * @param matrix_a,matrix_b Pointers to arrays of length n^2
 * @param n Stride of the matrix
 * @param transpose Transpose matrix_a
 */
void get_user_matrix(int *matrix_a, int * matrix_b, int n, int transpose);

/**
 * @brief Prints a usage message.
 *
 * @param prog Program name
 * @param error Error message to be displayed
 */
void print_usage(char *prog, char *error);

/**
 * @brief Prints the contents of a matrix of size n^2.
 *
 * @param matrix Pointer to array of size n^2
 * @param n Size Stride size of matrix
 */
void print_matrix(int *matrix, int n);

int main(int argc, char **argv) {
  int comm_sz, comm_rank;

  // Initialize MPI
	MPI_Init(&argc, &argv);

  // Determine number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  // Determine current processes rank
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  int sum;
  int local_n;
  char flag[2];
  char form[4];
  int *matrix_a;
  int *matrix_b;
  int *matrix_r;
  int n, i, j, k;
  int send, value;
  int transpose = FALSE;
  double start, end;

  get_user_input(comm_rank, form, flag, &n);

  // Create the input matrices
  matrix_a = malloc(n * n * sizeof(int));
  matrix_b = malloc(n * n * sizeof(int));
  matrix_r = malloc(n * n * sizeof(int)); 

  // Zero all matrices
  memset(matrix_a, 0, n * n * sizeof(int));
  memset(matrix_b, 0, n * n * sizeof(int));
  memset(matrix_r, 0, n * n * sizeof(int));
 
  // Determine the local item count 
  local_n = (n * n)/ comm_sz;
  
  if (comm_rank == 0) {
    // Transpose matrix_a if for is kij
    if (strcmp(form, "kij") == 0) {
      transpose = TRUE;
    }

    // Create the input matrices
    if (strcmp(flag, "I") == 0) {
      get_user_matrix(matrix_a, matrix_b, n, transpose);
    } else if (strcmp(flag, "R") == 0) {
      get_random_matrix(matrix_a, matrix_b, n, transpose);
    } else {
      printf("Invalid value for <flag> '%s'\n", flag);

      MPI_Finalize();

      return 1;
    }

#ifdef DEBUG
    printf("Matrix A\n");
    print_matrix(matrix_a, n);
    printf("Matrix B\n");
    print_matrix(matrix_b, n);
    printf("\n");
#endif 
  }

  // Synchronize processes for timing
  MPI_Barrier(MPI_COMM_WORLD);

  // Begin timing on root process
  if (comm_rank == 0) {
    start = MPI_Wtime();
  }

  if (strcmp(form, "ijk") == 0) { 
    // ijk form
    //
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
    // ikj form
    //
    MPI_Scatter(matrix_a, local_n, MPI_INT, matrix_a, local_n, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(matrix_b, n * n, MPI_INT, 0, MPI_COMM_WORLD);

    for (i = 0; i < local_n / n; ++i) {
      for (k = 0; k < n; ++k) {
        value = matrix_a[i * n + k];
      
        for (j = 0; j < n; ++j) {
          matrix_r[i * n + j] += value * matrix_b[k * n + j];    
        }
      }
    }
    
    MPI_Gather(matrix_r, local_n, MPI_INT, matrix_r, local_n, MPI_INT, 0, MPI_COMM_WORLD);
  } else if (strcmp(form, "kij") == 0) {
    // kij form
    // 
    MPI_Scatter(matrix_a, local_n, MPI_INT, matrix_a, local_n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(matrix_b, local_n, MPI_INT, matrix_b, local_n, MPI_INT, 0, MPI_COMM_WORLD);

    for (k = 0; k < local_n / n; ++k) {
      for (i = 0; i < n; ++i) {
        value = matrix_a[k * n + i];

        for (j = 0; j < n; ++j) {
          matrix_r[i * n + j] += value * matrix_b[k * n + j]; 
        }
      }
    }

    for (i = 0; i < n; ++i) {
      for (j = 0; j < n; ++j) {
        send = matrix_r[i * n + j];

        MPI_Reduce(&send, &matrix_r[i * n + j], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      }
    }

  } else {
    printf("Invalid value for <form> '%s'\n", form);
  }

  // Print results
  if (comm_rank == 0) {
    end = MPI_Wtime();

    printf("running on %i processors\n", comm_sz);
    printf("elapsed time = %f seconds\n", end - start); 

#ifdef DEBUG
    print_matrix(matrix_r, n);
#endif

    if (strcmp(flag, "I") == 0) {
      print_matrix(matrix_r, n);
    }
  }

  // Cleanup resources
  free(matrix_a);
  free(matrix_b);
  free(matrix_r);

	MPI_Finalize();

	return 0;
}

/**
 * @brieft Gets the users input.
 *
 * We get the needed information from the user and 
 * broadcast it to all processes. This is not included
 * in the timing. Valid values for @p form are ijk,
 * ikj, and kij. Valid values for @p flag are I for
 * input and R for random.
 */
void get_user_input(int rank, char form[], char flag[], int *n) {
  if (rank == 0) {
    scanf("%s", form);
    scanf("%s", flag);
    scanf("%i", n);
  }

  MPI_Bcast(form, 4, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(n, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

/**
 * @brief Generates two arrays of length n^2.
 *
 * We seed the random generated and populate two
 * arrays of size n^2 with ranom values between
 * 0 and MAX_RAND
 */
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

/**
 * @brief Read in two user input arrays length n^2.
 *
 * Read in two matrices nxn.
 */
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

/**
 * @brief Prints usage message
 */
void print_usage(char *prog, char *error) {
  printf("%s <form> <flag> <n> <A> <B>\n", prog);
  printf("\t%s\n\n", error);
  printf("\t<form> ijk, ikj, or kij\n");
  printf("\t<flag> R (Random maxtrices), or I (Input matrices)\n");
  printf("\t<n> size of the matrix n x n\n");
  printf("\t<A> n x n matrix (only required if flag == I)\n");
  printf("\t<B> n x n matrix (only required if flag == I)\n");
}

/**
 * @brief Prints contents of matrix
 */
void print_matrix(int *matrix, int n) {
  int x, y;

  for (x = 0; x < n; ++x) {
    for (y = 0; y < n; ++y) {
      printf("%i ", matrix[x * n + y]);
    }

    printf("\n");
  }
}
