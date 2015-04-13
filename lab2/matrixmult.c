#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

void print_usage(char *prog, char *error) {
  printf("%s <form> <flag> <n> <A> <B>\n", prog);
  printf("\t%s\n\n", error);
  printf("\t<form> ijk, ikj, or kij\n");
  printf("\t<flag> R (Random maxtrices), or I (Input matrices)\n");
  printf("\t<n> size of the matrix n x n\n");
  printf("\t<A> n x n matrix (only required if flag == I)\n");
  printf("\t<B> n x n matrix (only required if flag == I)\n");
}

int main(int argc, char **argv) {
  int comm_sz, comm_rank;

	MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  int n;
  char flag;
  char *form;
  int *matrix_a;
  int *matrix_b;
  
  form = malloc(3 * sizeof(char));
 
  scanf("%s", form);
  scanf("%s", &flag);
  scanf("%d", &n);

  if (flag == 'R') {
    printf("Random\n"); 
  } else if (flag == 'I') {
    printf("Input\n");
  } else {
    printf("Bad input\n");
  }  

  free(form);

	MPI_Finalize();

	return 0;
}
