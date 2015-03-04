#include <stdio.h>

#include <mpi.h>

int main(int argc, char **argv) {
	int comm_sz, comm_rank;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

	MPI_Finalize();

	return 0;
}
