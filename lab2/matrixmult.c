#include <stdio.h>

#include <mpi.h>

int main(int argc, char **argv) {
	if (argc < 4) {
		return 1;
	}

	MPI_Init(&argc, &argv);

	MPI_Finalize();

	return 0;
}
