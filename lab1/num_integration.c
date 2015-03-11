#include <math.h>
#include <stdio.h>

#include <mpi.h>

double calc_trap(double a, double b, int n, double h);
double func(double x);
void get_user_input(int rank, double *a, double *b, int *n);

int main(int argc, char **argv) {
	double a, b, h;
	double local_a, local_b, local_int, total_int;
	int n, local_n;
	int comm_sz, comm_rank;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

	get_user_input(comm_rank, &a, &b, &n);	

	h = (b-a)/n;

	local_n = n/comm_sz;	

	local_a = a+comm_rank*local_n*h;
		
	local_b = local_a+local_n*h;

	local_int = calc_trap(local_a, local_b, local_n, h);

	MPI_Reduce(&local_int, &total_int, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (comm_rank == 0) {
		printf("Running on %d processors.\n", comm_sz);
		printf("Elapsed time = %e seconds\n", 0.0);
		printf("With n = %d trapezoids,\n", n);
		printf("our estimate of the integral from %lf to %lf = %.13e\n", a, b, total_int);
		printf("absolute relative true error = %e is less than criteria = %e\n", 0.0, 0.0);
	}	

	MPI_Finalize();

	return 0;
}

double calc_trap(double a, double b, int n, double h) {
	int i;
	double approx = (func(a)+func(b))/2;

	for (i = 1; i <= n-1; ++i) {
		approx += func(a+i*h);
	}

	return h*approx;
}

double func(double x) {
	return 8+5*sin(x/4)-2*cos(x/5)+cos(x/3);
}

void get_user_input(int rank, double *a, double *b, int *n) {
	if (rank == 0) {
		printf("Enter a, b, and n\n");

		scanf("%lf%lf%d", a, b, n);
	}

	MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
