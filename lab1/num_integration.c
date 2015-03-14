/*
 * This code most comes from 
 * "An Introduction to Parallel Programming" 
 * By Pacheco
 *
 * Jason Boutte
 * CSCI 551
 * Lab #1
 * Sprin 2015 
 */

#include <math.h>
#include <stdio.h>

#include <mpi.h>

// Helper function for absolute values
#define ABS(x)(((x)<0)?-(x):(x))

// Maximum absolute relative true error for 
// 14 significant digits
#define MAX_ERROR 0.5*pow(10, -12)

// Function prototypes
double calc_rel_te(double actual, double approx);
double calc_trap(double a, double b, int n, double h);
double calc_actual(double a, double b);
double func(double x);
double func_int(double x);
void get_user_input(double *a, double *b, int *n);

int main(int argc, char **argv) {
	double a, b, h, rel_te;
	double local_a, local_b, local_int, total_int;
  double time_start, time_end;
	int n, local_n;
	int comm_sz, comm_rank;

	// Initialize MPI
	MPI_Init(&argc, &argv);

	// Get processes
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

	// Get current process rank
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

	// If root process get user input
  if (comm_rank == 0) {
    get_user_input(&a, &b, &n);	
  }

	// Synchronize all participating processes
  MPI_Barrier(MPI_COMM_WORLD);

	// Starting timing on root processes
  if (comm_rank == 0) {
    time_start = MPI_Wtime();
  }

	// Broadcast three values to all processes
	MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Calculate h value, this value is the same for 
	// all processes
	h = (b-a)/n;

	// Calculate our local n
	local_n = n/comm_sz;	

	// Calculate our starting endpoint
	local_a = a+comm_rank*local_n*h;
	
	// Calculate our ending endpoint	
	local_b = local_a+local_n*h;

	// Calculate our trapezoid
	local_int = calc_trap(local_a, local_b, local_n, h);

	// Send results to root process to be summed
	MPI_Reduce(&local_int, &total_int, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	// Root process prints results
	if (comm_rank == 0) {
    time_end = MPI_Wtime();

    rel_te = calc_rel_te(calc_actual(a, b), total_int);

		printf("Running on %d processors.\n", comm_sz);
		printf("Elapsed time = %e seconds\n", time_end-time_start);
		printf("With n = %d trapezoids,\n", n);
		printf("our estimate of the integral from %lf to %lf = %.13e\n", a, b, total_int);
		printf("absolute relative true error = %e %sis less than criteria = %e\n", rel_te, ((rel_te > MAX_ERROR) ? "NOT " : ""),  MAX_ERROR);
	}	

	// Clean up MPI
	MPI_Finalize();

	return 0;
}

// Calculates absolute relative true error as a percent
double calc_rel_te(double actual, double approx) {
  return ABS(actual-approx)*100/actual;
}

// Approximates integral using Trapezoidal Method
double calc_trap(double a, double b, int n, double h) {
	int i;
	double approx = (func(a)+func(b))/2;

	for (i = 1; i <= n-1; ++i) {
		approx += func(a+i*h);
	}

	return h*approx;
}

// Calculated definite integral
double calc_actual(double a, double b) {
  return func_int(b)-func_int(a);
}

// Original function
double func(double x) {
	return 8+5*sin(x/4)-2*cos(x/5)+cos(x/3);
}

// Indefinite integral
double func_int(double x) {
  return 8*x-10*sin(x/5)+3*sin(x/3)-20*cos(x/4);
}

// Get user input values for a, b, and n
void get_user_input(double *a, double *b, int *n) {
  printf("Enter a, b, and n\n");

  scanf("%lf%lf%d", a, b, n);
}
