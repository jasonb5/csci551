#include <math.h>
#include <stdio.h>

#define ABS(x)(((x)<0)?-(x):(x))
#define MAX_ERROR (0.5*pow(10, -12))

double calc_approx(double a, double b, int traps);
double calc_actual(double a, double b);
double func(double x);
double func_int(double x);

int main(int argc, char **argv) {
	int stride, index;
	double a, b, actual, approx, true_error;

	printf("Enter a and b\n");
	
	scanf("%lf%lf", &a, &b);

	actual = calc_actual(a, b);

	stride = 4096;

	// Repeat until we cannot get a smaller stride
	while (1) {
		if (stride == 1) {
      index += stride*2;

      break;
    }

		// Using current stride find an approximation
		// fitting our requirements of a absolute
		// relative true error less than equal to 
		// MAX_ERROR
		while (1) {
			index += stride;

			approx = calc_approx(a, b, index);

			true_error = ABS(actual-approx)*100/actual;

			if (true_error <= MAX_ERROR) {
				break;
			}
		}

		// Found a valid value step value back by strid
		index -= stride;

		// Reduce stride by 2 to repeat process
		stride /= 2;
	}

  printf("traps %d\n", index);
  printf("actual %.10lf\n", actual);
  printf("approx %.10lf\n", approx);
  printf("true error %.16e\n", true_error);

	return 0;
}

// Approximates integral using Trapezoidal Method
double calc_approx(double a, double b, int n) {
	int i;
	double h = (b-a)/n;
	double approx = (func(a)+func(b))/2;

	for (i = 1; i <= n-1; ++i) {
		approx += func(a+i*h);		
	}	

	return h * approx;
}

// Calculates definite integral
double calc_actual(double a, double b) {
	return func_int(b)-func_int(a);
}

// Original function
double func(double x) {
	return 8+5*sin(x/4)-2*cos(x/5)+cos(x/3);
}

// Indefinite integral function
double func_int(double x) {
	return 8*x-10*sin(x/5)+3*sin(x/3)-20*cos(x/4);
}
