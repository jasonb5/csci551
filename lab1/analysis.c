#include <math.h>
#include <stdio.h>

#ifdef DEBUG
#define debugf(x, ...) fprintf(stdout, x, ##__VA_ARGS__);
#else
#define debugf(x, ...)
#endif

#define ABS(x)(((x)<0)?-(x):(x))
#define MAX_ERROR 0.5*pow(10, -12)

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

	while (1) {
		if (stride == 1) break;

		printf("Index at %d\n", index);
		printf("Stride at %d\n", stride);

		while (1) {
			index += stride;

			approx = calc_approx(a, b, index);

			true_error = ABS(actual-approx)*100/actual;

			printf("traps %d\n", index);
			printf("actual %.10lf\n", actual);
			printf("approx %.10lf\n", approx);
			printf("true error %.16e\n", true_error);

			if (true_error <= MAX_ERROR) {
				printf("Found true error less than equal %.16e\n", MAX_ERROR);
		
				break;
			}
		}

		index -= stride;

		stride /= 2;
	}

	return 0;
}

double calc_approx(double a, double b, int n) {
	int i;
	double h = (b-a)/n;
	double approx = (func(a)+func(b))/2;

	for (i = 1; i <= n-1; ++i) {
		approx += func(a+i*h);		
	}	

	return h * approx;
}

double calc_actual(double a, double b) {
	return func_int(b)-func_int(a);
}

double func(double x) {
	return 8+5*sin(x/4)-2*cos(x/5)+cos(x/3);
}

double func_int(double x) {
	return 8*x-10*sin(x/5)+3*sin(x/3)-20*cos(x/4);
}
