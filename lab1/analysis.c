#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#define DEBUG

#ifdef DEBUG
#define Debug(x, ...) fprintf(stdout, x, ##__VA_ARGS__)
#else
#define Debug(x, ...)
#endif

#define DIVISIONS 20
#define SPREAD 2
#define MAX_TRUE_ERROR (0.5*pow(10,2-14))/100
#define ABS(x) (x<0) ? -x : x

int search(int value, int range, unsigned long int mid, unsigned long int min, unsigned long int max);
double true_error(int value, double h);
double func_def(double value);
double deriv_def(double value);

int main(int argc, char **argv) {
	int x, min, max, index, value;
	unsigned long int traps, mid, step, avg;

	if (argc != 4) {
		printf("Usage: %s min max trapezoids\n", argv[0]);

		exit(1);
	}

	min = atoi(argv[1]);
	max = atoi(argv[2]);
	traps = atol(argv[3]);

	mid = floor(traps/2);

	step = floor((max-min)/DIVISIONS);

	for (x = 0, avg = 0; x < DIVISIONS; ++x) {
		value = (int)min+x*step;

		index = search(value, max-min, mid, 0, traps);

		Debug("Returned %d trapezoids at value %d\n", index, value);

		avg += index;
	}

	Debug("Average Trapezoids %lu\n", (unsigned long int)floorl((long double)avg/DIVISIONS));
	Debug("Max True Error %.16f\n", MAX_TRUE_ERROR);
	
	return 0;
}

int search(int value, int range, unsigned long int mid, unsigned long int min, unsigned long int max) {
	double temp;
	unsigned long int x;
	double min_value = DBL_MAX;	
	int min_index = -1, child_index = -1;
	int lvalue, lmid, rvalue, rmid;

	for (x = mid-SPREAD; x < mid+SPREAD+1; ++x) {
		if (x < min || x > max) break;

		temp = true_error(value, (double)range/x);

		if (temp < min_value) {
			min_index = x;
			
			min_value = temp;
		}
	}

	if (min_value <= MAX_TRUE_ERROR) return min_index;

	lvalue = floor((mid-min)/2);
	rvalue = floor((max-mid)/2);

	lmid = min+lvalue;
	rmid = mid+rvalue;

	//Debug("Searching at %d trapezoids min %d max %d lvalue %d lmid %d rvalue %d rmid %d\n", mid, min, max, lvalue, lmid, rvalue, rmid);

	if (lvalue > 0 && min_index <= mid) {
		child_index = search(value, range, lmid, min, mid);
	}

	if (rvalue > 0 && min_index > mid) {
		child_index = search(value, range, rmid, mid, max);
	}	

	return (child_index == -1 && min_index != -1) ? min_index : child_index;
}

double true_error(int value, double h) {
	double val1 = func_def((double)value+h);
	double val2 = func_def((double)value);
	double approx_val = (val1-val2)/h;
	double true_val = deriv_def(value);

	//Debug("Approximate value %.8f\tTrue value %.8f\t True error %.16f\n", approx_val, true_val, (true_val-approx_val)/true_val);

	return ABS((true_val-approx_val)/true_val);
}

double func_def(double value) {
	return cos(value/3)-2*cos(value/5)+5*sin(value/4)+8;
}

double deriv_def(double value) {
	return ((double)2/5)*sin(value/5)-((double)1/3)*sin(value/3)+((double)5/4)*cos(value/4);
}
