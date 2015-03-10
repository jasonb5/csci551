#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#define SPREAD 4
#define ABS(x)(((x)<0)?-(x):(x))
#define MAX_ERROR 0.5*pow(10, -12)

int search(int traps, int tmin, int tmax, double a, double b);
double calc_true_error(double true_val, double approx_val);
double calc_actual(double a, double b);
double calc_approx(double a, double b, int traps, double h);
double func(double x);
double func_int(double x);

double current_min;

int main(int argc, char **argv) {
  int traps, mid, optimal_traps;
  double lend, rend, approx, actual;

  printf("Enter left, right, and trapezoids\n");

  scanf("%lf%lf%d", &lend, &rend, &traps);

  mid = (traps/2);

  optimal_traps = search(mid, 0, traps, lend, rend);

  printf("t = %d\n", optimal_traps);

  approx = calc_approx(lend, rend, optimal_traps, (rend-lend)/optimal_traps); 

  actual = calc_actual(lend, rend);

  printf("absolute relative true error %.16lf\n", calc_true_error(actual, approx));

	return 0;
}

int search(int traps, int tmin, int tmax, double a, double b) {
  int x, lvalue, rvalue, lmid, rmid, min_index;
  double h, approx, actual, true_error, min_value;

  if (traps <= tmin || traps > tmax) return traps;

  min_value = DBL_MAX;

  for (x = traps-SPREAD; x < traps+SPREAD+1; ++x) {
    h = (b-a)/x;

    approx = calc_approx(a, b, x, h); 

    actual = calc_actual(a, b);

    true_error = calc_true_error(actual, approx);

    if (true_error <= MAX_ERROR) {
      printf("traps %d\n", x);
      printf("h %.16lf\n", h);
      printf("approx %.16lf\n", approx);
      printf("actual %.16lf\n", actual);
      printf("true error %.13lf\n", true_error);
      printf("true error %.13lf\n", MAX_ERROR);

      return x;
    }

    if (true_error < min_value) {
      min_index = x;

      min_value = true_error;
    }

    printf("traps %d\n", x);
    printf("h %.16lf\n", h);
    printf("approx %.16lf\n", approx);
    printf("actual %.16lf\n", actual);
    printf("true error %.13lf\n", true_error);
    printf("true error %.13lf\n", MAX_ERROR);
  }

  printf("Done looking at values found min index %d\n", min_index);

  if (min_value < current_min) {
    current_min = min_value;
  } else {
    min_index = traps+1;
  }
  
  if (min_index <= traps) {
    lvalue = floor((traps-tmin)/2);
    lmid = tmin+lvalue;

    printf("We're going left\n");
    return search(lmid, tmin, traps, a, b);
  } else {
    rvalue = floor((tmax-traps)/2);
    rmid = traps+rvalue;

    printf("We're going right\n");
    return search(rmid, traps, tmax, a, b);
  }
}

double calc_true_error(double true_val, double approx_val) {
  return ABS(true_val-approx_val)/true_val;
}

double calc_actual(double a, double b) {
  return func_int(b)-func_int(a);
}

double calc_approx(double a, double b, int traps, double h) {
  int x;
  double approx = 0.0;

  for (x = 1; x < traps-1; ++x) {
    approx += func(a+x*h);
  }

  approx *= 2;
  approx += func(a)+func(b);
  approx *= h/2;

  return approx;
}

double func(double x) {
  return cos(x/3)-2*cos(x/5)+5*sin(x/4)+8;
}

double func_int(double x) {
  return 8*x-10*sin(x/5)+3*sin(x/3)-20*cos(x/4);
}
