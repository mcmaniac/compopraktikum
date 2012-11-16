#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct {
  double U;
  double I;
  double delta_I;
} measurement;

measurement* read_measurements(const char* fp, int *N);
double chi_square_(int N, measurement *m, double a, double b);
void linear_regression(int N, measurement *m, double *a, double *b, double *sig_a, double *sig_b);
