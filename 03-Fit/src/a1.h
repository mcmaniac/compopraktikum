#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct {
  double val;
  double err;
} data;

double average_value(int N, data *dat);
double chi_square(int N, data *dat);
double sigma_square_intern(int N, data *dat);
double sigma_square_extern(int M, int N, data *dat);
