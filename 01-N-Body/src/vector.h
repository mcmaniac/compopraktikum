#pragma once
#include <math.h>

/*
 * Vector calculation
 *
 */

typedef struct {
  double x;
  double y;
  double z;
} vector;

vector nullVector();

vector vector_add(const vector v1, const vector v2);
void   vector_add_to(vector* v1, const vector v2);

vector vector_diff(const vector v1, const vector v2);
void   vector_diff_from(vector* v1, const vector v2);

vector vector_cross_prod(const vector v1, const vector v2);
double vector_mult(const vector v1, const vector v2);

vector scalar_mult(double a, const vector v);
void   scalar_mult_to(double a, vector* v);

double vector_abs(const vector v);
