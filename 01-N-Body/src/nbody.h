#include <stdio.h>
#include <stdlib.h>

/*
 * Vector calculation
 *
 */

typedef struct {
  double x;
  double y;
  double z;
} vector;

vector vector_add(const vector v1, const vector v2);
void   vector_add_to(vector* v1, const vector v2);

vector vector_subtract(const vector v1, const vector v2);
void   vector_subtract_from(vector* v1, const vector v2);

vector scalar_mult(double a, const vector v);
void   scalar_mult_to(double a, vector* v);

/*
 * Project stuff
 *
 */

typedef struct {
  double mass;
  vector position;
  vector velocity;
} object;

int main ();

/*
 * IO - Reading & writing
 *
 */

object* read_data(const char* file, int* N, double* t_max, double* delta_t);

void print_object(const object o);
