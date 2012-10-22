#include <stdio.h>
#include <stdlib.h>
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

vector vector_add(const vector v1, const vector v2);
void   vector_add_to(vector* v1, const vector v2);

vector vector_diff(const vector v1, const vector v2);
void   vector_diff_from(vector* v1, const vector v2);

double vector_mult(const vector v1, const vector v2);

vector scalar_mult(double a, const vector v);
void   scalar_mult_to(double a, vector* v);

double vector_abs(const vector v);

/*
 * Project stuff
 *
 */

typedef struct {
  double mass;
  vector position;
  vector velocity;
} object;

typedef struct {
  int N;
  double t_max;
  double eta;
  object* objects;
} data;

void free_data(data* dat);

// gravitation constant
double G;

int main ();

void runge_kutta();

/*
 * IO - Reading & writing
 *
 */

data* read_data(const char* file); // , int* N, double* t_max, double* delta_t);

void print_object(const object o);
