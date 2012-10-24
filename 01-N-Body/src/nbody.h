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
// global configuration variables
const char* input;
const char* output;
FILE* file;

int main ();

void runge_kutta(data* dat, void (*output)(double time, double delta_t, const data* dat));

/*
 * Conserved quantities
 *
 */

double total_energy(const data* dat);
double total_momentum(const data* dat);
double total_angular_momentum(const data* dat);
double total_center_of_mass(const data* dat);
double total_runge_lenz(const data* dat);

/*
 * IO - Reading & writing
 *
 */

data* read_data(const char* file); // , int* N, double* t_max, double* delta_t);

void print_object(const object o);
void print_constants_to_file(double time, double delta_t, const data* dat);
