#pragma once
#include <stdlib.h>
#include <math.h>

/*
 * Vector calculation
 *
 */

typedef struct {
  double N;
  double* val;
} vector;

#define VectorSET(V, i, v) ((V).val[(i)] = v)
#define VectorGET(V, i)    ((V).val[(i)])

#define VectorX(V) (VectorGET((V), 0))
#define VectorY(V) (VectorGET((V), 1))
#define VectorZ(V) (VectorGET((V), 2))

#define VectorADD(V, i, v) (VectorSET((V), (i), VectorGET((V), (i)) + (v)))
#define VectorSUB(V, i, v) (VectorSET((V), (i), VectorGET((V), (i)) - (v)))
#define VectorMLT(V, i, v) (VectorSET((V), (i), VectorGET((V), (i)) * (v)))
#define VectorDIV(V, i, v) (VectorSET((V), (i), VectorGET((V), (i)) / (v)))

vector null_vector(int N);
vector vector_copy(const vector v);

void vector_destroy(vector v);

//vector vector_add(const vector v1, const vector v2);
void   vector_add(vector v1, const vector v2);

//vector vector_diff(const vector v1, const vector v2);
void   vector_diff(vector v1, const vector v2);

//vector vector_cross_prod(const vector v1, const vector v2);
double vector_mult(const vector v1, const vector v2);

//vector scalar_mult(double a, const vector v);
void   scalar_mult(double a, vector v);

double vector_abs(const vector v);

void vector_swap(vector v, int k, int l);
