#include "vector.h"

vector null_vector(int N)
{
  vector v = {
    .N   = N,
    .val = (double*) malloc(N*sizeof(double))
  };
  int i;
  for (i = 0; i < N; i++)
    VectorSET(v, i, 0);
  return v;
}

vector vector_copy(const vector v)
{
  vector v_c = {
    .N   = v.N,
    .val = (double*) malloc(v.N*sizeof(double))
  };
  int i;
  for (i = 0; i < v.N; i++)
    VectorSET(v_c, i, VectorGET(v, i));
  return v_c;
}

void vector_destroy(vector v)
{
  if (v.val)
    free(v.val);
}

/*
vector vector_add(const vector v1, const vector v2)
{
  vector v_res = null_vector(v1.N);
  int i;
  for (i = 0; i < v1.N; i++)
    VectorSET(v_res, i, VectorGET(v1, i) + VectorGET(v2, i));
  return v_res;
}
*/

void vector_add(vector v1, const vector v2)
{
  int i;
  for (i = 0; i < v1.N; i++)
    VectorSET(v1, i, VectorGET(v1, i) + VectorGET(v2, i));
}

/*
vector scalar_mult(double a, const vector v)
{
  vector v_res = null_vector(v.N);
  int i;
  for (i = 0; i < v.N; i++)
    VectorSET(v_res, i, a * VectorGET(v, i));
  return v_res;
}
*/

void scalar_mult(double a, vector v)
{
  int i;
  for (i = 0; i < v.N; i++)
    VectorSET(v, i, a * VectorGET(v, i));
}

/*
vector vector_diff(const vector v1, const vector v2)
{
  return vector_add(v1, scalar_mult(-1, v2));
}
*/

void vector_diff(vector v1, const vector v2)
{
  vector neg = vector_copy(v2);
  scalar_mult(-1, neg);
  vector_add(v1, neg);
  vector_destroy(neg);
}

double vector_abs(const vector v)
{
  double sum = 0;
  int i;
  for (i = 0; i < v.N; i++)
    sum += pow(VectorGET(v, i),2);
  return sqrt(sum);
}

double vector_mult(const vector v1, const vector v2)
{
  double sum = 0;
  int i;
  for (i = 0; i < v1.N; i++)
    sum += VectorGET(v1, i) * VectorGET(v2, i);
  return sum;
}

vector vector_cross_prod(const vector v1, const vector v2)
{
  vector v_res = null_vector(3);
  VectorX(v_res) = VectorY(v1) * VectorZ(v2) - VectorZ(v1) * VectorY(v2);
  VectorY(v_res) = VectorZ(v1) * VectorX(v2) - VectorX(v1) * VectorZ(v2);
  VectorZ(v_res) = VectorX(v1) * VectorY(v2) - VectorY(v1) * VectorX(v2);
  return v_res;
}

void vector_swap(vector v, int k, int l)
{
  double vk = VectorGET(v, k);
  VectorSET(v, k, VectorGET(v, l));
  VectorSET(v, l, vk);
}
