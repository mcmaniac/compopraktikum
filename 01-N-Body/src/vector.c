#include "vector.h"

vector nullVector()
{
  vector v = {0, 0, 0};
  return v;
}

vector vector_add(const vector v1, const vector v2)
{
  vector vres = {
    .x = v1.x + v2.x,
    .y = v1.y + v2.y,
    .z = v1.z + v2.z
  };
  return vres;
}

void vector_add_to(vector* v1, const vector v2)
{
  v1->x += v2.x;
  v1->y += v2.y;
  v1->z += v2.z;
}

vector scalar_mult(double a, const vector v)
{
  vector vres = {
    .x = a*v.x,
    .y = a*v.y,
    .z = a*v.z
  };
  return vres;
}

void scalar_mult_to(double a, vector* v)
{
  v->x = v->x * a;
  v->y = v->y * a;
  v->z = v->z * a;
}

vector vector_diff(const vector v1, const vector v2)
{
  return vector_add(v1, scalar_mult(-1, v2));
}

void vector_diff_from(vector* v1, const vector v2)
{
  vector_add_to(v1, scalar_mult(-1, v2));
}

double vector_abs(const vector v)
{
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

double vector_mult(const vector v1, const vector v2)
{
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}
vector vector_cross_prod(const vector v1, const vector v2)
{
  vector v;
  v.x = v1.y * v2.z - v1.z * v2.y;
  v.y = v1.z * v2.x - v1.x * v2.z;
  v.z = v1.x * v2.y - v1.y * v2.x;
  return v;
}
