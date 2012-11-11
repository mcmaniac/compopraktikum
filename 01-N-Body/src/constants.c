#include "main.h"

double total_energy(const data* dat, const vector* r, const vector* v)
{
  double E = 0;
  int i, j;
  for (i = 0; i < dat->N; i++)
  {
    double m = dat->objects[i].mass;
    // kinetic energy
    E += 0.5 * m * vector_mult(v[i], v[i]);
    for (j = i+1; j < dat->N; j++)
    {
      // potential energy
      E -= G * m * m / vector_abs( vector_diff(r[i], r[j]) );
    }
  }
  return E;
}

/*
double total_momentum(const data* dat, const vector* v)
{
  vector p = nullVector();
  int i;
  for (i = 0; i < dat->N; i++)
  {
    double m = dat->objects[i].mass;
    vector_add_to(&p, scalar_mult(m, v[i]));
  }
  return vector_abs(p);
}
*/

double total_angular_momentum(const data* dat, const vector* r, const vector* v)
{
  vector j = nullVector();
  int i;
  for (i = 0; i < dat->N; i++)
  {
    vector_add_to(&j, vector_cross_prod(r[i], v[i]));
  }
  return vector_abs(j);
}

/*
double total_center_of_mass(const data* dat)
{
  int i;
  vector com = nullVector();
  for (i = 0; i < dat->N; i++)
    vector_add_to(&com, dat->objects[i].position);
  return vector_abs(com);
}
*/

vector runge_lenz(const data* dat, const vector* r, const vector* v, vector* _j)
{
  if (dat->N != 2)
  {
    printf("Runge-Lenz-Vektor nur für N=2!");
    return nullVector();
  }
  vector r_rel = vector_diff(r[1], r[0]),
         v_rel = vector_diff(v[1], v[0]),
         j     = vector_cross_prod(r_rel, v_rel),
         e     = vector_diff(scalar_mult(1.0/(G*M), vector_cross_prod(v_rel,j)),
                             scalar_mult(1.0/vector_abs(r_rel), r_rel));
  if (_j)
    *_j = j;
  return e;
}

double total_runge_lenz(const data* dat, const vector* r, const vector* v)
{
  return vector_abs(runge_lenz(dat, r, v, NULL));
}

double semimajor_axis(const data* dat, const vector* r, const vector* v)
{
  vector j = nullVector(),
         e = runge_lenz(dat, r, v, &j);
  return vector_mult(j,j) / (G*M) / (1 - vector_mult(e,e));
}
