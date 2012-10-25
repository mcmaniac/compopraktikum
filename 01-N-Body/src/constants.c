#include "main.h"

double total_energy(const data* dat)
{
  double E = 0;
  int i, j;
  for (i = 0; i < dat->N; i++)
  {
    object o = dat->objects[i];
    // kinetic energy
    E += 0.5 * o.mass * vector_mult(o.velocity, o.velocity);
    for (j = i+1; j < dat->N; j++)
    {
      object u = dat->objects[j];
      // potential energy
      E -= G * o.mass * u.mass / vector_abs( vector_diff(o.position, u.position) );
    }
  }
  return E;
}

double total_momentum(const data* dat)
{
  vector p = nullVector();
  int i;
  for (i = 0; i < dat->N; i++)
  {
    object o = dat->objects[i];
    vector_add_to(&p, scalar_mult(o.mass, o.velocity));
  }
  return vector_abs(p);
}

double total_angular_momentum(const data* dat)
{
  vector j = nullVector();
  int i;
  for (i = 0; i < dat->N; i++)
  {
    object o = dat->objects[i];
    vector_add_to(&j, vector_cross_prod(o.position, o.velocity));
  }
  return vector_abs(j);
}

double total_center_of_mass(const data* dat)
{
  int i;
  vector com = nullVector();
  for (i = 0; i < dat->N; i++)
    vector_add_to(&com, dat->objects[i].position);
  return vector_abs(com);
}

vector runge_lenz(const data* dat, vector* _j)
{
  if (dat->N != 2)
  {
    printf("Runge-Lenz-Vektor nur für N=2!");
    return nullVector();
  }
  object o1 = dat->objects[0],
         o2 = dat->objects[1];
  vector r = vector_diff(o2.position, o1.position),
         v = vector_diff(o2.velocity, o1.velocity),
         j = vector_cross_prod(r, v),
         e = vector_diff(scalar_mult(1.0/(G*M), vector_cross_prod(v,j)),
                         scalar_mult(1.0/vector_abs(r), r));
  if (_j)
    *_j = j;
  return e;
}
double total_runge_lenz(const data* dat)
{
  return vector_abs(runge_lenz(dat, NULL));
}

double semimajor_axis(const data* dat)
{
  vector j = nullVector(),
         e = runge_lenz(dat, &j);
  return vector_mult(j,j) / (G*M) / (1 - vector_mult(e,e));
}
