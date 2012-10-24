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

double total_runge_lenz(const data* dat)
{
  int i;
  vector e = nullVector();
  for (i = 0; i < dat->N; i++)
  {
    object o = dat->objects[i];
    vector_add_to(&e, vector_diff( vector_cross_prod(o.velocity, vector_cross_prod(o.position, o.velocity))
                                 , scalar_mult(1.0/vector_abs(o.position), o.position) ));
  }
  return vector_abs(e);
}

double semimajor_axis(const data* dat)
{
  // find r_min and r_max
  double r0    = vector_abs(dat->objects[0].position),
         r_min = r0,
         r_max = r0;
  int i;
  for (i = 1; i < dat->N; i++)
  {
    double r = vector_abs(dat->objects[i].position);
    if (r_min > r)
      r_min = r;
    else if (r_max < r)
      r_max = r;
  }
  return 0.5 * (r_min + r_max);
}
