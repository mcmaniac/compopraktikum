#include "nbody.h"

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
  vector L = nullVector();
  int i;
  for (i = 0; i < dat->N; i++)
  {
    object o = dat->objects[i];
    vector_add_to(&L, vector_cross_prod(o.position, o.velocity));
  }
  return vector_abs(L);
}

double total_center_of_mass(const data* dat)
{
  int i;
  vector com = nullVector();
  for (i = 0; i < dat->N; i++)
    vector_add_to(&com, dat->objects[i].position);
  return vector_abs(com);
}
