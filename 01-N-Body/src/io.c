#include "nbody.h"

object* read_data(const char* filepath, int* N, double* t_max, double* delta_t)
{
  object* objs = (object*) malloc((*N)*sizeof(object));
  FILE* file = fopen(filepath, "r");
  int i;

  // Read basic data
  if (file)
    fscanf(file, "%i %lf %lf\n", N, t_max, delta_t);
  else
  {
    printf("File \"%s\" not found.\n", filepath);
    free(objs);
    free(file);
    return NULL;
  }

  // Read masses
  double M = 0;
  for (i = 0; i < *N; i++)
  {
    fscanf(file, "%lf\n", &objs[i].mass);
    M += objs[i].mass;
  }
  // Scale masses
  for (i = 0; i < *N; i++)
    objs[i].mass = objs[i].mass / M;

  // Read positions
  for (i = 0; i < *N; i++)
    fscanf(file, "%lf %lf %lf\n", &objs[i].position.x, &objs[i].position.y, &objs[i].position.z);
  // Find center of mass
  vector com = { 0, 0, 0 };
  for (i = 0; i < *N; i++)
    vector_add_to(&com, scalar_mult(objs[i].mass, objs[i].position));
  // Scale positions relative to COM
  for (i = 0; i < *N; i++)
    vector_subtract_from(&objs[i].position, com);

  // Read velocities
  for (i = 0; i < *N; i++)
    fscanf(file, "%lf %lf %lf\n", &objs[i].velocity.x, &objs[i].velocity.y, &objs[i].velocity.z);

  return objs;
}

void print_object(const object o)
{
  vector p = o.position, v = o.velocity;
  printf("M = %lf – P = (%.2lf,%.2lf,%.2lf) – V = (%.2lf,%.2lf,%.2lf)\n", o.mass, p.x,p.y,p.z, v.x,v.y,v.z);
}
