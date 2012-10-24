#include "main.h"

FILE* file;

void free_data(data* dat)
{
  free(dat->objects);
  free(dat);
}

data* read_data(const char* filepath)
{
  data* dat = (data*) malloc(sizeof(data));
  object* objs;
  int N;
  double t_max, eta;

  // Read basic data
  FILE* file = fopen(filepath, "r");
  if (file)
    fscanf(file, "%i %lf %lf\n", &N, &t_max, &eta);
  else
  {
    printf("File \"%s\" not found.\n", filepath);
    free(dat);
    free(file);
    return NULL;
  }

  objs = (object*) malloc(N*sizeof(object));

  int i;

  // Read masses
  double _M = 0;
  for (i = 0; i < N; i++)
  {
    fscanf(file, "%lf\n", &objs[i].mass);
    _M += objs[i].mass;
  }
  // Scale masses
  for (i = 0; i < N; i++)
    objs[i].mass = objs[i].mass / _M * M;

  // Read positions
  for (i = 0; i < N; i++)
    fscanf(file, "%lf %lf %lf\n", &objs[i].position.x, &objs[i].position.y, &objs[i].position.z);
  // Find center of mass
  vector com = nullVector();
  for (i = 0; i < N; i++)
    vector_add_to(&com, scalar_mult(objs[i].mass, objs[i].position));
  // Scale positions relative to COM
  for (i = 0; i < N; i++)
    vector_diff_from(&objs[i].position, com);

  // Read velocities
  for (i = 0; i < N; i++)
    fscanf(file, "%lf %lf %lf\n", &objs[i].velocity.x, &objs[i].velocity.y, &objs[i].velocity.z);

  dat->objects = objs;
  dat->N       = N;
  dat->t_max   = t_max;
  dat->eta     = eta;

  return dat;
}

void print_object(const object o)
{
  vector p = o.position, v = o.velocity;
  printf("M = %f - P = (%.2f,%.2f,%.2f) - V = (%.2f,%.2f,%.2f)\n", o.mass, p.x,p.y,p.z, v.x,v.y,v.z);
}

void print_constants_to_file(const char* filepath, double time, double delta_t, const data* dat)
{
  if (file)
  {
    double E   = total_energy(dat);
    double j   = total_angular_momentum(dat);
    double e   = total_runge_lenz(dat);
    double a_e = semimajor_axis(dat);
    fprintf(file, "%f %f %f %f %f %f\n", time, delta_t, E, j, e, a_e);
  }
}

void print_constants(double time, double delta_t, const data* dat)
{
  print_constants_to_file(output, time, delta_t, dat);
}
