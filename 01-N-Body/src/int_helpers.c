#include "main.h"

vector* init_r(const data* dat)
{
  vector *r = (vector*) malloc((dat->N)*sizeof(vector));
  int i;
  for (i = 0; i < dat->N; i++)
    r[i] = dat->objects[i].position;
  return r;
}

vector* init_v(const data* dat)
{
  vector *v = (vector*) malloc((dat->N)*sizeof(vector));
  int i;
  for (i = 0; i < dat->N; i++)
    v[i] = dat->objects[i].velocity;
  return v;
}

double delta_t(const data* dat, const vector* a, const vector* adot)
{
  double min = vector_abs(a[0])/vector_abs(adot[0]),
         min_;
  int i;
  for (i = 1; i < dat->N; i++)
  {
    min_ = vector_abs(a[i])/vector_abs(adot[i]);
    if (min_ < min)
      min = min_;
  }
  return delta_t_factor * dat->eta * min;
}

double delta_t_(const data* dat, const vector* a, const vector* a_1, const vector* a_2, const vector* a_3)
{
  double min = sqrt( (vector_abs(a[0]) * vector_abs(a_2[0]) + vector_abs(a_1[0]) * vector_abs(a_1[0]))
                   / (vector_abs(a_1[0]) * vector_abs(a_3[0]) + vector_abs(a_2[0]) * vector_abs(a_2[0])) ),
         min_;
  int i;
  for (i = 1; i < dat->N; i++)
  {
    min_ = sqrt( (vector_abs(a[i]) * vector_abs(a_2[i]) + vector_abs(a_1[i]) * vector_abs(a_1[i]))
               / (vector_abs(a_1[i]) * vector_abs(a_3[i]) + vector_abs(a_2[i]) * vector_abs(a_2[i])) );
    if (min_ < min)
      min = min_;
  }
  return delta_t_factor * dat->eta * min;
}

//void inc_time_(double* time, double* delta_t, const vector* a, const vector* adot, const vector* a_2, const vector* a_3);

/*
void update_time(double* time, double* delta_t, const data* dat, const vector)
{
  // find new delta_t
  set_new_delta_t(delta_t, dat);
  // increase time
  *time += *delta_t;
}

void update_and_output(double* time, double* delta_t, const data* dat, output_function output)
{
  update_time(time, delta_t, dat);
  // output values for current time
  output(*time, *delta_t, dat);
}

void set_new_delta_t(double* delta_t, const data* dat)
{
  vector *a    = (vector*) malloc((dat->N)*sizeof(vector)),
         *adot = (vector*) malloc((dat->N)*sizeof(vector));
  accelerations(a, dat);
  adots(adot, dat);
  double min = vector_abs(a[0]) / vector_abs(adot[0]);
  int i;
  for (i = 1; i < dat->N; i++)
  {
    double min_ = vector_abs(a[i]) / vector_abs(adot[i]);
    if (min_ < min)
      min = min_;
  }
  *delta_t = dat->eta * min * delta_t_factor;
  free(a); free(adot);
}
*/

/*
 * Acceleration according to (1.4)
 */
void accelerations(const data* dat, const vector* r, vector* a)
{
  int i, j;
  for (i = 0; i < dat->N; i++)
  {
    a[i] = null_vector(3);
    for (j = 0; j < dat->N; j++)
    {
      if (j != i)
      {
        vector ri  = r[i],
               rj  = r[j];
        double mj  = dat->objects[j].mass;
        vector rij = vector_diff(rj, ri);
        double r   = vector_abs(rij);
        vector_add_to(a[i], scalar_mult(G * mj / (r*r*r), rij));
      }
    }
  }
  return;
}

/*
 * Change of acceleration according to (1.5)
 */
void adots(const data* dat, const vector* r, const vector* v, vector* adot)
{
  int i, j;
  for (i = 0; i < dat->N; i++)
  {
    adot[i] = null_vector(3);
    for (j = 0; j < dat->N; j++)
    {
      if (j != i)
      {
        vector ri  = r[i],
               rj  = r[j],
               vi  = v[i],
               vj  = v[j];
        double mj  = dat->objects[j].mass;
        vector rij = vector_diff(rj, ri),
               vij = vector_diff(vj, vi);
        double r   = vector_abs(rij);
        vector_add_to(adot[i], vector_diff(
          scalar_mult(mj / (r*r*r), vij),
          scalar_mult(mj * 3 * vector_mult(vij,rij) / (r*r*r*r*r), rij)
        ));
      }
    }
  }
  return;
}
