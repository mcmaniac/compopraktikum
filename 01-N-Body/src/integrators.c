#include "integrators.h"

void runge_kutta(data* dat, void (*output)(double time, double delta_t, const data* dat))
{
  int i;
  int N            = dat->N;
  double delta_t   = dat->eta;

  object* objs     = dat->objects;

  vector* rn = (vector*) malloc(N*sizeof(vector));
  vector* vn = (vector*) malloc(N*sizeof(vector));

  vector* a1 = (vector*) malloc(N*sizeof(vector));
  vector* a2 = (vector*) malloc(N*sizeof(vector));
  vector* a3 = (vector*) malloc(N*sizeof(vector));
  vector* a4 = (vector*) malloc(N*sizeof(vector));

  vector* r1 = (vector*) malloc(N*sizeof(vector));
  vector* r2 = (vector*) malloc(N*sizeof(vector));
  vector* r3 = (vector*) malloc(N*sizeof(vector));
  vector* r4 = (vector*) malloc(N*sizeof(vector));

  vector* v1 = (vector*) malloc(N*sizeof(vector));
  vector* v2 = (vector*) malloc(N*sizeof(vector));
  vector* v3 = (vector*) malloc(N*sizeof(vector));
  vector* v4 = (vector*) malloc(N*sizeof(vector));

  double time = 0;
  while (time <= dat->t_max)
  {
    // accelerations for 1st terms
    accelerations(a1, dat);

    for (i = 0; i < N; i++)
    {
      // store "current" values
      rn[i] = dat->objects[i].position;
      vn[i] = dat->objects[i].velocity;

      // 1st terms
      v1[i] = scalar_mult(delta_t, a1[i]);
      r1[i] = scalar_mult(delta_t, v1[i]);

      // update position
      dat->objects[i].position = vector_add(rn[i], scalar_mult(0.5, r1[i]));
    }

    // accelerations for 2nd terms
    accelerations(a2, dat);

    for (i = 0; i < N; i++)
    {
      // 2nd terms
      v2[i] = scalar_mult(delta_t, a2[i]);
      r2[i] = scalar_mult(delta_t, vector_add(vn[i], scalar_mult(0.5, v1[i])));
      // update position
      dat->objects[i].position = vector_add(rn[i], scalar_mult(0.5, r2[i]));
    }

    // accelerations for 3rd terms
    accelerations(a3, dat);

    for (i = 0; i < N; i++)
    {
      // 3rd terms
      v3[i] = scalar_mult(delta_t, a3[i]);
      r3[i] = scalar_mult(delta_t, vector_add(vn[i], scalar_mult(0.5, v2[i])));
      // update position
      dat->objects[i].position = vector_add(rn[i], r3[i]);
    }

    // accelerations for 4th terms
    accelerations(a4, dat);

    for (i = 0; i < N; i++)
    {
      // 4th terms
      v4[i] = scalar_mult(delta_t, a4[i]);
      r4[i] = scalar_mult(delta_t, vector_add(vn[i], v3[i]));
      // update (final) position
      dat->objects[i].position = vector_add(vn[i],
                                 vector_add(scalar_mult(1.0/6.0, v1[i]),
                                 vector_add(scalar_mult(1.0/3.0, v2[i]),
                                 vector_add(scalar_mult(1.0/3.0, v3[i]),
                                            scalar_mult(1.0/6.0, v4[i])))));
      // update (final) velocity
      dat->objects[i].velocity = vector_add(rn[i],
                                 vector_add(scalar_mult(1.0/6.0, r1[i]),
                                 vector_add(scalar_mult(1.0/3.0, r2[i]),
                                 vector_add(scalar_mult(1.0/3.0, r3[i]),
                                            scalar_mult(1.0/6.0, r4[i])))));
    }

    // find new delta_t
    vector *a    = (vector*) malloc(N*sizeof(vector)),
           *adot = (vector*) malloc(N*sizeof(vector));
    accelerations(a, dat);
    adots(adot, dat);
    double min = vector_abs(a[0]) / vector_abs(adot[0]);
    double min_;
    for (i = 1; i < N; i++)
    {
      min_ = vector_abs(a[i]) / vector_abs(adot[i]);
      if (min_ < min)
        min = min_;
    }
    free(a); free(adot);
    delta_t = dat->eta * min;

    // output values for current time
    output(time, delta_t, dat);

    // increase time
    time += delta_t;
  }

  free(rn); free(vn);
  free(r1); free(r2); free(r3); free(r4);
  free(v1); free(v2); free(v3); free(v4);
  free(a1); free(a2); free(a3); free(a4);
}

/*
 * Acceleration according to (1.4)
 */
void accelerations(vector* a, const data* dat)
{
  object* objs = dat->objects;
  int N = dat->N, i, j;
  for (i = 0; i < N; i++)
  {
    a[i] = nullVector();
    for (j = 0; j < N; j++)
    {
      if (j != i)
      {
        vector ri  = objs[i].position,
               rj  = objs[j].position;
        double mj  = objs[j].mass;
        vector rij = vector_diff(rj, ri);
        double r   = vector_abs(rij);
        vector_add_to(&a[i], scalar_mult(G * mj / (r*r*r), rij));
      }
    }
  }
  return;
}

/*
 * Change of acceleration according to (1.5)
 */
void adots(vector* adot, const data* dat)
{
  object* objs = dat->objects;
  int N = dat->N, i, j;
  for (i = 0; i < N; i++)
  {
    adot[i] = nullVector();
    for (j = 0; j < N; j++)
    {
      if (j != i)
      {
        vector ri  = objs[i].position,
               rj  = objs[j].position,
               vi  = objs[i].velocity,
               vj  = objs[j].velocity;
        double mj  = objs[j].mass;
        vector rij = vector_diff(rj, ri),
               vij = vector_diff(vj, vi);
        double r   = vector_abs(rij);
        vector_add_to(&adot[i], vector_diff(
          scalar_mult(1/(r*r*r), vij),
          scalar_mult(3 * vector_mult(vij,rij) / (r*r*r*r*r), rij)
        ));
      }
    }
  }
  return;
}
