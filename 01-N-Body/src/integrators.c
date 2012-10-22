#include "integrators.h"

void runge_kutta(data* dat, void (*print_data)(data* dat))
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

  double time  = 0;
  while (time <= dat->t_max)
  {

    // TODO: delta_t bestimmen?

    // store "current" values
    for (i = 0; i < N; i++)
    {
      rn[i] = dat->objects[i].position;
      vn[i] = dat->objects[i].velocity;
    }

    // accelerations for 1st terms
    a1 = accelerations(dat);
    // 1st terms
    for (i = 0; i < N; i++)
    {
      v1[i] = scalar_mult(delta_t, a1[i]);
      r1[i] = scalar_mult(delta_t, v1[i]);
    }

    // accelerations for 2nd terms
    for (i = 0; i < N; i++)
      dat->objects[i].position = vector_add(rn[i], scalar_mult(0.5, r1[i]));
    a2 = accelerations(dat);
    // 2nd terms
    for (i = 0; i < N; i++)
    {
      v2[i] = scalar_mult(delta_t, a2[i]);
      r2[i] = scalar_mult(delta_t, vector_add(vn[i], scalar_mult(0.5, v1[i])));
    }

    // accelerations for 3rd terms
    for (i = 0; i < N; i++)
      dat->objects[i].position = vector_add(rn[i], scalar_mult(0.5, r2[i]));
    a3 = accelerations(dat);
    // 3rd terms
    for (i = 0; i < N; i++)
    {
      v3[i] = scalar_mult(delta_t, a3[i]);
      r3[i] = scalar_mult(delta_t, vector_add(vn[i], scalar_mult(0.5, v2[i])));
    }

    // accelerations for 4th terms
    for (i = 0; i < N; i++)
      dat->objects[i].position = vector_add(rn[i], r3[i]);
    a4 = accelerations(dat);
    // 4th terms
    for (i = 0; i < N; i++)
    {
      v4[i] = scalar_mult(delta_t, a4[i]);
      r4[i] = scalar_mult(delta_t, vector_add(vn[i], v3[i]));
    }

    // set new positions & velocities
    for (i = 0; i < N; i++)
    {
      dat->objects[i].position = vector_add(vn[i],
                                 vector_add(scalar_mult(1.0/6.0, v1[i]),
                                 vector_add(scalar_mult(1.0/3.0, v2[i]),
                                 vector_add(scalar_mult(1.0/3.0, v3[i]),
                                            scalar_mult(1.0/6.0, v4[i])))));
      dat->objects[i].velocity = vector_add(rn[i],
                                 vector_add(scalar_mult(1.0/6.0, r1[i]),
                                 vector_add(scalar_mult(1.0/3.0, r2[i]),
                                 vector_add(scalar_mult(1.0/3.0, r3[i]),
                                            scalar_mult(1.0/6.0, r4[i])))));
    }

    // print values
    print_data(dat);
  }

  free(rn);
  free(vn);
  free(a1);
  free(a2);
  free(a3);
  free(a4);
  free(v1);
  free(v2);
  free(v3);
  free(v4);
  free(r1);
  free(r2);
  free(r3);
  free(r4);
}

/*
 * Acceleration according to (1.4)
 */
vector* accelerations(const data* dat)
{
  object* objs = dat->objects;
  int N = dat->N, i, j;
  vector* a = (vector*) malloc(N*sizeof(vector));
  for (i = 0; i < N; i++)
  {
    a[i] = (vector) {0,0,0};
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
  return a;
}

/*
 * Change of acceleration according to (1.5)
 */
vector* adots(const data* dat)
{
  object* objs = dat->objects;
  int N = dat->N, i, j;
  vector* adot = (vector*) malloc(N*sizeof(vector));
  for (i = 0; i < N; i++)
  {
    adot[i] = (vector) {0,0,0};
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
  return adot;
}
