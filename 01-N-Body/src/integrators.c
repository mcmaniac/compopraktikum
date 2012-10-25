#include "main.h"

void euler(data* dat, output_function output)
{
  int i;
  double delta_t = dat->eta;
  vector  r, v;
  vector* a = (vector*) malloc((dat->N)*sizeof(vector));

  double time = 0;
  while (time <= dat->t_max)
  {
    accelerations(a, dat);
    for (i = 0; i < dat->N; i++)
    {
      object* o = &dat->objects[i];
      v = vector_add(o->velocity, scalar_mult(delta_t, a[i]));
      r = vector_add(o->position, scalar_mult(delta_t, o->velocity));
      o->position = r;
      o->velocity = v;
    }
    update_and_output(&time, &delta_t, dat, output);
  }
  free(a);
}

void euler_cromer(data* dat, output_function output)
{
  int i;
  double delta_t = dat->eta;
  vector r, v;
  vector* a = (vector*) malloc((dat->N)*sizeof(vector));

  double time = 0;
  while (time <= dat->t_max)
  {
    accelerations(a, dat);
    for (i = 0; i < dat->N; i++)
    {
      object* o = &dat->objects[i];
      v = vector_add(o->velocity, scalar_mult(delta_t, a[i]));
      r = vector_add(o->position, scalar_mult(0.5 * delta_t, vector_add(o->velocity, v)));
      o->position = r;
      o->velocity = v;
    }
    update_and_output(&time, &delta_t, dat, output);
  }
  free(a);
}

void leap_frog(data* dat, output_function output)
{
  int i;
  vector* a = (vector*) malloc((dat->N)*sizeof(vector));
  double delta_t = dat->eta;
  double time = 0;

  // first step
  for (i = 0; i < dat->N; i++)
  {
    vector_add_to(&dat->objects[i].position, scalar_mult(0.5 * delta_t, dat->objects[i].velocity));
    // manually increase time by 1/2*delta_t
    set_new_delta_t(&delta_t, dat);
    time += delta_t * 0.5;
    output(time, delta_t, dat);
  }

  // regular time steps
  while (time < (dat->t_max - 0.5*delta_t))
  {
    accelerations(a, dat);
    for (i = 0; i < dat->N; i++)
    {
      object* o = &dat->objects[i];
      vector_add_to(&o->velocity, scalar_mult(delta_t, a[i]));
      vector_add_to(&o->position, scalar_mult(delta_t, o->velocity));
    }
    update_and_output(&time, &delta_t, dat, output);
  }

  // final time step
  for (i = 0; i < dat->N; i++)
    vector_add_to(&dat->objects[i].position, scalar_mult(0.5 * delta_t, dat->objects[i].velocity));
  // manually increase time by 1/2*delta_t
  set_new_delta_t(&delta_t, dat);
  time += delta_t * 0.5;
  output(time, delta_t, dat);

  free(a);
}

// funktioniert nicht :)
void verlet(data* dat, output_function output)
{
  int i;
  double delta_t = dat->eta;

  vector* a      = (vector*) malloc((dat->N)*sizeof(vector));
  vector* r_prev = (vector*) malloc((dat->N)*sizeof(vector));
  vector* r_cur  = (vector*) malloc((dat->N)*sizeof(vector));

  // set r_(-1)
  accelerations(a, dat);
  for (i = 0; i < dat->N; i++)
  {
    object o = dat->objects[i];
    r_prev[i] = vector_add(o.position,
                vector_add(scalar_mult(-delta_t, o.velocity),
                           scalar_mult(delta_t*delta_t/2.0, a[i])));
  }

  double time = 0;
  while(time < dat->t_max)
  {
    for (i = 0; i < dat->N; i++)
    {
      object* o = &dat->objects[i];
      r_cur[i]  = o->position;
      vector_add_to(&o->position, vector_add(scalar_mult(2, o->position),
                                  vector_add(scalar_mult(-1, r_prev[i]),
                                             scalar_mult(delta_t*delta_t, a[i]))));
      o->velocity = scalar_mult(0.5/delta_t, vector_diff(o->position, r_prev[i]));
      r_prev[i] = r_cur[i];
    }
    update_and_output(&time, &delta_t, dat, output);
  }

  // last time step
  for (i = 0; i < dat->N; i++)
    dat->objects[i].velocity = scalar_mult(0.5/delta_t, vector_diff(dat->objects[i].position, r_prev[i]));
  update_and_output(&time, &delta_t, dat, output);

  free(a); free(r_prev); free(r_cur);
}

void runge_kutta(data* dat, output_function output)
{
  int i;
  int N            = dat->N;
  double delta_t   = dat->eta;

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
      dat->objects[i].position = vector_add(rn[i],
                                 vector_add(scalar_mult(1.0/6.0, r1[i]),
                                 vector_add(scalar_mult(1.0/3.0, r2[i]),
                                 vector_add(scalar_mult(1.0/3.0, r3[i]),
                                            scalar_mult(1.0/6.0, r4[i])))));
      // update (final) velocity
      dat->objects[i].velocity = vector_add(vn[i],
                                 vector_add(scalar_mult(1.0/6.0, v1[i]),
                                 vector_add(scalar_mult(1.0/3.0, v2[i]),
                                 vector_add(scalar_mult(1.0/3.0, v3[i]),
                                            scalar_mult(1.0/6.0, v4[i])))));
    }
    update_and_output(&time, &delta_t, dat, output);
  }

  free(rn); free(vn);
  free(r1); free(r2); free(r3); free(r4);
  free(v1); free(v2); free(v3); free(v4);
  free(a1); free(a2); free(a3); free(a4);
}

/*
 * Update & output current values
 */
void update_and_output(double* time, double* delta_t, const data* dat, output_function output)
{
  // find new delta_t
  set_new_delta_t(delta_t, dat);

  // increase time
  *time += *delta_t;

  // output values for current time
  output(*time, *delta_t, dat);
}

/*
 * Find new delta_t
 */
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
  *delta_t = dat->eta * min;
  free(a); free(adot);
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
          scalar_mult(mj / (r*r*r), vij),
          scalar_mult(mj * 3 * vector_mult(vij,rij) / (r*r*r*r*r), rij)
        ));
      }
    }
  }
  return;
}
