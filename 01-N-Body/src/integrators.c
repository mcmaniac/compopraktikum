#include "main.h"

void euler(const data* dat, output_function output)
{
  int i;
  double dt = dat->eta * delta_t_factor;

  vector *a   = (vector*) malloc((dat->N)*sizeof(vector)),
         *a_1 = (vector*) malloc((dat->N)*sizeof(vector)),
         *r   = init_r(dat),
         *v   = init_v(dat);

  double time = 0;
  output(time, dt, dat, r, v);

  while (time < dat->t_max)
  {
    accelerations(dat, r, a);
    adots(dat, r, v, a_1);
    for (i = 0; i < dat->N; i++)
    {
      vector v_neu = vector_add(v[i], scalar_mult(dt, a[i])),
             r_neu = vector_add(r[i], scalar_mult(dt, v[i]));
      r[i] = r_neu;
      v[i] = v_neu;
    }
    dt = delta_t(dat, a, a_1);
    time += dt;
    output(time, dt, dat, r, v);
  }
  free(a); free(r); free(v);
}

void euler_cromer(const data* dat, output_function output)
{
  int i;
  double dt = dat->eta * delta_t_factor;
  vector *a    = (vector*) malloc((dat->N)*sizeof(vector)),
         *a_1 = (vector*) malloc((dat->N)*sizeof(vector)),
         *r    = init_r(dat),
         *v    = init_v(dat);

  double time = 0;
  output(time, dt, dat, r, v);

  while (time < dat->t_max)
  {
    accelerations(dat, r, a);
    adots(dat, r, v, a_1);
    for (i = 0; i < dat->N; i++)
    {
      vector v_neu = vector_add(v[i], scalar_mult(dt, a[i])),
             r_neu = vector_add(r[i], scalar_mult(0.5 * dt, vector_add(v[i], v_neu)));
      r[i] = r_neu;
      v[i] = v_neu;
    }
    dt = delta_t(dat, a, a_1);
    time += dt;
    output(time, dt, dat, r, v);
  }
  free(a); free(r); free(v);
}

void leap_frog(const data* dat, output_function output)
{
  int i;
  double dt = dat->eta * delta_t_factor;

  vector *a    = (vector*) malloc((dat->N)*sizeof(vector)),
         *a_1 = (vector*) malloc((dat->N)*sizeof(vector)),
         *r_p  = (vector*) malloc((dat->N)*sizeof(vector)), // r_(n+1/2), not initialized
         *r    = init_r(dat),                               // r_(n+1)
         *r_n  = (vector*) malloc((dat->N)*sizeof(vector)), // r_(n+3/2), not initialized
         *v    = init_v(dat);                               // v_(n+1)

  double time = 0;
  output(time, dt, dat, r, v);

  // initial step
  for (i = 0; i < dat->N; i++)
    // set r_(1/2)
    r_n[i] = vector_add(r[i], scalar_mult(0.5 * dt, v[i]));

  // regular timesteps
  while (time < dat->t_max)
  {
    // a_(n+1/2)
    accelerations(dat, r_n, a);
    adots(dat, r_n, v, a_1);
    for (i = 0; i < dat->N; i++)
    {
      // store previous values as r_(n+1/2)
      r_p[i] = r_n[i];
      // v_(n+1)
      vector_add_to(&v[i], scalar_mult(dt, a[i]));
      // r_(n+3/2)
      vector_add_to(&r_n[i], scalar_mult(dt, v[i]));
      // build r_(n+1)
      r[i] = scalar_mult(0.5, vector_add(r_p[i], r_n[i]));
    }
    dt = delta_t(dat, a, a_1);
    time += dt;
    output(time, dt, dat, r, v);
  }
  free(a); free(r_p); free(r); free(r_n); free(v);
}

void verlet(const data* dat, output_function output)
{
  int i;
  const double dt = dat->eta * delta_t_factor; // konstant

  vector *a    = (vector*) malloc((dat->N)*sizeof(vector)),
         *a_1 = (vector*) malloc((dat->N)*sizeof(vector)),
         *r_p  = (vector*) malloc((dat->N)*sizeof(vector)), // r_(n-1), not initialized
         *r    = init_r(dat),                               // r_(n)
         *r_n  = (vector*) malloc((dat->N)*sizeof(vector)), // r_(n+1), not initialized
         *v    = init_v(dat);                               // v_(n)

  double time = 0;
  output(time, dt, dat, r, v);

  // set initial data
  accelerations(dat, r, a); // a_0
  adots(dat, r, v, a_1);
  for (i = 0; i < dat->N; i++)
  {
    // set r_(-1)
    r_p[i] = vector_add(r[i],
             vector_add(scalar_mult(dt, v[i]),
                        scalar_mult(0.5 * dt * dt, a[i])));
    // set r_1
    r_n[i] = vector_add(scalar_mult(2, r[i]),
             vector_add(scalar_mult(-1, r_p[i]),
                        scalar_mult(dt * dt, a[i])));
  }
  //time += dt;

  // regular timesteps
  while (time < dat->t_max)
  {
    accelerations(dat, r_n, a); // a_n+1 (gets shifted to a_n)
    for (i = 0; i < dat->N; i++)
    {
      // shift indexes n+1 -> n
      r_p[i] = r[i];
      r[i]   = r_n[i];
      r_n[i] = vector_add(scalar_mult(2, r[i]),
               vector_add(scalar_mult(-1, r_p[i]),
                          scalar_mult(dt * dt, a[i])));
      v[i] = scalar_mult(0.5 / dt, vector_diff(r_n[i], r_p[i]));
    }
    adots(dat, r, v, a_1);
    time += dt;
    output(time, dt, dat, r, v);
  }
  free(a); free(r_p); free(r); free(r_n); free(v);
}

void hermite(const data* dat, output_function output)
{
  int i;
  double dt = dat->eta * delta_t_factor;

  vector *a      = (vector*) malloc((dat->N)*sizeof(vector)), // a_n,          not initialized
         *a_p    = (vector*) malloc((dat->N)*sizeof(vector)), // a_(n+1)^p,    not initialized
         *a_1   = (vector*) malloc((dat->N)*sizeof(vector)), // a_1_n,       not initialized
         *a_1_p = (vector*) malloc((dat->N)*sizeof(vector)), // a_1_(n+1)^p, not initialized
         *a_2    = (vector*) malloc((dat->N)*sizeof(vector)), // a_(n+1)^(2),  not initialized
         *a_3    = (vector*) malloc((dat->N)*sizeof(vector)), // a_(n+1)^(3),  not initialized
         *r      = init_r(dat),                               // r_n
         *r_p    = (vector*) malloc((dat->N)*sizeof(vector)), // r_(n+1)^p,    not initialized
         *v      = init_v(dat),                               // v_n
         *v_p    = (vector*) malloc((dat->N)*sizeof(vector)); // v_(n+1)^p,    not initialized

  double time = 0;

  accelerations(dat, r, a);
  adots(dat, r, v, a_1);
  dt = delta_t(dat, a, a_1);

  output(time, dt, dat, r, v);

  while (time < dat->t_max)
  {
    // prediction step
    for (i = 0; i < dat->N; i++)
    {
      v_p[i] = vector_add(v[i],
               vector_add(scalar_mult(dt, a[i]),
                          scalar_mult(0.5 * dt * dt, a_1[i])));
      r_p[i] = vector_add(r[i],
               vector_add(scalar_mult(dt, v[i]),
               vector_add(scalar_mult(0.5 * dt * dt, a[i]),
                          scalar_mult(dt * dt * dt / 6.0, a_1[i]))));
    }
    // predict a^p, a_1^p
    accelerations(dat, r_p, a_p);
    adots(dat, r_p, v_p, a_1_p);
    for (i = 0; i < dat->N; i++)
    {
      // predict 2nd and 3rd derivations
      a_2[i] = vector_add(scalar_mult(2 * (-3) / dt / dt, vector_diff(a[i], a_p[i])),
                          scalar_mult(2 * (-1) / dt, vector_add(scalar_mult(2, a_1[i]), a_1_p[i])));
      a_3[i] = vector_add(scalar_mult(6 * 2 / dt / dt / dt, vector_diff(a[i], a_p[i])),
                          scalar_mult(6 / dt / dt, vector_add(a_1[i], a_1_p[i])));
      // correction steps
      v[i] = vector_add(v_p[i],
             vector_add(scalar_mult(dt * dt * dt / 6.0, a_2[i]),
                        scalar_mult(dt * dt * dt / 24.0, a_3[i])));
      r[i] = vector_add(r_p[i],
             vector_add(scalar_mult(dt * dt * dt * dt / 24.0, a_2[i]),
                        scalar_mult(dt * dt * dt * dt * dt / 120.0, a_3[i])));
    }
    time += dt;
    // a/a_1 values for next step
    accelerations(dat, r, a);
    adots(dat, r, v, a_1);
    // new dt
    dt = delta_t_(dat, a, a_1, a_2, a_3);
    output(time, dt, dat, r, v);
  }

  free(a);    free(a_p);
  free(a_1); free(a_1_p);
  free(r);    free(r_p);
  free(v);    free(v_p);
}

void hermite_iterated_2(const data* dat, output_function output)
{
  hermite_iterated(dat, output, 2);
}

void hermite_iterated_3(const data* dat, output_function output)
{
  hermite_iterated(dat, output, 3);
}

void hermite_iterated(const data* dat, output_function output, int iterations)
{
  int i,j;
  double dt = dat->eta * delta_t_factor;

  vector *a      = (vector*) malloc((dat->N)*sizeof(vector)), // a_n,          not initialized
         *a_p    = (vector*) malloc((dat->N)*sizeof(vector)), // a_(n+1)^p,    not initialized
         *a_1   = (vector*) malloc((dat->N)*sizeof(vector)), // a_1_n,       not initialized
         *a_1_p = (vector*) malloc((dat->N)*sizeof(vector)), // a_1_(n+1)^p, not initialized
         *r      = init_r(dat),                               // r_n
         *r_p    = (vector*) malloc((dat->N)*sizeof(vector)), // r_(n+1)^p,    not initialized
         *v      = init_v(dat),                               // v_n
         *v_p    = (vector*) malloc((dat->N)*sizeof(vector)); // v_(n+1)^p,    not initialized

  double time = 0;

  accelerations(dat, r, a);
  adots(dat, r, v, a_1);
  dt = delta_t(dat, a, a_1);

  output(time, dt, dat, r, v);

  while (time < dat->t_max)
  {
    // prediction step
    for (i = 0; i < dat->N; i++)
    {
      v_p[i] = vector_add(v[i],
               vector_add(scalar_mult(dt, a[i]),
                          scalar_mult(0.5 * dt * dt, a_1[i])));
      r_p[i] = vector_add(r[i],
               vector_add(scalar_mult(dt, v[i]),
               vector_add(scalar_mult(0.5 * dt * dt, a[i]),
                          scalar_mult(dt * dt * dt / 6.0, a_1[i]))));
    }
    // iteration of correction step
    for (j = 0; j < iterations; j++)
    {
      // predict a, a_1
      accelerations(dat, r_p, a_p);
      adots(dat, r_p, v_p, a_1_p);
      for (i = 0; i < dat->N; i++)
      {
        // correction steps -> overwrite "prediction" variables for convenience,
        // will get written in r,v after iteration is done
        v_p[i] = vector_add(v[i],
               vector_add(scalar_mult(dt / 2.0, vector_add(a_p[i], a[i])),
                          scalar_mult(dt * dt / 12.0, vector_diff(a_1_p[i], a_1[i]))));
        r_p[i] = vector_add(r[i],
                 vector_add(scalar_mult(dt / 2.0, vector_add(v_p[i], v[i])),
                            scalar_mult(dt * dt / 12.0, vector_diff(a_p[i], a[i]))));
      }
    }
    // apply iterated correction steps
    for (i = 0; i < dat->N; i++)
    {
      r[i] = r_p[i];
      v[i] = v_p[i];
    }
    time += dt;
    // a/a_1 values for next step
    accelerations(dat, r, a);
    adots(dat, r, v, a_1);
    // new dt
    dt = delta_t(dat, a, a_1);
    output(time, dt, dat, r, v);
  }

  free(a);   free(a_p);
  free(a_1); free(a_1_p);
  free(r);   free(r_p);
  free(v);   free(v_p);
}

void runge_kutta(const data* dat, output_function output)
{
  int i;
  double dt = dat->eta * delta_t_factor;

  vector *r   = init_r(dat),
         *v   = init_v(dat),
         *a   = (vector*) malloc((dat->N)*sizeof(vector)),
         *a_1 = (vector*) malloc((dat->N)*sizeof(vector));

  // temporary r vector for acceleration calculations
  vector *rt = (vector*) malloc((dat->N)*sizeof(vector));

  vector *a1 = (vector*) malloc((dat->N)*sizeof(vector)),
         *a2 = (vector*) malloc((dat->N)*sizeof(vector)),
         *a3 = (vector*) malloc((dat->N)*sizeof(vector)),
         *a4 = (vector*) malloc((dat->N)*sizeof(vector));

  vector *r1 = (vector*) malloc((dat->N)*sizeof(vector)),
         *r2 = (vector*) malloc((dat->N)*sizeof(vector)),
         *r3 = (vector*) malloc((dat->N)*sizeof(vector)),
         *r4 = (vector*) malloc((dat->N)*sizeof(vector));

  vector *v1 = (vector*) malloc((dat->N)*sizeof(vector)),
         *v2 = (vector*) malloc((dat->N)*sizeof(vector)),
         *v3 = (vector*) malloc((dat->N)*sizeof(vector)),
         *v4 = (vector*) malloc((dat->N)*sizeof(vector));

  double time = 0;
  while (time < dat->t_max)
  {
    accelerations(dat, r, a1);
    for (i = 0; i < dat->N; i++)
    {
      v1[i] = scalar_mult(dt, a1[i]);
      r1[i] = scalar_mult(dt, v[i]);
      // temp r for a2
      rt[i] = vector_add(r[i], scalar_mult(0.5, r1[i]));
    }
    accelerations(dat, rt, a2);
    for (i = 0; i < dat->N; i++)
    {
      v2[i] = scalar_mult(dt, a2[i]);
      r2[i] = scalar_mult(dt, vector_add(v[i], scalar_mult(0.5, v1[i])));
      // temp r for a3
      rt[i] = vector_add(r[i], scalar_mult(0.5, r2[i]));
    }
    accelerations(dat, rt, a3);
    for (i = 0; i < dat->N; i++)
    {
      v3[i] = scalar_mult(dt, a3[i]);
      r3[i] = scalar_mult(dt, vector_add(v[i], scalar_mult(0.5, v2[i])));
      // temp r for a3
      rt[i] = vector_add(r[i], scalar_mult(0.5, r3[i]));
    }
    accelerations(dat, rt, a4);
    for (i = 0; i < dat->N; i++)
    {
      v4[i] = scalar_mult(dt, a4[i]);
      r4[i] = scalar_mult(dt, vector_add(v[i], v3[i]));
      // calculate v_(n+1) and r_(n+1)
      vector_add_to(&v[i], vector_add(scalar_mult(1.0/6.0, v1[i]),
                           vector_add(scalar_mult(1.0/3.0, v2[i]),
                           vector_add(scalar_mult(1.0/3.0, v3[i]),
                                      scalar_mult(1.0/6.0, v4[i])))));
      vector_add_to(&r[i], vector_add(scalar_mult(1.0/6.0, r1[i]),
                           vector_add(scalar_mult(1.0/3.0, r2[i]),
                           vector_add(scalar_mult(1.0/3.0, r3[i]),
                                      scalar_mult(1.0/6.0, r4[i])))));
    }
    // increase time
    adots(dat, r, v, a_1);
    dt = delta_t(dat, a1, a_1); // a1 = a(t_n), a_1 = a'(t_n)
    time += dt;
    output(time, dt, dat, r, v);
  }
  free(r);  free(v);  free(a);  free(a_1);
  free(rt);
  free(r1); free(r2); free(r3); free(r4);
  free(v1); free(v2); free(v3); free(v4);
  free(a1); free(a2); free(a3); free(a4);
}
