#pragma once
#include "vector.h"

/*
 * Typedefs
 *
 */

typedef struct {
  double mass;
  vector position;
  vector velocity;
} object;

typedef struct {
  int N;
  double t_max;
  double eta;
  object* objects;
} data;

typedef void (*output_function)(double time, double delta_t, const data* dat);
typedef void (*integrator)(data* dat, output_function output);
