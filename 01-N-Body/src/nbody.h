#pragma once
#include "integrators.h"
#include "io.h"

double G,
       M,
       delta_t_factor;

const char *input,
           *output;

integrator      current_integrator;
output_function current_output;
