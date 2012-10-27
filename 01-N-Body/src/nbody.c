#include "nbody.h"

// set gravitation constant & total mass to 1
double G = 1;
double M = 1;

// timestep modification
double delta_t_factor = 0.1;

// global configuration variables
const char* input  = "dat/in2.txt";
const char* output = "results/runge_kutta.txt";

// Current integrator
integrator      current_integrator = &runge_kutta;
output_function current_output     = &print_2body;
