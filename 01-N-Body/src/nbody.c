#include "nbody.h"

// set gravitation constant & total mass to 1
double G = 1;
double M = 1;

// timestep modification
double delta_t_factor = 0.1;

// global configuration variables
const char* input  = "dat/in2.txt";
const char* output = "results/hermite_iterated_3.txt";

// Current integrator
integrator      current_integrator = &hermite_iterated_3;
output_function current_output     = &print_2body;
