#include "nbody.h"

// gravitation constant
double G = 1;

// global configuration variables
const char* input  = "in2.txt";
const char* output = "2_runge_kutta.txt";

// Current integrator
integrator current_integrator = &runge_kutta;
