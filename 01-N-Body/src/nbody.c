#include "nbody.h"

// set gravitation constant & total mass to 1
double G = 1;
double M = 1;

// global configuration variables
const char* input  = "in2.txt";
const char* output = "2_verlet.txt";

// Current integrator
integrator current_integrator = &verlet;
