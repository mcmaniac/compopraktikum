#include "nbody.h"

// set gravitation constant & total mass to 1
double G = 1;
double M = 1;

// Output
const char* output;
output_function current_output = &print_1k;

// Input file
const char* input  = "dat/pl.100";

// timestep modification
double delta_t_factor = 1.0;

// Integrators at work
void integrate(const data* dat)
{
  set_output("results/100-1/euler.txt");               euler(dat, current_output);
  set_output("results/100-1/euler_cromer.txt");        euler_cromer(dat, current_output);
  set_output("results/100-1/runge_kutta.txt");         runge_kutta(dat, current_output);
  set_output("results/100-1/leap_frog.txt");           leap_frog(dat, current_output);
  set_output("results/100-1/verlet.txt");              verlet(dat, current_output);
  set_output("results/100-1/hermite.txt");             hermite(dat, current_output);
  set_output("results/100-1/hermite_iterated_2.txt");  hermite_iterated_2(dat, current_output);
  set_output("results/100-1/hermite_iterated_3.txt");  hermite_iterated_3(dat, current_output);
}
