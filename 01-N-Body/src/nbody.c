#include "nbody.h"

// set gravitation constant & total mass to 1
double G = 1;
double M = 1;


// Input file
const char* input = "dat/pl.100";

// Output
const char* output;
output_function current_output = &print_1k;

// timestep modification
double delta_t_factor = 10.0;

// Integrators at work
void integrate(const data* dat)
{
  //set_output("results/1k-10/euler.txt");               euler(dat, current_output);
  //set_output("results/1k-10/euler_cromer.txt");        euler_cromer(dat, current_output);
  //set_output("results/1k-10/runge_kutta.txt");         runge_kutta(dat, current_output);
  //set_output("results/1k-10/leap_frog.txt");           leap_frog(dat, current_output);
  //set_output("results/1k-10/verlet.txt");              verlet(dat, current_output);
  set_output("results/test/hermite.txt");             hermite(dat, current_output);
  set_output("results/test/hermite_iterated_2.txt");  hermite_iterated_2(dat, current_output);
  set_output("results/test/hermite_iterated_3.txt");  hermite_iterated_3(dat, current_output);
}
