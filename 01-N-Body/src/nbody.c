#include "nbody.h"

// gravitation constant
double G = 1;

// input file
const char* input =
  "in2.txt";
  //"pla3.txt";
  //"pl.100";
  //"pl.1k";

// output file
const char* output =
  "runge_kutta.2";

// integrator to run
void current_integrator(data* dat)
{
  runge_kutta(dat, &print_constants);
}
