#pragma once
#include "typedefs.h"

void euler(data* dat, output_function);
void euler_cromer(data* dat, output_function);
void leap_frog(data* dat, output_function);
void verlet(data* dat, output_function);
void runge_kutta(data* dat, output_function);
