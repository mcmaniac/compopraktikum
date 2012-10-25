#pragma once
#include "typedefs.h"

void euler(const data* dat, output_function);
void euler_cromer(const data* dat, output_function output);
void leap_frog(const data* dat, output_function output);
void verlet(const data* dat, output_function output);
void hermite(const data* dat, output_function output);
void hermite_iterated_2(const data* dat, output_function output);
void hermite_iterated_3(const data* dat, output_function output);
void hermite_iterated(const data* dat, output_function output, int iterations);
//void runge_kutta(data* dat, output_function);
