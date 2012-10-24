#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nbody.h"

int main();

/*
 * Typedef related functions
 *
 */

void free_data(data* dat);

/*
 * Integrator helper functions
 *
 */

void accelerations(vector* a, const data* dat);
void adots(vector* adot, const data* dat);

/*
 * Conserved quantities
 *
 */

double total_energy(const data* dat);
double total_momentum(const data* dat);
double total_angular_momentum(const data* dat);
double total_center_of_mass(const data* dat);
double total_runge_lenz(const data* dat);
double semimajor_axis(const data* dat);

/*
 * IO - Reading & writing
 *
 */

FILE* file;

data* read_data(const char* file); // , int* N, double* t_max, double* delta_t);

void print_object(const object o);

void print_constants(double time, double delta_t, const data* dat);
void print_constants_to_file(const char* filepath, double time, double delta_t, const data* dat);
