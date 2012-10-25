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

vector* init_r(const data* dat);
vector* init_v(const data* dat);

void inc_time(double* time, double* delta_t);

void accelerations(const data* dat, const vector* r, vector* a);
void adots(const data* dat, const vector* r, const vector* v, vector* adot);
void set_new_delta_t(double* delta_t, const data* dat);
void update_time(double* time, double* delta_t, const data* dat);
void update_and_output(double* time, double* delta_t, const data* dat, output_function output);

/*
 * Conserved quantities
 *
 */

double total_energy(const data* dat, const vector* r, const vector* v);
//double total_momentum(const data* dat);
double total_angular_momentum(const data* dat, const vector* r, const vector* v);
//double total_center_of_mass(const data* dat);
vector runge_lenz(const data* dat, const vector* r, const vector* v, vector* _j);
double total_runge_lenz(const data* dat, const vector* r, const vector* v);
double semimajor_axis(const data* dat, const vector* r, const vector* v);

/*
 * IO - Reading & writing
 *
 */

FILE* file;

data* read_data(const char* file); // , int* N, double* t_max, double* delta_t);

void print_object(const object o);

void print_2body_to_file(const char* filepath, double time, double delta_t, const data* dat, const vector* r, const vector* v);
