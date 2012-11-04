#pragma once
#include "typedefs.h"

/*
 * Output functions
 *
 */
void print_2body(double time, double delta_t, const data* dat, const vector* r, const vector* v);
void print_1k(double time, double delta_t, const data* dat, const vector* r, const vector* v);

void set_output(const char* fp);
