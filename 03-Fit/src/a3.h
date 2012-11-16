#pragma once

#include <stdlib.h>
#include <stdio.h>

#include "typedefs.h"

#include <vector.h>
#include <matrix.h>

data* read_data(const char* fp, int *N);

matrix init_F(int N, data *dat, int n, double (*f)(int, double));
vector init_b(int N, data *dat, int n, double (*f)(int, double));

double polynom(int l, double x);
double legendre_polynom(int l, double x);
