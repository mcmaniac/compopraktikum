#pragma once

#include <stdlib.h>
#include <stdio.h>

#include "typedefs.h"

#include <vector.h>
#include <matrix.h>

matrix read_data(const char* fp);

matrix init_F(const matrix dat, int N, base_funct f);
vector init_b(const matrix dat, int N, base_funct f);

double polynom(int l, double x);
double legendre_polynom(int l, double x);
