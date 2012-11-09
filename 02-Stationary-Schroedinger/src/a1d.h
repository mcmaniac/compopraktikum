#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bessel.h"
#include "typedefs.h"

void check_orthogonality(int l_max, double R);

// for a3
int get_numer_of_zeroes(int l);
double get_kjl(int j, int l, double R);
double norm_factor(double kjl, int j, int l, double R);
