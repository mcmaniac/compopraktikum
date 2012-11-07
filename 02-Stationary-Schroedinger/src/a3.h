#include <stdlib.h>
#include <math.h>

#include "a1d.h"
#include "matrix.h"

double potential(double r, double r_0, double a_0, double V_0);
double a3_pot(double r);

matrix calc_H(int l, double R, double M);
