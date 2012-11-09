#include <stdlib.h>
#include <math.h>

#include <matrix.h>
#include <vector.h>

#include "a1d.h"
#include "bessel.h"

double potential(double r, double r_0, double a_0, double V_0);
double a3_pot(double r);

matrix calc_H(int l, double R, double M);

matrix calc_T(int l, double M);
matrix calc_V(int l, double R);

//vector calc_eigenvector(matrix P, int i);
double calc_psi(int j, int l, double r, double R, const matrix P);
