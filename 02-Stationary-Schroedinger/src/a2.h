#pragma once

#include <stdlib.h>
#include <math.h>

#include "matrix.h"

matrix a2_matrix(int N);

//void calc_Dprime(matrix D, matrix P, int p, int q);
matrix jacobi_diag(const matrix A, matrix P, double epsilon);
