#include "a1ac.h"
#include "a1d.h"

void a1(void)
{
  double R = 5.0 * pow(10,-15);

  // max. order of bessel functions
  int l_max = 100;

  // Accuracy for zero point calculations
  double epsilon = 0.0001;

  double x_max = 75,
         x_min = 0;

  // find zero points in 0..75
  find_zero_points(l_max, x_max, x_min, epsilon, R);

  // check orthogonality of j_l(k_il * r) and j_l(k_jl * r)
  check_orthogonality(l_max, R);
}
