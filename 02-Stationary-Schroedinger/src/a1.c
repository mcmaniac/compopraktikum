#include "a1ac.h"
#include "a1d.h"

int a1()
{
  double R = 5.0 * pow(10,-15);

  // max. order of bessel functions
  int l_max = 125;

  // Accuracy for zero point calculations
  //double epsilon = 0.0001;

  //double x_max = 75,
  //       x_min = 0;
  //find_zero_points(l_max, x_max, x_min, epsilon, R);

  // change order of bessel functions back to 100
  l_max = 100;

  check_orthogonality(l_max, R);

  return 0;
}
