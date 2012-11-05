#include "a1ac.h"

bessel newton_raphson(double epsilon, bessel j)
{
  while (fabs(j.val) > epsilon)
  {
    // get new bessel function for x_(i+1)
    double y = j.val / derive(j).val;
    j = j_l(j.l, j.x - y);
  }
  return j;
}

// Find zero points according to chapter (2.4.1)
// Also stores wave numbers (zero-point / R)
void find_zero_points(int l_max, double x_max, double x_min, double epsilon, double R)
{
  int i, l;

  // Intervallgrenzen
  double width = fabs(x_max - x_min);

  // Anzahl Nullstellen von j_0 ( = sin(x)/x ) ist floor(width / pi):
  int N = floor(width / PI);

  // initialize zero points
  double zeros[N];

  FILE *file_zp = fopen("results/zero-points/0.txt", "w+"),
       *file_wn = fopen("results/wave-numbers/0.txt", "w+");
  for (i = 0; i < N; i++)
  {
    zeros[i] = (i+1) * PI; // i+1 = n
    fprintf(file_zp, "%f %e\n", zeros[i], 0.0);
    fprintf(file_wn, "%f %e\n", zeros[i]/R, 0.0);
  }
  fclose(file_zp);
  fclose(file_wn);

  bessel j;
  char fp_zp[100], fp_wn[100];

  for (l = 1; l <= l_max; l++)
  {
    // open zero-point file
    snprintf(fp_zp, sizeof(fp_zp), "results/zero-points/%i.txt", l);
    printf("Opening file %s...\n", fp_zp);
    file_zp = fopen(fp_zp, "w+");
    if (!file_zp)
      printf("ERROR: Cannot open file %s.", fp_zp);

    // open wave-numbers file
    snprintf(fp_wn, sizeof(fp_wn), "results/wave-numbers/%i.txt", l);
    printf("Opening file %s...\n", fp_wn);
    file_wn = fopen(fp_wn, "w+");
    if (!file_wn)
      printf("ERROR: Cannot open file %s.", fp_wn);

    for (i = 0; i < N-1; i++) // only go up to N-1 !
    {
      // start newton raphson at position x_start = (zeros[i] + zeros[i+1]) / 2
      j = newton_raphson(epsilon, j_l(l, (zeros[i] + zeros[i+1])/2.0));

      if (j.x < zeros[i+1] && j.x <= x_max)
      {
        fprintf(file_zp, "%f %e\n", j.x, j.val);
        fprintf(file_wn, "%f %e\n", j.x/R, j.val);
        zeros[i] = j.x;
      }
      else
      {
        printf("ERROR: Invalid value, j_%i(%f) = %f\n", j.l, j.x, j.val);
      }

    }

    // find last zero point
    if (N > 0)
    {
      // find approx. width of interval between two zero points
      double dN;
      if (N > 2)
        dN = zeros[N-1] - zeros[N-2];
      else
        dN = x_max - zeros[N-1]; // use x_max as highest value instead

      // start at z_N + 1/2 dN
      j = newton_raphson(epsilon, j_l(l, zeros[N-1] + dN/2.0));

      if (zeros[N-1] < j.x && j.x < x_max)
      {
        fprintf(file_zp, "%f %e\n", j.x, j.val);
        fprintf(file_wn, "%f %e\n", j.x/R, j.val);
        zeros[N-1] = j.x;
      }
      else
      {
        N--; // "drop" the last zero point
      }
    }

    fclose(file_zp);
    fclose(file_wn);
  }
}
