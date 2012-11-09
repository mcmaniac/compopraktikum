#include "a1ac.h"
#include "a1d.h"

#include "a2.h"

#include "a3.h"
#include "bessel.h"

#include <matrix.h>

void plot_bessel_functions(int l_max, double x_0, double x_R, int N)
{
  int l;
  for (l = 0; l <= l_max; l++)
  {
    // open file
    char fp[100];
    snprintf(fp, sizeof(fp), "results/bessel/%i.txt", l);
    printf("Opening %s...\n", fp);
    FILE *file = fopen(fp, "w+");

    double dx = (x_R - x_0) / N,
           x  = x_0 + dx; // avoid x=0 value
    while (x < x_R)
    {
      bessel j = j_l(l, x);
      fprintf(file, "%f %f\n", x, j.val);
      x += dx;
    }

    fclose(file);
  }
}

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
  //check_orthogonality(l_max, R);
}

void a2(void)
{
  // dimensions
  int Dimensions[4] = {2, 5, 12, 15};

  // diagonalisation accuracy
  double epsilon = 0.0001;

  int i;
  for (i = 0; i < sizeof(Dimensions)/sizeof(int); i++)
  {
    int N = Dimensions[i];
  
    matrix M = a2_matrix(N),
           P = unity_matrix(N),
           D = jacobi_diag(M, P, epsilon);

    char fp[100];
    snprintf(fp, sizeof(fp), "results/diag/%i.txt", N);
    printf("Opening %s...\n", fp);
    FILE *file = fopen(fp, "w+");
    matrix_fprint(file, M);
    fclose(file);

    snprintf(fp, sizeof(fp), "results/diag/%i-diag.txt", N);
    printf("Opening %s...\n", fp);
    file = fopen(fp, "w+");
    matrix_fprint(file, D);
    fclose(file);
  }
}

void a3(void)
{
  char fp[100];

  // settings
  const int N = 1000;
  const double R  = 10.0  * pow(10,-15),
               dr = R / N,
               M  = 939.0 * pow(10,6);

  // plot potential
  printf("\nOpening results/pot_plot.txt...\n");
  FILE *file = fopen("results/pot_plot.txt", "w+");

  double r = 0.0;
  while (r < R)
  {
    fprintf(file, "%e %e\n", r, a3_pot(r));
    r += dr;
  }

  fclose(file);

  // accuracy for diagonalisation
  double epsilon = 0.0001;

  // Bessel function order
  int l, l_max = 3;

  for (l = 0; l <= l_max; l++)
  {
    // H diag
    printf("\nDiagonalisation of H...");
    matrix H = calc_H(l, R, M),
           P = unity_matrix(get_numer_of_zeroes(l)),
           H_diag = jacobi_diag(H, P, epsilon);
    printf("DONE\n");

    // Energies are on (i,i) of H_diag
    snprintf(fp, sizeof(fp), "results/energy/%i.txt", l);
    printf("Opening %s...\n", fp);
    file = fopen(fp, "w+");

    int i;
    for (i = 0; i < H_diag.N; i++)
      fprintf(file, "%i %e\n", i, MatrixGET(H_diag, i, i));

    fclose(file);

    // Plot PSI
    int j, // order of psi
        j_max = get_numer_of_zeroes(l);

    printf("%i\n", j_max);

    for (j = 0; j < j_max; j++)
    {
      double psi_j;

      // open psi_j file
      snprintf(fp, sizeof(fp), "results/psi/%i-%i.txt", l, j);
      printf("Opening %s...\n", fp);
      file = fopen(fp, "w+");

      r = 0.0; // reset r
      while (r < R)
      {
        psi_j = calc_psi(j, l, r, R, P);
        fprintf(file, "%e %e %e\n", r, psi_j, pow(psi_j,2));
        r += dr;
      }
      fclose(file);
    }
  }

  printf("\na3 DONE\n");
}
