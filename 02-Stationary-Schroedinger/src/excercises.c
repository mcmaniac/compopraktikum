#include "a1ac.h"
#include "a1d.h"

#include "a2.h"

#include "a3.h"

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

void a2(void)
{
  int N = 2;
  double epsilon = 0.0001;

  matrix D      = a2_matrix(N),
         P      = unity_matrix(N);

  printf("D:\n");
  matrix_print(D);

  matrix D_diag = jacobi_diag(D, P, epsilon);

  printf("\nD':\n");
  matrix_print(D_diag);

  printf("\nP:\n");
  matrix_print(P);
}

void a3(void)
{
  // Bessel function order
  int l = 4;

  // settings
  double R = 10  * pow(10,-15),
         M = 939 * pow(10,6);

  // plot potential
  int N = 1000;
  double r  = 0,
         dr = R / N;
  printf("Opening results/pot.dat...\n");
  FILE *file = fopen("results/pot.dat", "w+");
  while (r < R)
  {
    fprintf(file, "%e %e\n", r, a3_pot(r));
    r += dr;
  }

  // accuracy for diagonalisation
  double epsilon = 0.0001;

  printf("\nDiagonalisation of H...");
  matrix H = calc_H(l, R, M),
         P = unity_matrix(get_numer_of_zeroes(l)),
         H_diag = jacobi_diag(H, P, epsilon);
  printf("DONE\n");

  printf("Opening results/hamilton.txt...\n");
  file = fopen("results/hamilton.txt", "w+");
  matrix_fprint(file, H);
  fclose(file);

  printf("Opening results/hamilton-diag.txt...\n");
  file = fopen("results/hamilton-diag.txt", "w+");
  matrix_fprint(file, H_diag);
  fclose(file);

  printf("Opening results/trans.txt...\n");
  file = fopen("results/trans.txt", "w+");
  matrix_fprint(file, P);
  fclose(file);

  // Energies are on (i,i) of H_diag
  printf("Opening results/energy.txt...\n");
  file = fopen("results/energy.txt", "w+");
  int i;
  for (i = 0; i < H_diag.N; i++)
    fprintf(file, "%i %e\n", i, MatrixGET(H_diag, i, i) * pow(10,-6));
  fclose(file);

  printf("\na3 DONE\n");
}
