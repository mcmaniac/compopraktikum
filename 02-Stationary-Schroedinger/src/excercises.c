#include "a1ac.h"
#include "a1d.h"

#include "a2.h"

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
  int N = 15;
  double epsilon = 0.0001;

  matrix D      = a2_matrix(N),
         P      = unity_matrix(N),
         Dprime = jacobi_diag(D, P, epsilon);

  printf("D:\n");
  matrix_print(D);

  printf("\nD':\n");
  matrix_print(Dprime);

  printf("\nP:\n");
  matrix_print(P);

  /* Matrix mult test
  matrix A = null_matrix(2, 2),
         B = null_matrix(2, 2);
  MatrixSET(A, 1, 2, 5);
  MatrixSET(A, 2, 2, 3);
  MatrixSET(B, 1, 1, 2);
  MatrixSET(B, 2, 1, 2);

  matrix C = matrix_mult(A,B);
  printf("C_1,1 = %f\n", MatrixGET(C, 1, 1));
  printf("C_1,2 = %f\n", MatrixGET(C, 1, 2));
  printf("C_2,1 = %f\n", MatrixGET(C, 2, 1));
  printf("C_2,2 = %f\n", MatrixGET(C, 2, 2));
  */
}
