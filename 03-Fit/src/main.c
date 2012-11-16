#include "a1.h"
#include "a2.h"
#include "a3.h"

#include <matrix.h>
#include <vector.h>

void printvector(const char *var_name, vector x)
{
  printf("%s = (", var_name);
  int i;
  for (i = 0; i < x.N; i++)
  {
    if (i == x.N-1)
      printf("%e", VectorGET(x,i));
    else
      printf("%e, ", VectorGET(x,i));
  }
  printf(")\n");
}

void a1(void)
{
  int N = 3,
      M = 1;
  data *dat = (data*) malloc(N*sizeof(data));
  dat[0] = (data) { .val = 299793.0, .delta = 2.0 };
  dat[1] = (data) { .val = 299792.0, .delta = 4.5 };
  dat[2] = (data) { .val = 299782.0, .delta = 25.0 };

  double avg   = average_value(N, dat),
         sig_i = sigma_square_intern(N, dat),
         sig_e = sigma_square_extern(M, N, dat);

  printf("c_avg   = %f\n", avg);
  printf("sigma_i = %f\n", sqrt(sig_i));
  printf("sigma_e = %f\n", sqrt(sig_e));
}

void a2(void)
{
  int N;
  measurement *m = read_measurements("dat/a2.txt", &N);
  double a, b, sig_a, sig_b, chi;

  printf("\nLinear regression I = b*U:\n");
  linear_regression(N, m, NULL, &b, NULL, &sig_b);
  chi = chi_square_(N, m, 0, b);
  printf("b = %f +/- %f, chi = %f\n", b, sqrt(sig_b)/2, sqrt(chi));

  printf("\nLinear regression I = a + b*U:\n");
  linear_regression(N, m, &a, &b, &sig_a, &sig_b);
  chi = chi_square_(N, m, a, b);
  printf("a = %f +/- %f, b = %f +/- %f, chi = %f\n", a, sqrt(sig_a)/2, b, sqrt(sig_b)/2, sqrt(chi));
}

void a3(void)
{
  // read data
  int N;
  data *dat = read_data("dat/a3.txt", &N);

  // order of polynoms
  int n = 2;

  matrix F = init_F(N, dat, n, &polynom);
  vector b = init_b(N, dat, n, &polynom);
  printvector("b", b);

  vector x = solve_gauss(F,b);
  printvector("x", x);
}

int main()
{

  //a1();
  //a2();
  a3();

  return 0;
}
