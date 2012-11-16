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
  // order of polynoms + 1
  int N[5] = {3, 5, 9, 13, 17};

  // use polynoms
  base_funct f = &polynom;

  // output function (%i = current N)
  char *output = "results/a3/poly-%i.gp";

  // read data
  matrix dat = read_data("dat/a3.txt");

  int i;
  for (i = 0; i < sizeof(N)/sizeof(int); i++)
  {
    matrix F = init_F(dat, N[i], f);
    vector b = init_b(dat, N[i], f);

    // solve system of linear equations
    vector x = solve_gauss(F,b);

    // store results
    char fp[100];
    snprintf(fp, sizeof(fp), output, N[i]-1);
    printf("Opening file %s...\n", fp);
    FILE *out = fopen(fp, "w+");

    // raw data output
    //vector_fprint(out, x);

    // gnuplot output
    fprintf(out, "plot ");
    int k;
    for (k = 0; k < x.N; k++)
    {
      fprintf(out, "%e * x**%i", VectorGET(x,k), k);
      if (k < x.N-1)
        fprintf(out, " + ");
    }
    fprintf(out, "\n");
  }
}

int main()
{

  //a1();
  //a2();
  a3();

  return 0;
}
