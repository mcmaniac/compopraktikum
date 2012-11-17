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

  // output function (%i = current N)
  char *output_p_gp  = "results/a3/poly-%i.gp",
       *output_p_txt = "results/a3/poly-%i.txt",
       *output_l_dat = "results/a3/legendre-%i.dat",
       *output_l_txt = "results/a3/legendre-%i.txt";

  // read data
  matrix dat = read_data("dat/a3.txt");

  int i;
  for (i = 0; i < sizeof(N)/sizeof(int); i++)
  {
    //
    // use polynoms
    //
    matrix F = init_F(dat, N[i], &polynom);
    vector b = init_b(dat, N[i], &polynom);

    // solve system of linear equations
    vector x = solve_gauss(F,b);

    char fp[100];

    // parameter output
    snprintf(fp, sizeof(fp), output_p_txt, N[i]-1);
    printf("Opening file %s...\n", fp);
    FILE *file = fopen(fp, "w+");
    vector_fprint(file, x);
    fclose(file);

    // gnuplot output
    snprintf(fp, sizeof(fp), output_p_gp, N[i]-1);
    printf("Opening file %s...\n", fp);
    file = fopen(fp, "w+");
    fprintf(file, "plot ");
    int k;
    for (k = 0; k < x.N; k++)
    {
      fprintf(file, "%e * x**%i", VectorGET(x,k), k);
      if (k < x.N-1)
        fprintf(file, " + ");
    }
    fprintf(file, "\n");
    fclose(file);

    //
    // Use legendre polynoms
    //
    F = init_F(dat, N[i], &legendre_polynom);
    b = init_b(dat, N[i], &legendre_polynom);

    // solve system of linear equations
    x = solve_gauss(F,b);

    // parameter output
    snprintf(fp, sizeof(fp), output_l_txt, N[i]-1);
    printf("Opening file %s...\n", fp);
    file = fopen(fp, "w+");
    vector_fprint(file, x);

    // "plot" data for legendre polynoms
    snprintf(fp, sizeof(fp), output_l_dat, N[i]-1);
    printf("Opening file %s...\n", fp);
    file = fopen(fp, "w+");
    double r  = -1.0,
           dr = 0.01;
    while (r <= 1.0)
    {
      double sum = 0;
      for (k = 0; k < x.N; k++)
      {
        sum += VectorGET(x,k) * legendre_polynom(k,r);
      }
      fprintf(file, "%f %e\n", r, sum);
      r += dr;
    }
    fclose(file);
    file = NULL;
  }
}

int main()
{

  //a1();
  //a2();
  a3();

  return 0;
}
