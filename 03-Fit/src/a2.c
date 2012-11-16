#include "a2.h"

measurement* read_measurements(const char* fp, int *N)
{
  FILE *file = fopen(fp, "r");
  if (!file)
  {
    printf("Couldn't open file %s.\n", fp);
    return NULL;
  }

  // get number of elements
  int i = 0,
      r = fscanf(file, "%*f %*f %*f\n");
  while (r != EOF)
  {
    i++;
    r = fscanf(file, "%*f %*f %*f\n");
  }
  *N = i;

  // load data
  measurement *dat = (measurement*) malloc((*N)*sizeof(measurement));
  rewind(file);
  for (i = 0; i < *N; i++)
  {
    fscanf(file, "%lf %lf %lf\n", &dat[i].U, &dat[i].I, &dat[i].delta_I);
  }
  return dat;
}

double chi_square_(int N, measurement *m, double a, double b)
{
  double sum = 0;
  int i;
  for (i = 0; i < N; i++)
  {
    // m[i].I = y_i
    // m[i].U = x_i
    sum += pow((m[i].I - a - b*m[i].U)/(m[i].delta_I*2), 2);
  }
  return sum;
}

double S(int N, measurement *m)
{
  double sum = 0;
  int i;
  for (i = 0; i < N; i++)
  {
    sum += 1/pow(m[i].delta_I*2, 2);
  }
  return sum;
}

double S_x(int N, measurement *m)
{
  double sum = 0;
  int i;
  for (i = 0; i < N; i++)
  {
    // x_i = m[i].U
    sum += m[i].U / pow(m[i].delta_I*2, 2);
  }
  return sum;
}

double S_y(int N, measurement *m)
{
  double sum = 0;
  int i;
  for (i = 0; i < N; i++)
  {
    // y_i = m[i].I
    sum += m[i].I / pow(m[i].delta_I*2, 2);
  }
  return sum;
}

double S_xx(int N, measurement *m)
{
  double sum = 0;
  int i;
  for (i = 0; i < N; i++)
  {
    // x_i = m[i].U
    sum += pow(m[i].U, 2) / pow(m[i].delta_I*2, 2);
  }
  return sum;
}

double S_yy(int N, measurement *m)
{
  double sum = 0;
  int i;
  for (i = 0; i < N; i++)
  {
    // y_i = m[i].I
    sum += pow(m[i].I, 2) / pow(m[i].delta_I*2, 2);
  }
  return sum;
}

double S_xy(int N, measurement *m)
{
  double sum = 0;
  int i;
  for (i = 0; i < N; i++)
  {
    // x_i = m[i].U
    // y_i = m[i].I
    sum += m[i].U * m[i].I / pow(m[i].delta_I*2, 2);
  }
  return sum;
}

double D(int N, measurement *m)
{
  return S(N,m) * S_xx(N,m) - pow(S_x(N,m),2);
}

void linear_regression(int N, measurement *m, double *a, double *b, double *sig_a, double *sig_b)
{
  if (a == NULL)
  {
    *b     = S_y(N,m) / S_x(N,m);
    *sig_b = S(N,m) / D(N,m); // passt nich man...
  }
  else
  {
    *a     = ( S_xx(N,m) * S_y(N,m)  - S_x(N,m) * S_xy(N,m) ) / D(N,m);
    *b     = ( S(N,m)    * S_xy(N,m) - S_x(N,m) * S_y(N,m)  ) / D(N,m);
    *sig_a = S_xx(N,m) / D(N,m);
    *sig_b = S(N,m)    / D(N,m);
  }
}
