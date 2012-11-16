#include "a3.h"

data* read_data(const char* fp, int *N)
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
  data *dat = (data*) malloc((*N)*sizeof(data));
  rewind(file);
  for (i = 0; i < *N; i++)
  {
    fscanf(file, "%lf %lf\n", &dat[i].val, &dat[i].delta);
  }
  return dat;
}

matrix init_F(int N, data *dat, int n, double (*f)(int, double))
{
  matrix F = null_matrix(n, n);
  int i, j;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      double sum = 0;
      int k;
      for (k = 0; k < N; k++)
      {
        sum += 1 / pow(dat[k].delta*2, 2) * f(i, dat[k].val) * f(j, dat[k].val);
      }
      MatrixSET(F, i, j, sum);
    }
  }
  return F;
}

vector init_b(int N, data *dat, int n, double (*f)(int, double))
{
  vector b = null_vector(n);
  int i;
  for (i = 0; i < n; i++)
  {
    double sum = 0;
    int k;
    for (k = 0; k < N; k++)
    {
      sum += 1 / pow(dat[k].delta*2, 2) * f(i, dat[k].val);
    }
    VectorSET(b, i, sum);
  }
  return b;
}

double polynom(int l, double x)
{
  return pow(x,l);
}

double legendre_polynom(int l, double x)
{
  if (l == 0)
    return x;
  else if (l == 1)
    return x;
  else
    // rekursion
    return (2*l-1)/l * x * legendre_polynom(l-1,x) - (l-1)/l * legendre_polynom(l-2,x);
}
