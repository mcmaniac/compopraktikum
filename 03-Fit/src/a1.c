#include "a1.h"

double average_value(int N, data *dat)
{
  double sum1 = 0,
         sum2 = 0;
  int i;
  for (i = 0; i < N; i++)
  {
    sum1 += dat[i].val / pow(dat[i].err, 2);
    sum2 += 1          / pow(dat[i].err, 2);
  }
  return sum1 / sum2;
}

double chi_square(int N, data *dat)
{
  double sum   = 0,
         c_avg = average_value(N, dat);
  int i;
  for (i = 0; i < N; i++)
  {
    sum += pow( (dat[i].val - c_avg) / dat[i].err , 2);
  }
  return sum;
}

double sigma(data dat)
{
  return dat.err * 2;
}

double sigma_square_intern(int N, data *dat)
{
  double sum = 0;
  int i;
  for (i = 0; i < N; i++)
  {
    sum += 1.0 / pow(sigma(dat[i]), 2);
  }
  return 1.0/sum;
}

double sigma_square_extern(int M, int N, data *dat)
{
  double chi_2 = chi_square(N, dat),
         sig_2 = sigma_square_intern(N, dat);
  return chi_2 / (N - M) * sig_2;
}
