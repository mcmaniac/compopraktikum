#include "a1d.h"

// Wave number file pointer
FILE *file_wn;

void open_wave_number_file(int l)
{
  // close open file first
  if (file_wn)
    fclose(file_wn);

  // open new file
  char fp[100];
  snprintf(fp, sizeof(fp), "results/wave-numbers/%i.txt", l);
  file_wn = fopen(fp, "r");
}

void close_wave_number_file()
{
  if (file_wn)
    fclose(file_wn);
}

// get number of zeroes from file
int get_numer_of_zeroes(int l)
{
  int i = 0, r = 0;
  double z;

  char fp[100];
  snprintf(fp, sizeof(fp), "results/zero-points/%i.txt", l);
  FILE *file = fopen(fp, "r");
  if (!file)
    printf("ERROR: Cannot open file %s.", fp);

  r = fscanf(file, "%lf %*e\n", &z);
  while (r != EOF)
  {
    i++;
    r = fscanf(file, "%lf %*e\n", &z);
  }

  fclose(file);

  return i;
}

// Load wave number from file to avoid re-calculation
double get_kjl(int j, int l)
{
  int i = 0, r = 0;
  double kjl;

  // reset position in file
  rewind(file_wn);

  // look for wave-number
  while (r != EOF && i <= j)
  {
    r = fscanf(file_wn, "%lf %*e\n", &kjl);
    i++;
  }

  // last stored value is what we're looking for
  return kjl;
}

// Factor to normalize j_l(k_jl * r) with r < R
double norm_factor(double kjl, int j, int l, double R)
{
  if (l == 0)
    return sqrt( 2.0 / pow(R,3) ) * (j+1) * PI;
  else
    return sqrt( 2.0 / pow(R,3) ) / j_l(l-1, kjl * R).val;
}

// integration from 0..R (R = 5fm) according to eq. (2.8) and (2.21)
double integrate(int i, int j, int l, double R)
{
  // wave numbers & norm factors
  double kil = get_kjl(i, l),
         kjl = get_kjl(j, l),
         ail = norm_factor(kil, i, l, R),
         ajl = norm_factor(kjl, j, l, R);

  // number of timesteps
  int    N  = 100;
  double dr = R / N,
         r  = 0;

  // r = 0 value
  double res = 0;

  // 0 < r < R values
  while (r+dr < R)
  {
    r   += dr;
    res += pow(r,2) * j_l(l, kil * r).val * j_l(l, kjl * r).val;
  }

  // r = R value
  res += 0.5 * pow(r,2) * j_l(l, kil * R).val * j_l(l, kjl * R).val;

  return res * dr * ail * ajl;
}

// check orthogonality for j_0 .. j_100 on 0..R
void check_orthogonality(int l_max, double R)
{
  int l, i, j, max;
  FILE *file;
  char fp[100];

  for (l = 0; l <= l_max; l++)
  {
    snprintf(fp, sizeof(fp), "results/ortho/%i.txt", l);
    printf("Opening file %s...\n", fp);
    file = fopen(fp, "w+");

    // load wavenumber file
    open_wave_number_file(l);

    max = get_numer_of_zeroes(l);

    // write A_ij to file where A = | M - 1 |
    for (i = 0; i < max; i++)
    {
      for (j = 0; j < max; j++)
      {
        if (i == j)
          fprintf(file, "%f ", fabs(integrate(i, j, l, R) - 1) * pow(10,4));
        else
          fprintf(file, "%f ", fabs(integrate(i, j, l, R)) * pow(10,4));
      }
      fprintf(file, "\n");
    }
    fclose(file);
    close_wave_number_file();
  }
}
