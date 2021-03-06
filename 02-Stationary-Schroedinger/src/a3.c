#include "a3.h"

double potential(double r, double r_0, double a_0, double V_0)
{
  return V_0 / (1 + exp((r-r_0)/a_0));
}

double a3_pot(double r)
{
  double V_0 = -60.0 * pow(10,6),
         r_0 = 3.0   * pow(10,-15),
         a_0 = 0.4   * pow(10,-15);
  return potential(r, r_0, a_0, V_0);
}

matrix calc_T(int l, double M, double R)
{
  double hbarc = 197.3 * pow(10,-9);

  int i,
      N = get_numer_of_zeroes(l);
  matrix T = null_matrix(N, N);
  for (i = 0; i < N; i++)
    MatrixSET(T, i, i, pow(hbarc,2) / (2 * M) * pow(get_kjl(i,l,R),2));

  return T;
}

// Basically same integration function as in "a1d.c", but with a V(r) factor
double integrate_V(int l, int i, int j, double R)
{
  // wave numbers & norm factors
  double kil = get_kjl(i, l, R),
         kjl = get_kjl(j, l, R);

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
    res += a3_pot(r) * pow(r,2) * j_l(l, kil * r).val * j_l(l, kjl * r).val;
  }

  // r = R value
  res += 0.5 * a3_pot(R) * pow(R,2) * j_l(l, kil * R).val * j_l(l, kjl * R).val;

  // norm values
  double ail = norm_factor(kil, i, l, R),
         ajl = norm_factor(kjl, j, l, R);

  return res * dr * ail * ajl; // multiply with norm factors
}

matrix calc_V(int l, double R)
{
  int N = get_numer_of_zeroes(l);
  matrix V = null_matrix(N, N);
  int i, j;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      MatrixSET(V, i, j, integrate_V(l, i, j, R));
  return V;
}

matrix calc_H(int l, double R, double M)
{
  matrix H = matrix_add(calc_V(l, R), calc_T(l, M, R));
  return H;
}

double calc_psi(int j, int l, double r, double R, const matrix P)
{
  int i;
  double kil, ail,
         psi = 0;

  for (i = 0; i < P.N; i++)
  {
    kil = get_kjl(i, l, R);
    ail = norm_factor(kil, i, l, R);
    psi += MatrixGET(P, i, j) * ail * j_l(l, kil * r).val;
  }

  return psi;
}
