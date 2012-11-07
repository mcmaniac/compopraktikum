#include "a2.h"

matrix a2_matrix(int N)
{
  matrix m = null_matrix(N, N);
  int i, j;
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      double v = ((double) N) / ((double) (i+1)+(j+1)) + (i+1) + (j+1);
      MatrixSET(m, i, j, v);
    }
  }
  return m;
}

// Signum function
int sign(double x)
{
  return (x > 0) - (x < 0);
}

void calc_Dprime(matrix D, matrix P, int p, int q)
{
  matrix Dprime = matrix_copy(D);

  double d_pq = MatrixGET(D, p, q),
         d_pp = MatrixGET(D, p, p),
         d_qq = MatrixGET(D, q, q);

  double theta = (d_qq - d_pp) / (2 * d_pq),                             // (11.1.8)
         t     = sign(theta) / (fabs(theta) + sqrt(pow(theta,2.0)+1.0)), // (11.1.10)
         c     = 1.0 / sqrt(pow(t,2.0)+1.0),                             // (11.1.11)
         s     = t * c,                                                  // (11.1.12)
         tau   = s / (1.0 + c);                                          // (11.1.18)

  // set to 0 (11.1.13)
  MatrixSET(Dprime, p, q, 0);
  MatrixSET(Dprime, q, p, 0);

  MatrixSET(Dprime, p, p, d_pp - t*d_pq); // (11.1.14)
  MatrixSET(Dprime, q, q, d_qq + t*d_pq); // (11.1.15)

  // (11.1.16)
  int r;
  for (r = 0; r < D.N; r++)
  {
    if (r != p && r != q)
    {
      double d_rp = MatrixGET(D, r, p),
             d_rq = MatrixGET(D, r, q);

      double dPrime_rp = d_rp - s * (d_rq + tau * d_rp);

      // D symmetric => also D'
      MatrixSET(Dprime, r, p, dPrime_rp);
      MatrixSET(Dprime, p, r, dPrime_rp);
    }
  }

  // (11.1.17)
  for (r = 0; r < D.N; r++)
  {
    if (r != p && r != q)
    {
      double d_rq = MatrixGET(D, r, q),
             d_rp = MatrixGET(D, r, p);

      double dPrime_rq = d_rq + s * (d_rp - tau * d_rq);

      // symmetric again
      MatrixSET(Dprime, r, q, dPrime_rq);
      MatrixSET(Dprime, q, r, dPrime_rq);
    }
  }

  matrix_copy_to(Dprime, D);

  // Transformtion matrix P

  // pq-rotation matrix with values c & s
  matrix P_pq = unity_matrix(D.N);
  MatrixSET(P_pq, p, p, c);
  MatrixSET(P_pq, p, q, s);
  MatrixSET(P_pq, q, p, -1.0*s);
  MatrixSET(P_pq, q, q, c);

  matrix P_n = matrix_mult(P, P_pq);
  matrix_copy_to(P_n, P);
}

void get_biggest_off_diag_index(const matrix A, double *max, int *p, int *q)
{
  *max = MatrixGET(A, 0, 1);
  *p   = 0;
  *q   = 1;
  int i, j;
  for (i = 0; i < A.N; i++)
  {
    for (j = i+1; j < A.M; j++)
    {
      double aij = MatrixGET(A, i, j);
      if (fabs(aij) > fabs(*max))
      {
        *p   = i;
        *q   = j;
        *max = aij;
      }
    }
  }
}

matrix jacobi_diag(const matrix A, matrix P, double epsilon)
{
  // copy A into D
  matrix D = matrix_copy(A);

  // set P to 1-matrix, just to make sure
  matrix_copy_to(unity_matrix(A.N), P);

  // current maximum + index
  double max;
  int p, q;
  get_biggest_off_diag_index(D, &max, &p, &q);

  while (fabs(max) > epsilon)
  {
    calc_Dprime(D, P, p, q);
    get_biggest_off_diag_index(D, &max, &p, &q);
  }

  return D;
}
