#include "matrix.h"

void matrix_destroy(matrix A)
{
  free(A.val);
}

matrix null_matrix(int N, int M)
{
  int i, j;
  matrix A = {
    .N   = N,
    .M   = M,
    .val = (double*) malloc(N*M*sizeof(double))
  };
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      MatrixSET(A, i, j, 0);
  return A;
}

matrix unity_matrix(int N)
{
  int i, j;
  matrix A = {
    .N   = N,
    .M   = N,
    .val = (double*) malloc(N*N*sizeof(double))
  };
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      if (i == j)
        MatrixSET(A, i, j, 1);
      else
        MatrixSET(A, i, j, 0);
    }
  }
  return A;
}

matrix matrix_copy(const matrix A)
{
  matrix B = null_matrix(A.N, A.M);
  matrix_copy_to(A, B);
  return B;
}

// Replace values of matrix B with values of matrix A
void matrix_copy_to(const matrix A, matrix B)
{
  int i, j;
  for (i = 0; i < A.N; i++)
    for (j = 0; j < A.M; j++)
      MatrixSET(B, i, j, MatrixGET(A, i, j));
}

matrix matrix_add(const matrix A, const matrix B)
{
  matrix C = null_matrix(A.N, A.M);
  int i, j;
  for (i = 0; i < A.N; i++)
    for (j = 0; j < A.M; j++)
      MatrixSET(C, i, j, MatrixGET(A, i, j) + MatrixGET(B, i, j));
  return C;
}

matrix matrix_mult(const matrix A, const matrix B)
{
  matrix C = null_matrix(A.N, B.M);
  int i, j, k;
  double sum;
  for (i = 0; i < A.N; i++)
  {
    for (j = 0; j < B.M; j++)
    {
      sum = 0;
      for (k = 0; k < A.M; k++)
      {
        sum += MatrixGET(A,i,k) * MatrixGET(B,k,j);
      }
      MatrixSET(C, i, j, sum);
    }
  }
  return C;
}

void matrix_print(const matrix A)
{
  int i, j;
  for (i = 0; i < A.N; i++)
  {
    for (j = 0; j < A.M; j++)
      printf("%f ", MatrixGET(A, i, j));
    printf("\n");
  }
}

void matrix_fprint(FILE *file, const matrix A)
{
  int i, j;
  for (i = 0; i < A.N; i++)
  {
    for (j = 0; j < A.M; j++)
      fprintf(file, "%e ", MatrixGET(A, i, j));
    fprintf(file, "\n");
  }
}

void matrix_swap_rows(matrix A, int k, int l)
{
  if (k == l)
    return;
  int j;
  for (j = 0; j < A.M; j++)
  {
    double tmp = MatrixGET(A, k, j);
    MatrixSET(A, k, j, MatrixGET(A, l, j));
    MatrixSET(A, l, j, tmp);
  }
}

vector solve_gauss(const matrix A, const vector b_)
{
  if (A.N != A.M)
    printf("Warning: Gauss only for NxN matrices\n");

  matrix L = null_matrix(A.N, A.M),
         R = matrix_copy(A);
  vector b = vector_copy(b_);

  // Find pivots
  int i, j, k;
  for (k = 0; k < A.N; k++)
  {
    double max = 0;
    int i_max = k;
    for (i = k; i < A.N; i++)
    {
      if (fabs(MatrixGET(R,i,k)) > max)
      {
        max   = fabs(MatrixGET(R,i,k));
        i_max = i;
      }
    }
    matrix_swap_rows(R, k, i_max);
    vector_swap(b, k, i_max);
  }

  // LU decomposition
  for (k = 0; k < A.N-1; k++)
  {
    for (i = k+1; i < A.N; i++)
    {
      MatrixSET(L,i,k, MatrixGET(R,i,k) / MatrixGET(R,k,k));
    }
    for (j = k+1; j < A.N; j++)
    {
      for (i = k+1; i < A.N; i++)
      {
        MatrixSET(R,i,j, MatrixGET(R,i,j) - MatrixGET(L,i,k) * MatrixGET(R,k,j));
      }
    }
  }
  // clean up R
  for (i = 1; i < A.N; i++)
  {
    for (j = 0; j < i; j++)
    {
      MatrixSET(R,i,j, 0);
    }
  }
  // set L_ii to "1"
  for (i = 0; i < A.N; i++)
  {
    MatrixSET(L,i,i, 1);
  }

  // forwards
  vector y = null_vector(A.N);
  for (i = 0; i < A.N; i++)
  {
    double sum = 0;
    for (k = 0; k < i; k++)
      sum += MatrixGET(L,i,k) * VectorGET(y,k);
    VectorSET(y,i, 1 / MatrixGET(L,i,i) * (VectorGET(b,i) - sum));
  }

  // backwards
  vector x = null_vector(A.N);
  for (i = A.N-1; i >= 0; i--)
  {
    double sum = 0;
    for (k = i+1; k < A.N; k++)
      sum += MatrixGET(R,i,k) * VectorGET(x,k);
    VectorSET(x,i, 1 / MatrixGET(R,i,i) * (VectorGET(y,i) - sum));
  }

  matrix_destroy(R);
  matrix_destroy(L);
  vector_destroy(y);

  return x;
}
