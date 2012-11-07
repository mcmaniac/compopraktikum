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
      fprintf(file, "%f ", MatrixGET(A, i, j));
    fprintf(file, "\n");
  }
}
