#pragma once

#include <stdlib.h>
#include <stdio.h>

typedef struct {
  int N;
  int M;
  double *val;
} matrix;

#define MatrixGET(A, i, j)    ((A).val[(A).M*(i)+(j)])
#define MatrixSET(A, i, j, v) ((A).val[(A).M*(i)+(j)] = (v))

matrix null_matrix(int N, int M);
matrix unity_matrix(int N);
void   matrix_destroy(matrix A);

matrix matrix_copy(const matrix A);
void   matrix_copy_to(const matrix A, matrix B);

matrix matrix_add(const matrix A, const matrix B);
matrix matrix_mult(const matrix A, const matrix B);

void matrix_print(const matrix A);
