#pragma once
#include <stdio.h>

#define PI 3.14159265358979

typedef struct {
  int l;
  double x;
  double val;
} bessel;

typedef struct {
  int l;
  FILE *file;
} wavenumber_file;
