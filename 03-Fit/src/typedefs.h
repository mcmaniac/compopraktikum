#pragma once

typedef struct {
  double val;
  double delta;
} data;

typedef double (*base_funct)(int l, double x);
