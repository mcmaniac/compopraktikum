#include "bessel.h"

bessel j_plus(const bessel j, const bessel j_m)
{
  bessel j_p = {
    .l   = j.l + 1,
    .x   = j.x,
    .val = (2.0 * j.l + 1.0) / j.x * j.val - j_m.val
  };
  return j_p;
}

bessel j_minus(const bessel j, const bessel j_p)
{
  bessel j_m = {
    .l   = j.l - 1,
    .x   = j.x,
    .val = (2.0 * j.l + 1) / j.x * j.val - j_p.val
  };
  return j_m;
}

bessel j_0(double x)
{
  bessel j = {
    .l   = 0,
    .x   = x,
    .val = sin(x) / x
  };
  return j;
}

bessel j_1(double x)
{
  bessel j = {
    .l   = 1,
    .x   = x,
    .val = sin(x) / (x*x) - cos(x) / x
  };
  return j;
}

// j_M(x)
bessel j_M(int M, double x)
{
  bessel j = {
    .l   = M,
    .x   = x,
    .val = 0
  };
  return j;
}

// j_(M-1)(x)
bessel j_M_m(int M, double x)
{
  bessel j = {
    .l   = M-1,
    .x   = x,
    .val = 1
  };
  return j;
}

bessel j_l_upwards(int l, double x)
{
  int i;
  bessel j_m = j_0(x),
         j   = j_1(x),
         j_p;
  for (i = 1; i < l; i++)
  {
    // new value
    j_p = j_plus(j, j_m);
    // shift
    j_m = j;
    j   = j_p;
  }
  return j;
}

bessel j_l_downwards(int l, int M, double x)
{
  int i;
  bessel j_p = j_M(M, x),
         j   = j_M_m(M, x),
         j_m,
         j_l,
         j0  = j_0(x);
  for (i = M-2; i >= 0; i--)
  {
    j_m = j_minus(j, j_p);
    // store j_l
    if (i == l)
      j_l = j_m;
    // shift
    j_p = j;
    j   = j_m;
  }
  j_l.val = j_l.val * j0.val / j.val; // " j_l / c "
  return j_l;
}

bessel j_l(int l, double x)
{
  if (l == 0)
    return j_0(x);
  else if (l == 1)
    return j_1(x);
  else if (x > l)
    return j_l_upwards(l, x);
  else
    return j_l_downwards(l, l*l, x);
}

// Calculate j'_l(x) with j_l and j_(l-1)
bessel deriv_up(const bessel j, const bessel j_m)
{
  bessel j_d = {
    .l   = j.l,
    .x   = j.x,
    .val = j_m.val - (j.l + 1) / j.x * j.val
  };
  return j_d;
}

// Calculate j'_l(x) with j_l and j_(l+1)
bessel deriv_down(const bessel j, const bessel j_p)
{
  bessel j_d = {
    .l   = j.l,
    .x   = j.x,
    .val = j.l / j.x * j.val - j_p.val
  };
  return j_d;
}

bessel derive(const bessel j)
{
  bessel j_d = deriv_down(j, j_l(j.l + 1, j.x));
  return j_d;
}
