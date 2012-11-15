#include "a1.h"

void a1(void)
{
  int N = 3,
      M = 1;
  data *dat = (data*) malloc(N*sizeof(data));
  dat[0] = (data) { .val = 299793.0, .err = 2.0 };
  dat[1] = (data) { .val = 299792.0, .err = 4.5 };
  dat[2] = (data) { .val = 299782.0, .err = 25.0 };

  double avg   = average_value(N, dat),
         sig_i = sigma_square_intern(N, dat),
         sig_e = sigma_square_extern(M, N, dat);

  printf("c_avg   = %f\n", avg);
  printf("sigma_i = %f\n", sqrt(sig_i));
  printf("sigma_e = %f\n", sqrt(sig_e));
}

int main()
{
  a1();
  return 0;
}
