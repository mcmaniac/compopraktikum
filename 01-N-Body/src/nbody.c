#include "nbody.h"

int main ()
{
  int N;
  double t_max, delta_t;
  object* objs = read_data("pla3.txt", &N, &t_max, &delta_t);
  if (!objs)
    return -1;

  printf("N = %i\nt_max = %lf\ndelta_t = %lf\n\n", N, t_max, delta_t);

  int i;
  for (i = 0; i < N; i++)
    print_object(objs[i]);

  free(objs);
  return 0;
}
