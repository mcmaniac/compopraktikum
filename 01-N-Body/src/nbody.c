#include "nbody.h"

// gravitation constant
double G = 1;

int main()
{
  data* dat = read_data("in2.txt");
  if (!dat)
    return -1;

  printf("N = %i\nt_max = %lf\neta = %lf\n\n", dat->N, dat->t_max, dat->eta);

  int i;
  for (i = 0; i < dat->N; i++)
    print_object(dat->objects[i]);

  free_data(dat);
  return 0;
}
