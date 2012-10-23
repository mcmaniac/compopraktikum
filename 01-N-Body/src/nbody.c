#include "nbody.h"

// gravitation constant
double G = 1;

int main()
{
  data* dat = read_data("in2.txt");
  //data* dat = read_data("pla3.txt");
  //data* dat = read_data("pl.100");
  //data* dat = read_data("pl.1k");

  if (!dat)
    return -1;

  printf("N = %i\nt_max = %lf\neta = %lf\n\n", dat->N, dat->t_max, dat->eta);

  int i;
  for (i = 0; i < dat->N; i++)
    print_object(dat->objects[i]);

  printf("\nCleaning runge kutta file…\n");
  FILE *file = fopen("runge_kutta_2.txt", "w");
  //FILE *file = fopen("runge_kutta_3.txt", "w");
  //FILE *file = fopen("runge_kutta_100.txt", "w");
  //FILE *file = fopen("runge_kutta_1k.txt", "w");
  fprintf(file, "");

  printf("Starting runge kutta…\n");
  runge_kutta(dat, &print_runge_kutta);

  printf("DONE.\n");

  free_data(dat);
  return 0;
}

void print_runge_kutta(double time, double delta_t, const data* dat)
{
  print_constants_to_file("runge_kutta_2.txt", time, delta_t, dat);
  //print_constants_to_file("runge_kutta_3.txt", time, delta_t, dat);
  //print_constants_to_file("runge_kutta_100.txt", time, delta_t, dat);
  //print_constants_to_file("runge_kutta_1k.txt", time, delta_t, dat);
}
