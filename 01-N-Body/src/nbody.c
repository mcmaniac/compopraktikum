#include "nbody.h"

// gravitation constant
double G = 1;

//const char* input  = "in2.txt";
//const char* input  = "pla3.txt";
//const char* input  = "pl.100";
const char* input  = "pl.1k";

//const char* output = "runge_kutta_2.txt";
//const char* output = "runge_kutta_3.txt";
//const char* output = "runge_kutta_100.txt";
const char* output = "runge_kutta_1k.txt";

FILE* file;

int main()
{
  data* dat = read_data(input);

  /*
  if (!dat)
  {
    printf("Error reading input data from file \"%s\"\n", input);
    return -1;
  }

  printf("N = %i\nt_max = %lf\neta = %lf\n\n", dat->N, dat->t_max, dat->eta);


  int i;
  for (i = 0; i < dat->N; i++)
    print_object(dat->objects[i]);
  */

  double a_e = semimajor_axis(dat);
  printf("a_e = %lf\n", a_e);

  /*
  printf("\nE   = %lf\n|j| = %lf\n|e| = %lf\na_e = %lf\n\n",
         total_energy(dat), total_angular_momentum(dat),
         total_runge_lenz(dat), semimajor_axis(dat));

  // open file
  file = fopen(output, "w+");
  if (file)
  {
    //printf("Starting runge kutta...\n");
    //runge_kutta(dat, &print_constants_to_file);
    //printf("DONE.\n");
    free_data(dat);
    fclose(file);
    return 0;
  }
  else
  {
    printf("Error opening output file \"%s\"\n", output);
    return -1;
  }
  */

  return 0;
}
