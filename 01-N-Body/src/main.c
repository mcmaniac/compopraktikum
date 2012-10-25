#include "main.h"

int main()
{
  data* dat = read_data(input);
  if (!dat)
  {
    printf("Error reading input data from file \"%s\"\n", input);
    return -1;
  }

  printf("N = %i\nt_max = %f\neta = %f\n\n", dat->N, dat->t_max, dat->eta);

  int i;
  for (i = 0; i < dat->N; i++)
    print_object(dat->objects[i]);

  printf("\nE   = %f\n|j| = %f\n|e| = %f\na_e = %f\n\n",
         total_energy(dat), total_angular_momentum(dat),
         total_runge_lenz(dat), semimajor_axis(dat));

  // start integration method
  printf("Starting integrator...\n");
  file = fopen(output, "w+");
  if (file)
  {
    current_integrator(dat, current_output);
    printf("DONE. Results stored in \"%s\"\n", output);
    fclose(file);
  }
  else
  {
    printf("Error opening output file \"%s\"\n", output);
    return -1;
  }

  free_data(dat);
  return 0;
}
