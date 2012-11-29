#include <stdlib.h>

#include "random.h"

/*
 * Spin value +/- 1
 *
 */
int get_random_spin()
{
  if (rand() % 2 - 1)
    return 1;
  else
    return -1;
}
