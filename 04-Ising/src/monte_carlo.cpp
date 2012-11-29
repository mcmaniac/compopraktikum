#include <stdlib.h>
#include <math.h>

#include <iostream>

#include "ising.h"
#include "monte_carlo.h"

using namespace std;

void monte_carlo_step(IsingModell &IM, const double k_B, const double T)
{
  int N = IM.getN(),
      M = IM.getM();
  for (int n = 0; n < N*M; n++)
  {
    // pick random coordinate
    int i = rand() % N,
        j = rand() % M;

    // get current energy at (i,j)
    double E_cur = IM.energy(i,j);

    // get energy at (i,j) with (i,j) flipped
    IM.flip(i,j);
    double E_flipped = IM.energy(i,j);

    double delta_E = E_cur - E_flipped;

    // check whether we have to revert the flip
    if (delta_E > 0)
    {
      // pick a random number between 0 and 1
      double R = rand() % 1000 / 1000.0;
      // calculate P
      double P = exp(- delta_E / (k_B * T));
      // revert the flip if not R < P
      if (R >= P)
        IM.flip(i,j);
    }
  }
}
