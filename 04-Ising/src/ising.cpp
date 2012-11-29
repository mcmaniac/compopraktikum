#include <stdlib.h>

#include "ising.h"
#include "random.h"

// for internal use only
#define Val(i,j) (spins[M*((i) % N) + ((j) % M)])

IsingModell::IsingModell(int N_, int M_, double J_)
{
  N = N_;
  M = M_;
  J = J_;
  spins = (int*) malloc(N*M*sizeof(int));

  // initialize with random values
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < M; j++)
    {
      Val(i,j) = get_random_spin();
    }
  }
}

IsingModell::~IsingModell()
{
  free(spins);
}

int IsingModell::getN()
{
  return N;
}

int IsingModell::getM()
{
  return M;
}

int& IsingModell::operator()(int i, int j)
{
  return Val(i,j);
}

double IsingModell::energy(int i, int j)
{
  double s_ij = Val(i,j),
         sum  = 0;
  sum += Val(i-1,j) * s_ij;
  sum += Val(i+1,j) * s_ij;
  sum += Val(i,j-1) * s_ij;
  sum += Val(i,j+1) * s_ij;
  return - J * sum;
}

void IsingModell::flip(int i, int j)
{
  Val(i,j) *= -1;
}


