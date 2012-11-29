#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <iostream>
#include <fstream>

#include "ising.h"
#include "random.h"
#include "monte_carlo.h"

using namespace std;

const double J = 1;
const double k_B = 1;

const int N = 50;
const int M = 50;

// number of monte carlo steps
const int mc_steps = 100;

void print_ising(IsingModell &IM, const double T)
{
  char fp[100];
  snprintf(fp, sizeof(fp), "results/a1/T-%.1f.txt", T);

  // open file
  cout << "Opening " << fp << " ... ";
  ofstream output;
  output.open(fp);

  // write values to file
  for (int i = 0; i < IM.getN(); i++)
  {
    for (int j = 0; j < IM.getM(); j++)
      output << IM(i,j) << " ";
    output << endl;
  }

  output.close();
  cout << "DONE" << endl;
}

int main()
{
  // init random seed
  srand(time(NULL));

  for (double T = 0.1; T <= 6.0; T += 0.1)
  {
    // initialize new ising modell representation
    IsingModell IM(N,M,J);

    for (int n = 0; n < mc_steps; n++)
    {
      monte_carlo_step(IM, k_B, T);
    }

    print_ising(IM, T);
  }
}
