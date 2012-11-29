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

void print_ising(const char *fp, IsingModell &IM, const double T)
{

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

void a1(int mc_steps)
{
  for (double T = 0.0; T <= 6.0; T += 0.1)
  {
    // initialize new ising modell representation
    IsingModell IM(N,M,J);

    for (int n = 0; n < mc_steps; n++)
      monte_carlo_step(IM, k_B, T);

    char fp[100];
    snprintf(fp, sizeof(fp), "results/a1/T-%.1f.txt", T);
    print_ising(fp, IM, T);
  }
}

void a2(int mc_steps)
{
  // use same ising modell
  IsingModell IM(N,M,J);
  for (double T = 6.0; T > 0; T -= 0.1)
  {
    for (int n = 0; n < mc_steps; n++)
      monte_carlo_step(IM, k_B, T);

    char fp[100];
    snprintf(fp, sizeof(fp), "results/a2/T-%.1f.txt", T);
    print_ising(fp, IM, T);
  }
}

int main()
{
  // init random seed
  srand(time(NULL));

  a1(100);
  a2(100);
  return 0;
}
