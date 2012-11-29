#pragma once

class IsingModell
{
  private:
    int N, M; // dimensions
    int* spins;
    double J;

  public:
    IsingModell(int, int, double);
    ~IsingModell();

    int getN();
    int getM();

    int& operator()(int,int);

    double energy(int i, int j);
    void flip(int i, int j);
};
