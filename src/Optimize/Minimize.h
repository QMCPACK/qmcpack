#ifndef MINIMIZE_H
#define MINIMIZE_H

//#include "../Blitz.h"
#include <blitz/array.h>
class MinimizeFunction
{


public:
  typedef double scalar; 
  scalar dummy;  // for satisfying the compiler.
  virtual int NumParams() = 0;

  virtual scalar &Params(int i) = 0;

  virtual scalar Params(int i) const = 0;

  virtual scalar Cost() = 0;

  virtual void WriteStuff() = 0;

  void ReadParameters(char *FileName);
  void WriteParameters(char *FileName);
};



class Minimizer
{
public:
  virtual void Minimize(MinimizeFunction &MinFunc) = 0;
};


class ConjugateGradient : public Minimizer
{
public:
  typedef double scalar; 
  scalar epsilon;
  scalar Tolerance;
  scalar StepSize;
  MinimizeFunction *MinFunc;
  void Minimize (MinimizeFunction &MinimFunc);
};


class AnnealingSchedule
{
public:
  typedef double scalar; 
  scalar StartTemp, EndTemp;
  int NumTemps;
  int StepsPerTemp;
  virtual scalar Temp(int TempNum) = 0;
};

class ExponentialSchedule : public AnnealingSchedule
{
public:
  typedef double scalar; 
  scalar Chi;
  scalar Temp(int TempNum);
};


class LinearSchedule : public AnnealingSchedule
{
public:
  typedef double scalar; 
  scalar Temp(int TempNum);
};


class VanderbiltAnnealer : public Minimizer
{
public:
  typedef double scalar; 
  AnnealingSchedule *Schedule;
  // The Q matrix gives the step for a random uniform vector, u.
  // \Delta x = Q*u

  //  std::vector<scalar> Q(2);
  blitz::Array<scalar,2> Q;

  // xMean holds the mean position for the last M steps;
  blitz::Array<scalar,1> xMean;
  // Covariance holds the covariance matrix for the last block of steps.
  blitz::Array<scalar,2> Covariance;
  // Holds the Minimum vector found
  blitz::Array<scalar,1> MinParams;
  // Holds the minimum cost found
  scalar MinCost;

  scalar GrowthFactor;
  scalar AcceptRatio;
  scalar MeanCost;
  blitz::Array<scalar,2> s;
  scalar kT;
  scalar HeatCapacity;
  MinimizeFunction *MinFunc;

  void Metropolis(int NumSteps);

  void Minimize (MinimizeFunction &MinimFunc);
  VanderbiltAnnealer(AnnealingSchedule &ASchedule)
  {    
    Schedule = &ASchedule;
  }
};



#endif
