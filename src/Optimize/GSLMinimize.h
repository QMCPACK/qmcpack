//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef GSL_MINIMIZE_H
#define GSL_MINIMIZE_H

#include <vector>
#include "Optimize/OptimizeBase.h"
#include "OhmmsData/ParameterSet.h"
/** minimizer using GSL library
 *
 * @authors K. Esler, J. Vincent and J. Kim
 *
 * This class is only available for internal use by the developers
 * and cannot be distributed with other codes.
 */
class GSLConjugateGradient : public MinimizerBase<double>
{
public:
  typedef double Return_t;
  int MaxCGStep;
  Return_t epsilon;
  Return_t Tolerance;
  Return_t StepSize;
  ObjectFuncType *MinFunc;
  GSLConjugateGradient();
  bool optimize(ObjectFuncType *atarget);
  bool put(xmlNodePtr cur);
};


//class AnnealingSchedule
//{
//public:
//  typedef double scalar;
//  scalar StartTemp, EndTemp;
//  int NumTemps;
//  int StepsPerTemp;
//  virtual scalar Temp(int TempNum) = 0;
//};
//
//class ExponentialSchedule : public AnnealingSchedule
//{
//public:
//  typedef double scalar;
//  scalar Chi;
//  scalar Temp(int TempNum);
//};
//
//
//class LinearSchedule : public AnnealingSchedule
//{
//public:
//  typedef double scalar;
//  scalar Temp(int TempNum);
//};


//class VanderbiltAnnealer : public Minimizer
//{
//public:
//  typedef double scalar;
//  AnnealingSchedule *Schedule;
//  // The Q matrix gives the step for a random uniform vector, u.
//  // \Delta x = Q*u
//
//  //  std::vector<scalar> Q(2);
//  std::vector<scalar> Q;
//
//  // xMean holds the mean position for the last M steps;
//  //blitz::Array<scalar,1> xMean;
//  std::vector<scalar> xMean;
//  // Covariance holds the covariance matrix for the last block of steps.
//  //blitz::Array<scalar,2> Covariance;
//  std::vector<scalar> Covariance;
//  // Holds the Minimum vector found
//  std::vector<scalar> MinParams;
//  // Holds the minimum cost found
//  scalar MinCost;
//
//  scalar GrowthFactor;
//  scalar AcceptRatio;
//  scalar MeanCost;
//  //blitz::Array<scalar,2> s;
//  std::vector<scalar> s;
//  scalar kT;
//  scalar HeatCapacity;
//  MinimizeFunction *MinFunc;
//
//  void Metropolis(int NumSteps);
//
//  void Minimize (MinimizeFunction &MinimFunc);
//  VanderbiltAnnealer(AnnealingSchedule &ASchedule)
//  {
//    Schedule = &ASchedule;
//  }
//};
//
#endif
