//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "CubicSpline.h"

void 
CubSpline::Update()
{
  SimpleGrid &x = grid;
  int N = x.NumPoints();

  std::vector<double> U(N);

  if (StartDeriv > 0.99e30) {
    d2y[0] = 0.0;  // Use "Natural" boundary conditions--ie d^2y/dr^2 = 0 
    U[0] = 0.0;
  }
  else {
    d2y[0] = -0.5;
    U[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-StartDeriv);
  }

  d2y[x.NumPoints()-1] = 0.0;
  
  // This part translated from Numerical Recipes.
  for (int i=1; i<x.NumPoints()-1; i++) {
    double sig = (x[i] - x[i-1]) / (x[i+1]-x[i-1]);
    double p = sig *d2y[i-1]+2.0;
    d2y[i] = (sig-1.0)/p;
    U[i] = (6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-
		 (y[i]-y[i-1])/(x[i]-x[i-1])) /
	    (x[i+1]-x[i-1])-sig*U[i-1])/p;   
  }

  double Un, Qn;
  if (EndDeriv > 0.99e30) {
    Un = 0.0;
    Qn = 0.0;
  }
  else {
    Qn = 0.5;
    Un = (3.0/(x[N-1]-x[N-2]))*(EndDeriv-(y[N-1]-y[N-2])/(x[N-1]-x[N-2]));
  }
  
  d2y[N-1] = (Un-Qn*U[N-2])/(Qn*d2y[N-2]+1.0);
  
  for (int k=N-2; k>=0; k--)
    d2y[k] = d2y[k]*d2y[k+1] + U[k];

  UpToDate = true;
}
