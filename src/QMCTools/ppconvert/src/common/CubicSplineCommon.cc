//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "CubicSplineCommon.h"
#include <iostream>

void CubicSplineCommon::Update()
{
  Grid &x = *grid;
  int N = x.NumPoints;

  Array<double,1> U(N);

  if (StartDeriv > 0.99e30)
    {
      d2y(0) = 0.0;  // Use "Natural" boundary conditions--ie d^2y/dr^2 = 0 
      U(0) = 0.0;
    }
  else
    {
      d2y(0) = -0.5;
      U(0) = (3.0/(x(1)-x(0)))*((y(1)-y(0))/(x(1)-x(0))-StartDeriv);
    }

  d2y(x.NumPoints-1) = 0.0;
  
  // This part translated from Numerical Recipes.
  for (int i=1; i<x.NumPoints-1; i++)
    {
      double sig = (x(i) - x(i-1)) / (x(i+1)-x(i-1));
      double p = sig *d2y(i-1)+2.0;
      d2y(i) = (sig-1.0)/p;
      U(i) = (6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-
		   (y(i)-y(i-1))/(x(i)-x(i-1))) /
	      (x(i+1)-x(i-1))-sig*U(i-1))/p;   
    }

  double Un, Qn;
  if (EndDeriv > 0.99e30)
    {
      Un = 0.0;
      Qn = 0.0;
    }
  else
    {
      Qn = 0.5;
      Un = (3.0/(x(N-1)-x(N-2)))*(EndDeriv-(y(N-1)-y(N-2))/(x(N-1)-x(N-2)));
    }
  
  d2y(N-1) =
    (Un-Qn*U(N-2))/(Qn*d2y(N-2)+1.0);
  
  for (int k=N-2; k>=0; k--)
    d2y(k) = d2y(k)*d2y(k+1) + U(k);

  UpToDate = 1;
}






void
MultiCubicSplineCommon::Update(int i)
{
  Grid &x = *grid;
  int N = x.NumPoints;

  Array<double,1> U(N);

  if (StartDeriv(i) > 0.99e30)
    {
      d2y(0,i) = 0.0;  // Use "Natural" boundary conditions--ie d^2y/dr^2 = 0 
      U(0) = 0.0;
    }
  else
    {
      d2y(0,i) = -0.5;
      U(0) = (3.0/(x(1)-x(0)))*((y(1,i)-y(0,i))/(x(1)-x(0))-StartDeriv(i));
    }

  d2y(NumGridPoints-1,i) = 0.0;
  
  // This part translated from Numerical Recipes.
  for (int j=1; j<NumGridPoints-1; j++)
    {
      double sig = (x(j) - x(j-1)) / (x(j+1)-x(j-1));
      double p = sig *d2y(j-1,i)+2.0;
      d2y(j,i) = (sig-1.0)/p;
      U(j) = (6.0*((y(j+1,i)-y(j,i))/(x(j+1)-x(j))-
		   (y(j,i)-y(j-1,i))/(x(j)-x(j-1))) /
	      (x(j+1)-x(j-1))-sig*U(j-1))/p;   
    }

  double Un, Qn;
  if (EndDeriv(i) > 0.99e30)
    {
      Un = 0.0;
      Qn = 0.0;
    }
  else
    {
      Qn = 0.5;
      Un = (3.0/(x(N-1)-x(N-2)))*(EndDeriv(i)-(y(N-1,i)-y(N-2,i))/(x(N-1)-x(N-2)));
    }
  
  d2y(N-1,i) =
    (Un-Qn*U(N-2))/(Qn*d2y(N-2,i)+1.0);
  
  for (int k=N-2; k>=0; k--)
    d2y(k,i) = d2y(k,i)*d2y(k+1,i) + U(k);

  UpToDate(i) = true;
}



void Test2()
{
  LinearGrid x(0.0, 5.0, 6);
  Array<double,1> y(6);
  
  y = 0.0, 2.0, -2.0, 3.0, -3.0, 1.0;
  CubicSplineCommon spline(&x, y, 1.0, -1.0);
  for (int i=0; i<=5000; i++)
    {
      double r = ((double)i/5000.01) * 5.0;
      fprintf (stderr, "%15.12f %15.12f %15.12f %15.12f %15.12f\n",
	       r, spline(r), spline.Deriv(r), spline.Deriv2(r),
	       spline.Deriv3(r));  
    }
  
}


void TestSpline()
{
  LinearGrid x(0.0, 20.0*M_PI, 501);
  Array<double,1> SinVals(501);

  for (int i=0; i<501; i++)
    SinVals(i) = sin(x(i));
  
  CubicSplineCommon spline (&x, SinVals);

  //  for (int i=0; i<500; i++)
  //  fprintf (stderr, "%15.12f %15.12f\n", spline.y(i), spline.d2y(i));

  for (int i=0; i<=5000; i++)
    {
      double r = ((double)i/5000.0)*2.0*M_PI;
      fprintf (stderr, "%15.16f %1.16f %1.16f %15.16f %15.16f %15.16f %15.16f\n",	       r, sin(r), cos(r), spline(r), spline.Deriv(r), spline.Deriv2(r),
	       spline.Deriv3(r));  
    }
  
}


void TestMultiCubicSplineCommon()
{
  LinearGrid x(0.0, 2.0*M_PI, 51);
  Array<double,2> SinCosVals(51,2);

  for (int i=0; i<51; i++)
    {
      SinCosVals(i,0) = sin(x(i));
      SinCosVals(i,1) = cos(x(i));
    }
  
  Array<double,1> StartDeriv(2);
  Array<double,1> EndDeriv(2);
  StartDeriv(0) = 1.0;
  StartDeriv(1) = 0.0;
  EndDeriv(0) = 1.0;
  EndDeriv(1) = 0.0;
  MultiCubicSplineCommon spline (&x, SinCosVals, StartDeriv, EndDeriv);

  //  for (int i=0; i<500; i++)
  //  fprintf (stderr, "%15.12f %15.12f\n", spline.y(i), spline.d2y(i));

  Array<double,1> yVal(2), yDeriv(2), yDeriv2(2), yDeriv3(2);
  for (int i=0; i<=5000; i++)
    {
      double r = ((double)i/5000.0)*2.0*M_PI;
      spline(r, yVal);
      spline.Deriv (r, yDeriv);
      spline.Deriv2 (r, yDeriv2);
      spline.Deriv3 (r, yDeriv3);

      fprintf (stderr,"%15.16f %1.16f %1.16f %15.16f %15.16f %15.16f %15.16f ",
	       r, sin(r), cos(r), yVal(0), yDeriv(0), yDeriv2(0), yDeriv3(0));
      fprintf (stderr,"%15.16f %15.16f %15.16f %15.16f\n",
	       yVal(1), yDeriv(1), yDeriv2(1), yDeriv3(1));
    }
  
}


CubicSplineCommon&
CubicSplineCommon::operator=(const CubicSplineCommon &spline)
{
  UpToDate = spline.UpToDate;
  y.resize(spline.y.shape());     y   = spline.y;
  d2y.resize(spline.d2y.shape()); d2y = spline.d2y;
  NumParams = spline.NumParams;
  grid = spline.grid;
  StartDeriv = spline.StartDeriv;
  EndDeriv   = spline.EndDeriv;
  
  return *this;
}


//main()
//{
  //Test2();
  //TestSpline();
//  TestMultiCubicSplineCommon();
//}

